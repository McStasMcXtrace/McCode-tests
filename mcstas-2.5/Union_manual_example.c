/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: Union_manual_example.instr (Union_manual_example)
 * Date:       Fri Oct 11 20:27:06 2019
 * File:       Union_manual_example.c
 * Compile:    cc -o Union_manual_example.out Union_manual_example.c  -I@MCCODE_LIB@/share/
 * CFLAGS= -I@MCCODE_LIB@/share/
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

#line 695 "Union_manual_example.c"

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

#line 928 "Union_manual_example.c"

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

#line 4947 "Union_manual_example.c"

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

#line 5307 "Union_manual_example.c"
#ifdef MC_TRACE_ENABLED
int mctraceenabled = 1;
#else
int mctraceenabled = 0;
#endif
#define MCSTAS "/usr/share/mcstas/2.5/"
int mcdefaultmain = 1;
char mcinstrument_name[] = "Union_manual_example";
char mcinstrument_source[] = "Union_manual_example.instr";
char *mcinstrument_exe=NULL; /* will be set to argv[0] in main */
int main(int argc, char *argv[]){return mccode_main(argc, argv);}
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);

/* Shared user declarations for all components 'Incoherent_process'. */
#line 59 "/usr/share/mcstas/2.5/contrib/union/Incoherent_process.comp"
#ifndef Union
#define Union $Revision: 0.8 $

#include "Union_functions.c"
#include "Union_initialization.c"

#endif


struct Incoherent_physics_storage_struct{
// Variables that needs to be transfered between any of the following places:
// The initialize in this component
// The function for calculating my
// The function for calculating scattering

double my_scattering;
};

// Function for calculating my in Incoherent case
int Incoherent_physics_my(double *my,double *k_initial, union data_transfer_union data_transfer, struct focus_data_struct *focus_data) {
    *my = data_transfer.pointer_to_a_Incoherent_physics_storage_struct->my_scattering;
    return 1;
};

// Function for basic incoherent scattering event
int Incoherent_physics_scattering(double *k_final, double *k_initial, double *weight, union data_transfer_union data_transfer, struct focus_data_struct *focus_data) {

    //New version of incoherent scattering
    double k_length = sqrt(k_initial[0]*k_initial[0]+k_initial[1]*k_initial[1]+k_initial[2]*k_initial[2]);

    Coords k_out;
    // Here is the focusing system in action, get a vector
    double solid_angle;
    focus_data->focusing_function(&k_out,&solid_angle,focus_data);
    NORM(k_out.x,k_out.y,k_out.z);
    *weight *= solid_angle*0.25/PI;

    k_final[0] = k_out.x*k_length; k_final[1] = k_out.y*k_length; k_final[2] = k_out.z*k_length;
    return 1;
};

#line 5368 "Union_manual_example.c"

/* Shared user declarations for all components 'Powder_process'. */
#line 56 "/usr/share/mcstas/2.5/contrib/union/Powder_process.comp"
#ifndef Union
#define Union $Revision: 0.8 $

#include "Union_functions.c"
#include "Union_initialization.c"

#endif



// Share section of PowderN 8/3 2016 from McStas.org

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

/* Declare structures and functions only once in each instrument. */
#ifndef POWDERN_DECL_UNION
#define POWDERN_DECL_UNION
/* format definitions in the order {j d F2 DW Dd inv2d q F strain} */
#ifndef Crystallographica
#define Crystallographica { 4,5,7,0,0,0,0,0,0 }
#define Fullprof          { 4,0,8,0,0,5,0,0,0 }
#define Lazy              {17,6,0,0,0,0,0,13,0 }
#define Undefined         { 0,0,0,0,0,0,0,0,0 }
#endif

struct line_data_union
{
double F2;                  /* Value of structure factor */
double q;                   /* Qvector */
int j;                      /* Multiplicity */
double DWfactor;            /* Debye-Waller factor */
double w;                   /* Intrinsic line width */
double Epsilon;             /* Strain=delta_d_d/d shift in ppm */
};

struct line_info_struct_union
{
struct line_data_union *list;     /* Reflection array */
int  count;                  /* Number of reflections */
double Dd;
double DWfactor;
double V_0;
double rho;
double at_weight;
double at_nb;
double sigma_a; // should not be used
double sigma_i; // should not be used
char   compname[256];
double flag_barns;
int    shape; /* 0 cylinder, 1 box, 2 sphere, 3 OFF file */
int    column_order[9]; /* column signification */
int    flag_warning;
char   type;  /* interaction type of event t=Transmit, i=Incoherent, c=Coherent */
double dq;    /* wavevector transfer [Angs-1] */
double Epsilon; /* global strain in ppm */
double XsectionFactor;
double my_s_v2_sum;
double my_a_v;
double my_inc;
double *w_v,*q_v, *my_s_v2;
double radius_i,xwidth_i,yheight_i,zdepth_i; // not to be used, but still here
double v; /* last velocity (cached) */
      double Nq;
      int    nb_reuses, nb_refl, nb_refl_count;
      double v_min, v_max;
      double xs_Nq[CHAR_BUF_LENGTH];
      double xs_sum[CHAR_BUF_LENGTH];
      double neutron_passed;
      long   xs_compute, xs_reuse, xs_calls;
    };

  off_struct offdata_union;
  
  // PN_list_compare *****************************************************************

  int PN_list_compare_union (void const *a, void const *b)
  {
     struct line_data_union const *pa = a;
     struct line_data_union const *pb = b;
     double s = pa->q - pb->q;
     
     if (!s) return 0;
     else    return (s < 0 ? -1 : 1);
  } /* PN_list_compare */

  int read_line_data_union(char *SC_file, struct line_info_struct_union *info)
  {
    struct line_data_union *list = NULL;
    int    size = 0;
    t_Table sTable; /* sample data table structure from SC_file */
    int    i=0;
    int    mult_count  =0;
    char   flag=0;
    double q_count=0, j_count=0, F2_count=0;
    char **parsing;
    int    list_count=0;

    if (!SC_file || !strlen(SC_file) || !strcmp(SC_file, "NULL")) {
      printf("PowderN: %s: Using incoherent elastic scattering only\n",info->compname);
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
    else size = sTable.rows;
    Table_Info(sTable);
    printf("PowderN: %s: Reading %d rows from %s\n",
          info->compname, size, SC_file);

    if (info->column_order[0] == 4 && info->flag_barns !=0)
    printf("PowderN: %s: Powder file probably of type Crystallographica/Fullprof (lau)\n"
           "WARNING: but F2 unit is set to barns=1 (barns). Intensity might be 100 times too high.\n",
           info->compname);
  if (info->column_order[0] == 17 && info->flag_barns == 0)
    printf("PowderN: %s: Powder file probably of type Lazy Pulver (laz)\n"
           "WARNING: but F2 unit is set to barns=0 (fm^2). Intensity might be 100 times too low.\n",
           info->compname);
    /* allocate line_data array */
    list = (struct line_data_union*)malloc(size*sizeof(struct line_data_union));

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
        printf("PowderN: %s: line %i has invalid definition\n"
               "         (mult=0 or q=0 or d=0)\n", info->compname, i);
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
          printf("PowderN: %s: Set multiplicity to 1 for lines [%i:%i]\n"
                  "         (d-spacing %g is duplicated %i times)\n",
            info->compname, list_count-mult_count, list_count-1, list[list_count-1].q, mult_count);
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
    qsort(list, list_count, sizeof(struct line_data_union),  PN_list_compare_union);
    
    printf("PowderN: %s: Read %i reflections from file '%s'\n", 
      info->compname, list_count, SC_file);
    
    info->list  = list;
    info->count = list_count;

    return(list_count);
  } /* read_line_data_union */


/* computes the number of possible reflections (return value), and the total xsection 'sum' */
/* this routine looks for a pre-computed value in the Nq and sum cache tables               */
/* when found, the earch starts from the corresponding lower element in the table           */
int calc_xsect_union(double v, double *qv, double *my_sv2, int count, double *sum,
  struct line_info_struct_union *line_info) {
  int    Nq = 0, line=0, line0=0;
  *sum=0;
  
  //printf("Line_info when entering cross_section calculation\n");
  //printf("v = %f, qv = %f, my_sv2 = %f, count = %d, sum = %f\n",v,*qv,*my_sv2,count,*sum);
  //printf("v = %f\n",v);
  //printf("line_info->v = %f, line_info->v_min = %f, line_info->v_max = %f, line_info->neutron_passed = %f\n",line_info->v,line_info->v_min,line_info->v_max,line_info->neutron_passed);
  //printf("line_info->xs_reuses = %d, line_info->xs_compute = %d\n",line_info->xs_reuse,line_info->xs_compute);
  
  
  /* check if a line_info element has been recorded already */
  if (v >= line_info->v_min && v <= line_info->v_max && line_info->neutron_passed >= CHAR_BUF_LENGTH) {
    line = (int)floor(v - line_info->v_min)*CHAR_BUF_LENGTH/(line_info->v_max - line_info->v_min);
    Nq    = line_info->xs_Nq[line];
    *sum  = line_info->xs_sum[line];
    if (!Nq && *sum == 0) {
      /* not yet set: we compute the sum up to the corresponding speed in the table cache */
      //printf("Nq and sum not yet set, have to do this calculation now\n");
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
      //printf("line_info->xs_Nq[line] = %f, line_info->xs_sum[line] = %f, line_info->xs_compute = %d\n",line_info->xs_Nq[line],line_info->xs_sum[line],line_info->xs_compute);
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
  
  //printf("cross_section function to return %d lines to scatter with, with cross section sum %f \n",Nq,*sum);
  return(Nq);
} /* calc_xsect_union */

#endif /* !POWDERN_DECL */



struct Powder_physics_storage_struct{
    // Variables that needs to be transfered between any of the following places:
    // The initialize in this component
    // The function for calculating my
    // The function for calculating scattering
    
    struct line_info_struct_union *line_info_storage;
    double my_scattering;
    double vertical_angular_limit;
};

// Obsolete: Function for initializing test_physics. Done in component instead.
int Powder_physics_initialize(union data_transfer_union data_transfer) {
      // Obsolte
      return 1;
};

// Function for calculating my in a test case.
int Powder_physics_my(double *my,double *k_initial, union data_transfer_union data_transfer, struct focus_data_struct *focus_data) {
    //*my = data_transfer.pointer_to_a_Powder_physics_storage_struct->my_scattering;
    
    
    int method_switch = 1;
    // For test
    int line_v,line0,line,count;
    
    // Should not interfer with the global variables
    double vx = k_initial[0]*K2V;
    double vy = k_initial[1]*K2V;
    double vz = k_initial[2]*K2V;
    
    // Not sure one can do this, but I do not see why not
    struct line_info_struct_union *line_info = data_transfer.pointer_to_a_Powder_physics_storage_struct->line_info_storage;
    
    double v = sqrt(vx*vx + vy*vy + vz*vz);
    //printf("Velocity = %f \n",v);
    
    //printf("line_info->v = %f, line_info->v_min = %f, line_info->v_max = %f, line_info->neutron_passed = %f\n",line_info->v,line_info->v_min,line_info->v_max,line_info->neutron_passed);
    // Here the maximum and minimum v is recorded, should this be for scattering events or cross section calculations?
    if (line_info->neutron_passed < CHAR_BUF_LENGTH) {
      if (v < line_info->v_min) line_info->v_min = v;
      if (v > line_info->v_max) line_info->v_max = v;
      line_info->neutron_passed++;
    }
    
    if (method_switch == 1) {
    // Here the cross section is calculated and stored
    if ( fabs(v - line_info->v) < 1e-6) {
        line_info->nb_reuses++;
      } else {
        //printf("calling crosssection calculation \n");
        // int calc_xsect_union(double v, double *qv, double *my_sv2, int count, double *sum, struct line_info_struct *line_info)
        line_info->Nq = calc_xsect_union(v, line_info->q_v, line_info->my_s_v2, line_info->count, &line_info->my_s_v2_sum, line_info);
        line_info->v = v;
        line_info->nb_refl += line_info->Nq;
        line_info->nb_refl_count++;
      }
    } else {
    if ( fabs(v - line_info->v) < 1e-6) {
        line_info->nb_reuses++;
      } else {
        //printf("calling crosssection calculation \n");
        if (v >= line_info->v_min && v <= line_info->v_max && line_info->neutron_passed >= CHAR_BUF_LENGTH) {
        line = (int)floor(v - line_info->v_min)*CHAR_BUF_LENGTH/(line_info->v_max - line_info->v_min);
        line_info->Nq = line_info->xs_Nq[line];
        line_info->my_s_v2_sum  = line_info->xs_sum[line];
        if (!line_info->Nq && line_info->my_s_v2_sum == 0) {
          /* not yet set: we compute the sum up to the corresponding speed in the table cache */
          //printf("Nq and sum not yet set, have to do this calculation now\n");
          double line_v = line_info->v_min + line*(line_info->v_max - line_info->v_min)/CHAR_BUF_LENGTH;
          for(line0=0; line0<count; line0++) {
            if (line_info->q_v[line0] <= 2*line_v) { /* q < 2*kf: restrict structural range */
              line_info->my_s_v2_sum += line_info->my_s_v2[line0];
              if (line_info->Nq < line0+1) line_info->Nq=line0+1; /* determine maximum line index which can scatter */
            } else break;
          }
          line_info->xs_Nq[line] = line_info->Nq;
          line_info->xs_sum[line]= line_info->my_s_v2_sum;
          line_info->xs_compute++;
          //printf("line_info->xs_Nq[line] = %f, line_info->xs_sum[line] = %f, line_info->xs_compute = %d\n",line_info->xs_Nq[line],line_info->xs_sum[line],line_info->xs_compute);
        } else line_info->xs_reuse++;
        line0 = line_info->Nq;
        }
          
        line_info->xs_calls++;
          
        for(line=line0; line<count; line++) {
          if (line_info->q_v[line] <= 2*v) { /* q < 2*kf: restrict structural range */
            line_info->my_s_v2_sum += line_info->my_s_v2[line];
            if (line_info->Nq < line+1) line_info->Nq=line+1; /* determine maximum line index which can scatter */
          } else break;
        }
        line_info->v = v;
        line_info->nb_refl += line_info->Nq;
        line_info->nb_refl_count++;
      }
    }
    
     *my = line_info->my_s_v2_sum/(v*v);
    //printf("Returned my scattering of %f \n",*my);
    //printf("compute = %d and reuse = %d \n",line_info->xs_compute,line_info->xs_reuse);
    
    return 1;
};

// Function that provides a basic nonuniform elastic scattering. Unphysical for testing purposes.
int Powder_physics_scattering(double *k_final, double *k_initial, double *weight, union data_transfer_union data_transfer, struct focus_data_struct *focus_data) {

    // This component need to write to its storage transfer for each event, is that possible with this structure?
    struct line_info_struct_union *line_info = data_transfer.pointer_to_a_Powder_physics_storage_struct->line_info_storage;
    double vertical_angular_limit = data_transfer.pointer_to_a_Powder_physics_storage_struct->vertical_angular_limit;
    
    
    // Should not interfer with the global variables
    double vx = k_initial[0]*K2V;
    double vy = k_initial[1]*K2V;
    double vz = k_initial[2]*K2V;
    
    double v = sqrt(vx*vx + vy*vy + vz*vz);
    
    int line;
    double arg;
    double theta;
    double alpha,alpha0;
    
    double vout_x,vout_y,vout_z;
    double tmp_vx,tmp_vy,tmp_vz;
    double nx,ny,nz;
    double my_s_n;
    
    // copy from PowderN component
    if (line_info->count > 0) {
          /* choose line */
          if (line_info->Nq > 1) line=floor(line_info->Nq*rand01());  /* Select between Nq powder lines */
          else line = 0;
          if (line_info->w_v[line])
            arg = line_info->q_v[line]*(1+line_info->w_v[line]*randnorm())/(2.0*v);
          else
            arg = line_info->q_v[line]/(2.0*v);
            my_s_n = line_info->my_s_v2[line]/(v*v);
          if(fabs(arg) > 1) {
            //printf("Powder scattering function returned 0, should not happen\n");
            return 0; /* No bragg scattering possible (was absorb)*/
          }
          theta = asin(arg);          /* Bragg scattering law */

          /* Choose point on Debye-Scherrer cone */
          if (vertical_angular_limit)
          { /* relate height of detector to the height on DS cone */
            arg = sin(vertical_angular_limit*DEG2RAD/2)/sin(2*theta);
            /* If full Debye-Scherrer cone is within d_phi, don't focus */
            if (arg < -1 || arg > 1) vertical_angular_limit = 0;
            /* Otherwise, determine alpha to rotate from scattering plane
               into vertical_angular_limit focusing area*/
            else alpha = 2*asin(arg);
          }
          if (vertical_angular_limit) {
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

          /* now find a nearly vertical rotation axis:
           * Either
           *  (v along Z) x (X axis) -> nearly Y axis
           * Or
           *  (v along X) x (Z axis) -> nearly Y axis
           */
          if (fabs(scalar_prod(1,0,0,vx/v,vy/v,vz/v)) < fabs(scalar_prod(0,0,1,vx/v,vy/v,vz/v))) {
            nx = 1; ny = 0; nz = 0;
          } else {
            nx = 0; ny = 0; nz = 1;
          }
          vec_prod(tmp_vx,tmp_vy,tmp_vz, vx,vy,vz, nx,ny,nz);

          /* v_out = rotate 'v' by 2*theta around tmp_v: Bragg angle */
          rotate(vout_x,vout_y,vout_z, vx,vy,vz, 2*theta, tmp_vx,tmp_vy,tmp_vz);

          /* tmp_v = rotate v_out by alpha0 around 'v' (Debye-Scherrer cone) */
          rotate(tmp_vx,tmp_vy,tmp_vz, vout_x,vout_y,vout_z, alpha0, vx, vy, vz);
          vx = tmp_vx;
          vy = tmp_vy;
          vz = tmp_vz;
        
          k_final[0] = V2K*vx; k_final[1] = V2K*vy; k_final[2] = V2K*vz;
        
          //*weight *= line_info->Nq*my_s_n; I believe my_s_n is part of the correction for sampling posistion, not to be done here
          *weight *= line_info->Nq*my_s_n/line_info->my_s_v2_sum*v*v;
        
          //printf("my_s_n = %f \n",my_s_n);
        
          // What to do with my_s_n ?
          /*
          pmul  = line_info->Nq*l_full*my_s_n*exp(-(line_info->my_a_v/v+my_s)*(l+l_1))
                                  /(1-(p_inc+p_transmit));
          */
          // Correction in case of vertical_angular_limit focusing - BUT only when d_phi != 0
          if (vertical_angular_limit) *weight *= alpha/PI;
          
        
          line_info->type = 'c';
          line_info->dq = line_info->q_v[line]*V2K;
          
          
        } else {
        /* else transmit <-- No powder lines in file */
        printf("Error, need lines in the PowderN input file\n");
        }


    //printf("Powder scattering function returned 1\n");
    return 1;
};

#line 8356 "Union_manual_example.c"

/* Shared user declarations for all components 'Union_make_material'. */
#line 57 "/usr/share/mcstas/2.5/contrib/union/Union_make_material.comp"
#ifndef Union
#define Union $Revision: 0.8 $

#include "Union_functions.c"
#include "Union_initialization.c"

#endif

// This function checks if global_process_element should be included in this material when using automatic linking, returns 1 if yes, 0 if no.
int automatic_linking_materials_function(struct global_process_element_struct global_process_element, struct pointer_to_global_material_list global_material_list,int current_index) {
    // Remember this function is used before the current material is added to global_material_list
    // debug info
    //MPI_MASTER(
    //printf("Checking if process with index %d should be automatically linked to material with index %d\n",global_process_element.component_index,current_index);
    //)

    // Check if this is the first make_material, which makes the problem simpler.
    if (global_material_list.num_elements == 0) {
       if (global_process_element.component_index < current_index) return 1;
       else return 0;
    }
    // In case there are more than 1 make_material, global_material_list.elements[global_material_list.num_elements-1].component_index makes sense.
    if (global_process_element.component_index < current_index && global_process_element.component_index > global_material_list.elements[global_material_list.num_elements-1].component_index) return 1;
    else return 0;
}

void manual_linking_function_material(char *input_string, struct pointer_to_global_process_list *global_process_list, struct pointer_to_1d_int_list *accepted_processes, char *component_name) {
    // Need to check a input_string of text for an occurance of name. If it is in the inputstring, yes return 1, otherwise 0.
   char *token;
   int loop_index;
   char local_string[256];
   
   strcpy(local_string,input_string);
   // get the first token
   token = strtok(local_string,",");
   
   // walk through other tokens
   while( token != NULL ) 
   {
      //printf( " %s\n", token );
      for (loop_index=0;loop_index<global_process_list->num_elements;loop_index++) {
        if (strcmp(token,global_process_list->elements[loop_index].name) == 0) {
          add_element_to_int_list(accepted_processes,loop_index);
          break;
        }
        
        if (loop_index == global_process_list->num_elements - 1) {
          // All possible process names have been looked through, and the break was not executed.
          // Alert the user to this problem by showing the process name that was not found and the currently available processes
            printf("\n");
            printf("ERROR: The process string \"%s\" in Union material \"%s\" had an entry that did not match a specified process. \n",input_string,component_name);
            printf("       The unrecoignized process name was: \"%s\" \n",token);
            printf("       The processes available at this point (need to be defined before the material): \n");
            for (loop_index=0;loop_index<global_process_list->num_elements;loop_index++)
              printf("         %s\n",global_process_list->elements[loop_index].name);
            exit(1);
        }
      }
      
      // Updates the token
      token = strtok(NULL,",");
   }
}

// This function is needed in initialize of all geometry components
// Possible to insert these functions in make material, as they are only compiled once instead of many times
int manual_linking_function(char *name, char *input_string) {
    // Need to check a input_string of text for an occurance of name. If it is in the inputstring, yes return 1, otherwise 0.
   char *token;
   int return_integer=0;
   char local_string[124];
   
   strcpy(local_string,input_string);
   /* get the first token */
   token = strtok(local_string,",");
   
   /* walk through other tokens */
   while( token != NULL ) 
   {
      //printf( " %s\n", token );
      if (strcmp(token,name) == 0) return_integer=1;
      
      token = strtok(NULL,",");
   }
   
   return return_integer;
}



/*
int count_commas(char *string) {
  int return_value = 0;
  
  int index;
  for (index=0;index<strlen(string);index++) {
    printf("%c \n",string[index]);
    if (string[index]==',') return_value++;
  }
    
  //printf("number_of_commas = %d \n",return_value);
  return return_value;
}
*/


#line 8466 "Union_manual_example.c"

/* Shared user declarations for all components 'Guide_gravity'. */
#line 124 "/usr/share/mcstas/2.5/optics/Guide_gravity.comp"
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
#line 8822 "Union_manual_example.c"

/* Shared user declarations for all components 'Union_box'. */
#line 77 "/usr/share/mcstas/2.5/contrib/union/Union_box.comp"
#ifndef Union
#define Union $Revision: 0.8 $

#include "Union_functions.c"
#include "Union_initialization.c"

#endif


void mcdisplay_box_function(struct lines_to_draw *lines_to_draw_output,int index, struct geometry_struct **Geometries,int number_of_volumes) {
    // Function to call in mcdisplay section of the sample component for this volume
    // One can assume that Volumes[index] refers to a volume with the geometry described in this file
    
    double depth = Geometries[index]->geometry_parameters.p_box_storage->z_depth;
    double width1 = Geometries[index]->geometry_parameters.p_box_storage->x_width1;
    double width2 = Geometries[index]->geometry_parameters.p_box_storage->x_width2;
    double height1 = Geometries[index]->geometry_parameters.p_box_storage->y_height1;
    double height2 = Geometries[index]->geometry_parameters.p_box_storage->y_height2;
    
    Coords x_vector = Geometries[index]->geometry_parameters.p_box_storage->x_vector;
    Coords y_vector = Geometries[index]->geometry_parameters.p_box_storage->y_vector;
    Coords z_vector = Geometries[index]->geometry_parameters.p_box_storage->z_vector;
    
    Coords center = Geometries[index]->center;
    
    Coords square1[4],square2[4];
    
    square1[0] = coords_add(coords_add(coords_add(center,coords_scalar_mult(z_vector,-0.5*depth)),coords_scalar_mult(x_vector,-0.5*width1)),coords_scalar_mult(y_vector,-0.5*height1));
    
    square1[1] = coords_add(square1[0],coords_scalar_mult(x_vector,width1));
    square1[2] = coords_add(square1[1],coords_scalar_mult(y_vector,height1));
    square1[3] = coords_add(square1[0],coords_scalar_mult(y_vector,height1));
    
    square2[0] = coords_add(coords_add(coords_add(center,coords_scalar_mult(z_vector,0.5*depth)),coords_scalar_mult(x_vector,-0.5*width2)),coords_scalar_mult(y_vector,-0.5*height2));
    
    square2[1] = coords_add(square2[0],coords_scalar_mult(x_vector,width2));
    square2[2] = coords_add(square2[1],coords_scalar_mult(y_vector,height2));
    square2[3] = coords_add(square2[0],coords_scalar_mult(y_vector,height2));
    
    struct lines_to_draw lines_to_draw_temp;
    lines_to_draw_temp.number_of_lines = 0;
    
    int iterate;
    for (iterate=0;iterate<3;iterate++) {
        lines_to_draw_temp = draw_line_with_highest_priority(square1[iterate],square1[iterate+1],index,Geometries,number_of_volumes,2);
        merge_lines_to_draw(lines_to_draw_output,&lines_to_draw_temp);
    }
    lines_to_draw_temp = draw_line_with_highest_priority(square1[3],square1[0],index,Geometries,number_of_volumes,2);
    merge_lines_to_draw(lines_to_draw_output,&lines_to_draw_temp);

    for (iterate=0;iterate<3;iterate++) {
        lines_to_draw_temp = draw_line_with_highest_priority(square2[iterate],square2[iterate+1],index,Geometries,number_of_volumes,2);
        merge_lines_to_draw(lines_to_draw_output,&lines_to_draw_temp);
    }
    lines_to_draw_temp = draw_line_with_highest_priority(square2[3],square2[0],index,Geometries,number_of_volumes,2);
    merge_lines_to_draw(lines_to_draw_output,&lines_to_draw_temp);

    for (iterate=0;iterate<4;iterate++) {
        lines_to_draw_temp = draw_line_with_highest_priority(square1[iterate],square2[iterate],index,Geometries,number_of_volumes,2);
        merge_lines_to_draw(lines_to_draw_output,&lines_to_draw_temp);
    }
};

void initialize_box_geometry_from_main_component(struct geometry_struct *box) {
    // Function to be called in initialize of the main component
    // This is done as the rotation matrix needs to be relative to the main component instead of global
    // Everything done in initialize in this component file has the rotation matrix relative to global
    Coords simple_vector = coords_set(1,0,0);
    Coords rotated_vector;
    
    rotated_vector = rot_apply(box->rotation_matrix,simple_vector);
    NORM(rotated_vector.x,rotated_vector.y,rotated_vector.z);
    box->geometry_parameters.p_box_storage->x_vector = rotated_vector;
    
    simple_vector = coords_set(0,1,0);
    rotated_vector = rot_apply(box->rotation_matrix,simple_vector);
    NORM(rotated_vector.x,rotated_vector.y,rotated_vector.z);
    box->geometry_parameters.p_box_storage->y_vector = rotated_vector;
    
    simple_vector = coords_set(0,0,1);
    rotated_vector = rot_apply(box->rotation_matrix,simple_vector);
    NORM(rotated_vector.x,rotated_vector.y,rotated_vector.z);
    box->geometry_parameters.p_box_storage->z_vector = rotated_vector;
};

struct pointer_to_1d_coords_list box_shell_points(struct geometry_struct *geometry,int max_number_of_points) {
  // This function returns an array of corner positions for the box in the main coordinate system.
  // Normally one would limit it to a maximum number of points, but as there are only 8 for the box,
  //  it is hardcoded to 8. Other geometries can be approximated with a variable number of points.
  
  struct pointer_to_1d_coords_list corner_points;
  corner_points.elements = malloc(8*sizeof(Coords));
  corner_points.num_elements = 8;
  
  double depth = geometry->geometry_parameters.p_box_storage->z_depth;
  double width1 = geometry->geometry_parameters.p_box_storage->x_width1;
  double width2 = geometry->geometry_parameters.p_box_storage->x_width2;
  double height1 = geometry->geometry_parameters.p_box_storage->y_height1;
  double height2 = geometry->geometry_parameters.p_box_storage->y_height2;
    
  Coords x_vector = geometry->geometry_parameters.p_box_storage->x_vector;
  Coords y_vector = geometry->geometry_parameters.p_box_storage->y_vector;
  Coords z_vector = geometry->geometry_parameters.p_box_storage->z_vector;
    
  Coords center = geometry->center;
    
  corner_points.elements[0] = coords_add(coords_add(coords_add(center,coords_scalar_mult(z_vector,-0.5*depth)),coords_scalar_mult(x_vector,-0.5*width1)),coords_scalar_mult(y_vector,-0.5*height1));
    
  corner_points.elements[1] = coords_add(corner_points.elements[0],coords_scalar_mult(x_vector,width1));
  corner_points.elements[2] = coords_add(corner_points.elements[1],coords_scalar_mult(y_vector,height1));
  corner_points.elements[3] = coords_add(corner_points.elements[0],coords_scalar_mult(y_vector,height1));
    
  corner_points.elements[4] = coords_add(coords_add(coords_add(center,coords_scalar_mult(z_vector,0.5*depth)),coords_scalar_mult(x_vector,-0.5*width2)),coords_scalar_mult(y_vector,-0.5*height2));
    
  corner_points.elements[5] = coords_add(corner_points.elements[4],coords_scalar_mult(x_vector,width2));
  corner_points.elements[6] = coords_add(corner_points.elements[5],coords_scalar_mult(y_vector,height2));
  corner_points.elements[7] = coords_add(corner_points.elements[4],coords_scalar_mult(y_vector,height2));
  
  return corner_points;

}

#line 8948 "Union_manual_example.c"

/* Shared user declarations for all components 'Union_cylinder'. */
#line 73 "/usr/share/mcstas/2.5/contrib/union/Union_cylinder.comp"
#ifndef Union
#define Union $Revision: 0.8 $

#include "Union_functions.c"
#include "Union_initialization.c"

#endif


void mcdisplay_cylinder_function(struct lines_to_draw *lines_to_draw_output,int index, struct geometry_struct **Geometries,int number_of_volumes) {
    // Function to call in mcdisplay section of the sample component for this volume
    // One can assume that Geometries[index] refers to a geometry as described in this file
    // The 4 lines describin the cylinders side are aligned to the local frame of the cylinder,
    //   it would be nicer to have them alligned with the global frame so that they show up nicely in
    //   pgplotters on mcdisplay.
    // One could get the current global rotation and use this to counteract this effect.
    
    double height = Geometries[index]->geometry_parameters.p_cylinder_storage->height;
    double radius = Geometries[index]->geometry_parameters.p_cylinder_storage->cyl_radius;
    Coords direction = Geometries[index]->geometry_parameters.p_cylinder_storage->direction_vector;
    Coords center = Geometries[index]->center;
    
    Coords bottom_point = coords_add(center,coords_scalar_mult(direction,0.5*height));
    Coords top_point = coords_add(center,coords_scalar_mult(direction,-0.5*height));
    
    struct lines_to_draw lines_to_draw_temp;
    lines_to_draw_temp.number_of_lines = 0;
    
    lines_to_draw_temp = draw_circle_with_highest_priority(top_point,direction,radius,index,Geometries,number_of_volumes,2);
    merge_lines_to_draw(lines_to_draw_output,&lines_to_draw_temp);

    lines_to_draw_temp = draw_circle_with_highest_priority(bottom_point,direction,radius,index,Geometries,number_of_volumes,2);
    merge_lines_to_draw(lines_to_draw_output,&lines_to_draw_temp);
    
    Coords point1,point2;
    int iterate,number_of_points=4;
    
    for (iterate=0;iterate<number_of_points;iterate++) {
        point1 = point_on_circle(top_point,direction,radius,iterate,number_of_points);
        point2 = point_on_circle(bottom_point,direction,radius,iterate,number_of_points);
        lines_to_draw_temp = draw_line_with_highest_priority(point1,point2,index,Geometries,number_of_volumes,2);
        merge_lines_to_draw(lines_to_draw_output,&lines_to_draw_temp);
    }
};

void initialize_cylinder_geometry_from_main_component(struct geometry_struct *cylinder) {
    // Function to be called in initialize of the main component
    // This is done as the rotation matrix needs to be relative to the main component instead of global
    // Everything done in initialize in this component file has the rotation matrix relative to global
    
    Coords simple_vector;
    Coords cyl_vector;
    
    // Start with vector that points along the cylinder in the local frame
    simple_vector = coords_set(0,1,0);

    // Rotate the direction vector of the cylinder to the master component frame of reference
    cyl_vector = rot_apply(cylinder->rotation_matrix,simple_vector);
    NORM(cyl_vector.x,cyl_vector.y,cyl_vector.z);
    cylinder->geometry_parameters.p_cylinder_storage->direction_vector.x = cyl_vector.x;
    cylinder->geometry_parameters.p_cylinder_storage->direction_vector.y = cyl_vector.y;
    cylinder->geometry_parameters.p_cylinder_storage->direction_vector.z = cyl_vector.z;
    // if (verbal == 1) printf("Cords vector1 = (%f,%f,%f)\n",cyl_vector.x,cyl_vector.y,
}

struct pointer_to_1d_coords_list cylinder_shell_points(struct geometry_struct *geometry,int max_number_of_points) {
  // Function that returns a number (less than max) of points on the geometry surface
  // If used, remember to free the space allocated.
  int points_per_circle = floor(max_number_of_points/2);
  
  struct pointer_to_1d_coords_list cylinder_shell_array;
  cylinder_shell_array.elements = malloc(2*points_per_circle*sizeof(Coords));
  cylinder_shell_array.num_elements = 2*points_per_circle;
  
  Coords cyl_direction = geometry->geometry_parameters.p_cylinder_storage->direction_vector;
  Coords center = geometry->center;
  double radius = geometry->geometry_parameters.p_cylinder_storage->cyl_radius;
  double height = geometry->geometry_parameters.p_cylinder_storage->height;
  
  Coords cyl_top_point = coords_add(center,coords_scalar_mult(cyl_direction,0.5*height));
  Coords cyl_bottom_point = coords_add(center,coords_scalar_mult(cyl_direction,-0.5*height));
  
  points_on_circle(cylinder_shell_array.elements,cyl_top_point,cyl_direction,radius,points_per_circle);
  // Need to verify this pointer arithimatic works as intended
  points_on_circle(cylinder_shell_array.elements+points_per_circle,cyl_bottom_point,cyl_direction,radius,points_per_circle);
  
  return cylinder_shell_array;
}

#line 9041 "Union_manual_example.c"

/* Shared user declarations for all components 'Union_master'. */
#line 53 "/usr/share/mcstas/2.5/contrib/union/Union_master.comp"
#ifndef Union
#define Union $Revision: 0.9 $

#include "Union_functions.c"
#include "Union_initialization.c"

#endif
#line 9052 "Union_manual_example.c"

/* Shared user declarations for all components 'Monitor_nD'. */
#line 214 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
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
          { Set_Vars_Coord_Type = DEFS->COORD_V; strcpy(Set_Vars_Coord_Label,"Velocity [m/s]"); strcpy(Set_Vars_Coord_Var,"v"); lmin = 0; lmax = 10000; }
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



#line 10996 "Union_manual_example.c"

/* Instrument parameters. */

#define mcNUMIPAR 0
int mcnumipar = 0;
struct mcinputtable_struct mcinputtable[mcNUMIPAR+1] = {
  NULL, NULL, instr_type_double, ""
};

/* User declarations from instrument definition. */
#define mccompcurname  Union_manual_example
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaUnion_manual_example coords_set(0,0,0)
#undef mcposaUnion_manual_example
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

/* Setting parameters for component 'Al_inc' [1]. */
MCNUM mccAl_inc_sigma;
MCNUM mccAl_inc_packing_factor;
MCNUM mccAl_inc_unit_cell_volume;
MCNUM mccAl_inc_interact_fraction;

/* Definition parameters for component 'Al_powder' [2]. */
#define mccAl_powder_format Undefined
/* Setting parameters for component 'Al_powder' [2]. */
char mccAl_powder_reflections[16384];
MCNUM mccAl_powder_packing_factor;
MCNUM mccAl_powder_Vc;
MCNUM mccAl_powder_delta_d_d;
MCNUM mccAl_powder_DW;
MCNUM mccAl_powder_nb_atoms;
MCNUM mccAl_powder_d_phi;
MCNUM mccAl_powder_density;
MCNUM mccAl_powder_weight;
MCNUM mccAl_powder_barns;
MCNUM mccAl_powder_Strain;
MCNUM mccAl_powder_interact_fraction;

/* Setting parameters for component 'Al' [3]. */
char mccAl_process_string[16384];
MCNUM mccAl_my_absorption;
MCNUM mccAl_absorber;

/* Setting parameters for component 'Fe_inc' [4]. */
MCNUM mccFe_inc_sigma;
MCNUM mccFe_inc_packing_factor;
MCNUM mccFe_inc_unit_cell_volume;
MCNUM mccFe_inc_interact_fraction;

/* Definition parameters for component 'Fe_powder' [5]. */
#define mccFe_powder_format Undefined
/* Setting parameters for component 'Fe_powder' [5]. */
char mccFe_powder_reflections[16384];
MCNUM mccFe_powder_packing_factor;
MCNUM mccFe_powder_Vc;
MCNUM mccFe_powder_delta_d_d;
MCNUM mccFe_powder_DW;
MCNUM mccFe_powder_nb_atoms;
MCNUM mccFe_powder_d_phi;
MCNUM mccFe_powder_density;
MCNUM mccFe_powder_weight;
MCNUM mccFe_powder_barns;
MCNUM mccFe_powder_Strain;
MCNUM mccFe_powder_interact_fraction;

/* Setting parameters for component 'Fe' [6]. */
char mccFe_process_string[16384];
MCNUM mccFe_my_absorption;
MCNUM mccFe_absorber;

/* Setting parameters for component 'Progress' [7]. */
char mccProgress_profile[16384];
MCNUM mccProgress_percent;
MCNUM mccProgress_flag_save;
MCNUM mccProgress_minutes;

/* Setting parameters for component 'Source' [8]. */
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

/* Setting parameters for component 'Guide' [9]. */
MCNUM mccGuide_w1;
MCNUM mccGuide_h1;
MCNUM mccGuide_w2;
MCNUM mccGuide_h2;
MCNUM mccGuide_l;
MCNUM mccGuide_R0;
MCNUM mccGuide_Qc;
MCNUM mccGuide_alpha;
MCNUM mccGuide_m;
MCNUM mccGuide_W;
MCNUM mccGuide_nslit;
MCNUM mccGuide_d;
MCNUM mccGuide_mleft;
MCNUM mccGuide_mright;
MCNUM mccGuide_mtop;
MCNUM mccGuide_mbottom;
MCNUM mccGuide_nhslit;
MCNUM mccGuide_G;
MCNUM mccGuide_aleft;
MCNUM mccGuide_aright;
MCNUM mccGuide_atop;
MCNUM mccGuide_abottom;
MCNUM mccGuide_wavy;
MCNUM mccGuide_wavy_z;
MCNUM mccGuide_wavy_tb;
MCNUM mccGuide_wavy_lr;
MCNUM mccGuide_chamfers;
MCNUM mccGuide_chamfers_z;
MCNUM mccGuide_chamfers_lr;
MCNUM mccGuide_chamfers_tb;
MCNUM mccGuide_nelements;
MCNUM mccGuide_nu;
MCNUM mccGuide_phase;
char mccGuide_reflect[16384];

/* Setting parameters for component 'Hexagonal_container_1' [11]. */
char mccHexagonal_container_1_material_string[16384];
MCNUM mccHexagonal_container_1_priority;
MCNUM mccHexagonal_container_1_xwidth;
MCNUM mccHexagonal_container_1_yheight;
MCNUM mccHexagonal_container_1_zdepth;
MCNUM mccHexagonal_container_1_xwidth2;
MCNUM mccHexagonal_container_1_yheight2;
MCNUM mccHexagonal_container_1_visualize;
int mccHexagonal_container_1_target_index;
MCNUM mccHexagonal_container_1_target_x;
MCNUM mccHexagonal_container_1_target_y;
MCNUM mccHexagonal_container_1_target_z;
MCNUM mccHexagonal_container_1_focus_aw;
MCNUM mccHexagonal_container_1_focus_ah;
MCNUM mccHexagonal_container_1_focus_xw;
MCNUM mccHexagonal_container_1_focus_xh;
MCNUM mccHexagonal_container_1_focus_r;
MCNUM mccHexagonal_container_1_p_interact;
char mccHexagonal_container_1_mask_string[16384];
char mccHexagonal_container_1_mask_setting[16384];
MCNUM mccHexagonal_container_1_number_of_activations;

/* Setting parameters for component 'Hexagonal_container_2' [12]. */
char mccHexagonal_container_2_material_string[16384];
MCNUM mccHexagonal_container_2_priority;
MCNUM mccHexagonal_container_2_xwidth;
MCNUM mccHexagonal_container_2_yheight;
MCNUM mccHexagonal_container_2_zdepth;
MCNUM mccHexagonal_container_2_xwidth2;
MCNUM mccHexagonal_container_2_yheight2;
MCNUM mccHexagonal_container_2_visualize;
int mccHexagonal_container_2_target_index;
MCNUM mccHexagonal_container_2_target_x;
MCNUM mccHexagonal_container_2_target_y;
MCNUM mccHexagonal_container_2_target_z;
MCNUM mccHexagonal_container_2_focus_aw;
MCNUM mccHexagonal_container_2_focus_ah;
MCNUM mccHexagonal_container_2_focus_xw;
MCNUM mccHexagonal_container_2_focus_xh;
MCNUM mccHexagonal_container_2_focus_r;
MCNUM mccHexagonal_container_2_p_interact;
char mccHexagonal_container_2_mask_string[16384];
char mccHexagonal_container_2_mask_setting[16384];
MCNUM mccHexagonal_container_2_number_of_activations;

/* Setting parameters for component 'Sample' [13]. */
char mccSample_material_string[16384];
MCNUM mccSample_priority;
MCNUM mccSample_radius;
MCNUM mccSample_yheight;
MCNUM mccSample_visualize;
int mccSample_target_index;
MCNUM mccSample_target_x;
MCNUM mccSample_target_y;
MCNUM mccSample_target_z;
MCNUM mccSample_focus_aw;
MCNUM mccSample_focus_ah;
MCNUM mccSample_focus_xw;
MCNUM mccSample_focus_xh;
MCNUM mccSample_focus_r;
MCNUM mccSample_p_interact;
char mccSample_mask_string[16384];
char mccSample_mask_setting[16384];
MCNUM mccSample_number_of_activations;

/* Setting parameters for component 'Master' [14]. */
MCNUM mccMaster_allow_inside_start;
MCNUM mccMaster_history_limit;

/* Definition parameters for component 'Banana_monitor' [15]. */
#define mccBanana_monitor_user1 FLT_MAX
#define mccBanana_monitor_user2 FLT_MAX
#define mccBanana_monitor_user3 FLT_MAX
/* Setting parameters for component 'Banana_monitor' [15]. */
MCNUM mccBanana_monitor_xwidth;
MCNUM mccBanana_monitor_yheight;
MCNUM mccBanana_monitor_zdepth;
MCNUM mccBanana_monitor_xmin;
MCNUM mccBanana_monitor_xmax;
MCNUM mccBanana_monitor_ymin;
MCNUM mccBanana_monitor_ymax;
MCNUM mccBanana_monitor_zmin;
MCNUM mccBanana_monitor_zmax;
MCNUM mccBanana_monitor_bins;
MCNUM mccBanana_monitor_min;
MCNUM mccBanana_monitor_max;
MCNUM mccBanana_monitor_restore_neutron;
MCNUM mccBanana_monitor_radius;
char mccBanana_monitor_options[16384];
char mccBanana_monitor_filename[16384];
char mccBanana_monitor_geometry[16384];
char mccBanana_monitor_username1[16384];
char mccBanana_monitor_username2[16384];
char mccBanana_monitor_username3[16384];
int mccBanana_monitor_nowritefile;

/* User component declarations. */

/* User declarations for component 'Al_inc' [1]. */
#define mccompcurname  Al_inc
#define mccompcurtype  Incoherent_process
#define mccompcurindex 1
#define This_process mccAl_inc_This_process
#define Incoherent_storage mccAl_inc_Incoherent_storage
#define effective_my_scattering mccAl_inc_effective_my_scattering
#define sigma mccAl_inc_sigma
#define packing_factor mccAl_inc_packing_factor
#define unit_cell_volume mccAl_inc_unit_cell_volume
#define interact_fraction mccAl_inc_interact_fraction
#line 104 "/usr/share/mcstas/2.5/contrib/union/Incoherent_process.comp"
// Needed for transport to the main component
struct global_process_element_struct global_process_element;
struct scattering_process_struct This_process;

#ifndef PROCESS_DETECTOR
	//struct pointer_to_global_process_list global_process_list = {0,NULL};
	#define PROCESS_DETECTOR dummy
#endif

// Declare for this component, to do calculations on the input / store in the transported data
struct Incoherent_physics_storage_struct Incoherent_storage;
double effective_my_scattering;

#line 11270 "Union_manual_example.c"
#undef interact_fraction
#undef unit_cell_volume
#undef packing_factor
#undef sigma
#undef effective_my_scattering
#undef Incoherent_storage
#undef This_process
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Al_powder' [2]. */
#define mccompcurname  Al_powder
#define mccompcurtype  Powder_process
#define mccompcurindex 2
#define format mccAl_powder_format
#define This_process mccAl_powder_This_process
#define Powder_storage mccAl_powder_Powder_storage
#define effective_my_scattering mccAl_powder_effective_my_scattering
#define d_phi mccAl_powder_d_phi
#define line_info mccAl_powder_line_info
#define columns mccAl_powder_columns
#define reflections mccAl_powder_reflections
#define packing_factor mccAl_powder_packing_factor
#define Vc mccAl_powder_Vc
#define delta_d_d mccAl_powder_delta_d_d
#define DW mccAl_powder_DW
#define nb_atoms mccAl_powder_nb_atoms
#define d_phi mccAl_powder_d_phi
#define density mccAl_powder_density
#define weight mccAl_powder_weight
#define barns mccAl_powder_barns
#define Strain mccAl_powder_Strain
#define interact_fraction mccAl_powder_interact_fraction
#line 617 "/usr/share/mcstas/2.5/contrib/union/Powder_process.comp"
// Needed for transport to the main component
struct global_process_element_struct global_process_element;
struct scattering_process_struct This_process;

#ifndef PROCESS_DETECTOR
	//struct pointer_to_global_process_list global_process_list = {0,NULL};
	#define PROCESS_DETECTOR dummy
#endif

// Declare for this component, to do calculations on the input / store in the transported data
struct Powder_physics_storage_struct Powder_storage;
struct line_info_struct_union line_info;
double effective_my_scattering;

int columns[9] = format;


#line 11323 "Union_manual_example.c"
#undef interact_fraction
#undef Strain
#undef barns
#undef weight
#undef density
#undef d_phi
#undef nb_atoms
#undef DW
#undef delta_d_d
#undef Vc
#undef packing_factor
#undef reflections
#undef columns
#undef line_info
#undef d_phi
#undef effective_my_scattering
#undef Powder_storage
#undef This_process
#undef format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Al' [3]. */
#define mccompcurname  Al
#define mccompcurtype  Union_make_material
#define mccompcurindex 3
#define loop_index mccAl_loop_index
#define this_material mccAl_this_material
#define accepted_processes mccAl_accepted_processes
#define global_material_element mccAl_global_material_element
#define process_string mccAl_process_string
#define my_absorption mccAl_my_absorption
#define absorber mccAl_absorber
#line 167 "/usr/share/mcstas/2.5/contrib/union/Union_make_material.comp"
// Needed for transport to the main component
struct global_material_element_struct global_material_element;
struct physics_struct this_material;

#ifndef MATERIAL_DETECTOR
	//struct pointer_to_global_material_list global_material_list = {0,NULL};
	#define MATERIAL_DETECTOR dummy
#endif

int loop_index;
int found_process;
int specified_processes;
char local_string[256];

struct pointer_to_1d_int_list accepted_processes = {0,NULL};

// Add setup for loggers since make_material is called before any volume / master
#ifndef UNION_LOGGER_DECLARE
	//struct pointer_to_global_logger_list global_all_volume_logger_list = {0,NULL};
    //struct pointer_to_global_logger_list global_specific_volumes_logger_list = {0,NULL};
	#define UNION_LOGGER_DECLARE dummy
#endif

#line 11382 "Union_manual_example.c"
#undef absorber
#undef my_absorption
#undef process_string
#undef global_material_element
#undef accepted_processes
#undef this_material
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Fe_inc' [4]. */
#define mccompcurname  Fe_inc
#define mccompcurtype  Incoherent_process
#define mccompcurindex 4
#define This_process mccFe_inc_This_process
#define Incoherent_storage mccFe_inc_Incoherent_storage
#define effective_my_scattering mccFe_inc_effective_my_scattering
#define sigma mccFe_inc_sigma
#define packing_factor mccFe_inc_packing_factor
#define unit_cell_volume mccFe_inc_unit_cell_volume
#define interact_fraction mccFe_inc_interact_fraction
#line 104 "/usr/share/mcstas/2.5/contrib/union/Incoherent_process.comp"
// Needed for transport to the main component
struct global_process_element_struct global_process_element;
struct scattering_process_struct This_process;

#ifndef PROCESS_DETECTOR
	//struct pointer_to_global_process_list global_process_list = {0,NULL};
	#define PROCESS_DETECTOR dummy
#endif

// Declare for this component, to do calculations on the input / store in the transported data
struct Incoherent_physics_storage_struct Incoherent_storage;
double effective_my_scattering;

#line 11419 "Union_manual_example.c"
#undef interact_fraction
#undef unit_cell_volume
#undef packing_factor
#undef sigma
#undef effective_my_scattering
#undef Incoherent_storage
#undef This_process
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Fe_powder' [5]. */
#define mccompcurname  Fe_powder
#define mccompcurtype  Powder_process
#define mccompcurindex 5
#define format mccFe_powder_format
#define This_process mccFe_powder_This_process
#define Powder_storage mccFe_powder_Powder_storage
#define effective_my_scattering mccFe_powder_effective_my_scattering
#define d_phi mccFe_powder_d_phi
#define line_info mccFe_powder_line_info
#define columns mccFe_powder_columns
#define reflections mccFe_powder_reflections
#define packing_factor mccFe_powder_packing_factor
#define Vc mccFe_powder_Vc
#define delta_d_d mccFe_powder_delta_d_d
#define DW mccFe_powder_DW
#define nb_atoms mccFe_powder_nb_atoms
#define d_phi mccFe_powder_d_phi
#define density mccFe_powder_density
#define weight mccFe_powder_weight
#define barns mccFe_powder_barns
#define Strain mccFe_powder_Strain
#define interact_fraction mccFe_powder_interact_fraction
#line 617 "/usr/share/mcstas/2.5/contrib/union/Powder_process.comp"
// Needed for transport to the main component
struct global_process_element_struct global_process_element;
struct scattering_process_struct This_process;

#ifndef PROCESS_DETECTOR
	//struct pointer_to_global_process_list global_process_list = {0,NULL};
	#define PROCESS_DETECTOR dummy
#endif

// Declare for this component, to do calculations on the input / store in the transported data
struct Powder_physics_storage_struct Powder_storage;
struct line_info_struct_union line_info;
double effective_my_scattering;

int columns[9] = format;


#line 11472 "Union_manual_example.c"
#undef interact_fraction
#undef Strain
#undef barns
#undef weight
#undef density
#undef d_phi
#undef nb_atoms
#undef DW
#undef delta_d_d
#undef Vc
#undef packing_factor
#undef reflections
#undef columns
#undef line_info
#undef d_phi
#undef effective_my_scattering
#undef Powder_storage
#undef This_process
#undef format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Fe' [6]. */
#define mccompcurname  Fe
#define mccompcurtype  Union_make_material
#define mccompcurindex 6
#define loop_index mccFe_loop_index
#define this_material mccFe_this_material
#define accepted_processes mccFe_accepted_processes
#define global_material_element mccFe_global_material_element
#define process_string mccFe_process_string
#define my_absorption mccFe_my_absorption
#define absorber mccFe_absorber
#line 167 "/usr/share/mcstas/2.5/contrib/union/Union_make_material.comp"
// Needed for transport to the main component
struct global_material_element_struct global_material_element;
struct physics_struct this_material;

#ifndef MATERIAL_DETECTOR
	//struct pointer_to_global_material_list global_material_list = {0,NULL};
	#define MATERIAL_DETECTOR dummy
#endif

int loop_index;
int found_process;
int specified_processes;
char local_string[256];

struct pointer_to_1d_int_list accepted_processes = {0,NULL};

// Add setup for loggers since make_material is called before any volume / master
#ifndef UNION_LOGGER_DECLARE
	//struct pointer_to_global_logger_list global_all_volume_logger_list = {0,NULL};
    //struct pointer_to_global_logger_list global_specific_volumes_logger_list = {0,NULL};
	#define UNION_LOGGER_DECLARE dummy
#endif

#line 11531 "Union_manual_example.c"
#undef absorber
#undef my_absorption
#undef process_string
#undef global_material_element
#undef accepted_processes
#undef this_material
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Progress' [7]. */
#define mccompcurname  Progress
#define mccompcurtype  Progress_bar
#define mccompcurindex 7
#define IntermediateCnts mccProgress_IntermediateCnts
#define StartTime mccProgress_StartTime
#define EndTime mccProgress_EndTime
#define CurrentTime mccProgress_CurrentTime
#define profile mccProgress_profile
#define percent mccProgress_percent
#define flag_save mccProgress_flag_save
#define minutes mccProgress_minutes
#line 44 "/usr/share/mcstas/2.5/misc/Progress_bar.comp"
#ifndef PROGRESS_BAR
#define PROGRESS_BAR
#else
#error Only one Progress_bar component may be used in an instrument definition.
#endif

double IntermediateCnts;
time_t StartTime;
time_t EndTime;
time_t CurrentTime;
#line 11566 "Union_manual_example.c"
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

/* User declarations for component 'Source' [8]. */
#define mccompcurname  Source
#define mccompcurtype  Source_simple
#define mccompcurindex 8
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
#line 60 "/usr/share/mcstas/2.5/sources/Source_simple.comp"
double pmul, srcArea;
int square;
double tx,ty,tz;
#line 11603 "Union_manual_example.c"
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

/* User declarations for component 'Guide' [9]. */
#define mccompcurname  Guide
#define mccompcurtype  Guide_gravity
#define mccompcurindex 9
#define GVars mccGuide_GVars
#define pTable mccGuide_pTable
#define w1 mccGuide_w1
#define h1 mccGuide_h1
#define w2 mccGuide_w2
#define h2 mccGuide_h2
#define l mccGuide_l
#define R0 mccGuide_R0
#define Qc mccGuide_Qc
#define alpha mccGuide_alpha
#define m mccGuide_m
#define W mccGuide_W
#define nslit mccGuide_nslit
#define d mccGuide_d
#define mleft mccGuide_mleft
#define mright mccGuide_mright
#define mtop mccGuide_mtop
#define mbottom mccGuide_mbottom
#define nhslit mccGuide_nhslit
#define G mccGuide_G
#define aleft mccGuide_aleft
#define aright mccGuide_aright
#define atop mccGuide_atop
#define abottom mccGuide_abottom
#define wavy mccGuide_wavy
#define wavy_z mccGuide_wavy_z
#define wavy_tb mccGuide_wavy_tb
#define wavy_lr mccGuide_wavy_lr
#define chamfers mccGuide_chamfers
#define chamfers_z mccGuide_chamfers_z
#define chamfers_lr mccGuide_chamfers_lr
#define chamfers_tb mccGuide_chamfers_tb
#define nelements mccGuide_nelements
#define nu mccGuide_nu
#define phase mccGuide_phase
#define reflect mccGuide_reflect
#line 334 "/usr/share/mcstas/2.5/optics/Guide_gravity.comp"
  Gravity_guide_Vars_type GVars;
  t_Table pTable;
#line 11667 "Union_manual_example.c"
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

/* User declarations for component 'sample_position' [10]. */
#define mccompcurname  sample_position
#define mccompcurtype  Arm
#define mccompcurindex 10
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Hexagonal_container_1' [11]. */
#define mccompcurname  Hexagonal_container_1
#define mccompcurtype  Union_box
#define mccompcurindex 11
#define loop_index mccHexagonal_container_1_loop_index
#define this_box_volume mccHexagonal_container_1_this_box_volume
#define global_geometry_element mccHexagonal_container_1_global_geometry_element
#define this_box_storage mccHexagonal_container_1_this_box_storage
#define material_string mccHexagonal_container_1_material_string
#define priority mccHexagonal_container_1_priority
#define xwidth mccHexagonal_container_1_xwidth
#define yheight mccHexagonal_container_1_yheight
#define zdepth mccHexagonal_container_1_zdepth
#define xwidth2 mccHexagonal_container_1_xwidth2
#define yheight2 mccHexagonal_container_1_yheight2
#define visualize mccHexagonal_container_1_visualize
#define target_index mccHexagonal_container_1_target_index
#define target_x mccHexagonal_container_1_target_x
#define target_y mccHexagonal_container_1_target_y
#define target_z mccHexagonal_container_1_target_z
#define focus_aw mccHexagonal_container_1_focus_aw
#define focus_ah mccHexagonal_container_1_focus_ah
#define focus_xw mccHexagonal_container_1_focus_xw
#define focus_xh mccHexagonal_container_1_focus_xh
#define focus_r mccHexagonal_container_1_focus_r
#define p_interact mccHexagonal_container_1_p_interact
#define mask_string mccHexagonal_container_1_mask_string
#define mask_setting mccHexagonal_container_1_mask_setting
#define number_of_activations mccHexagonal_container_1_number_of_activations
#line 203 "/usr/share/mcstas/2.5/contrib/union/Union_box.comp"
// Needed for transport to the main component
struct global_geometry_element_struct global_geometry_element;

#ifndef ANY_GEOMETRY_DETECTOR_DECLARE
    #define ANY_GEOMETRY_DETECTOR_DECLARE dummy
	//struct pointer_to_global_geometry_list global_geometry_list = {0,NULL};
#endif

int loop_index,found_geometries;

double x_component;
double y_component;
double z_component;

struct Volume_struct this_box_volume;
struct box_storage this_box_storage;

#line 11763 "Union_manual_example.c"
#undef number_of_activations
#undef mask_setting
#undef mask_string
#undef p_interact
#undef focus_r
#undef focus_xh
#undef focus_xw
#undef focus_ah
#undef focus_aw
#undef target_z
#undef target_y
#undef target_x
#undef target_index
#undef visualize
#undef yheight2
#undef xwidth2
#undef zdepth
#undef yheight
#undef xwidth
#undef priority
#undef material_string
#undef this_box_storage
#undef global_geometry_element
#undef this_box_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Hexagonal_container_2' [12]. */
#define mccompcurname  Hexagonal_container_2
#define mccompcurtype  Union_box
#define mccompcurindex 12
#define loop_index mccHexagonal_container_2_loop_index
#define this_box_volume mccHexagonal_container_2_this_box_volume
#define global_geometry_element mccHexagonal_container_2_global_geometry_element
#define this_box_storage mccHexagonal_container_2_this_box_storage
#define material_string mccHexagonal_container_2_material_string
#define priority mccHexagonal_container_2_priority
#define xwidth mccHexagonal_container_2_xwidth
#define yheight mccHexagonal_container_2_yheight
#define zdepth mccHexagonal_container_2_zdepth
#define xwidth2 mccHexagonal_container_2_xwidth2
#define yheight2 mccHexagonal_container_2_yheight2
#define visualize mccHexagonal_container_2_visualize
#define target_index mccHexagonal_container_2_target_index
#define target_x mccHexagonal_container_2_target_x
#define target_y mccHexagonal_container_2_target_y
#define target_z mccHexagonal_container_2_target_z
#define focus_aw mccHexagonal_container_2_focus_aw
#define focus_ah mccHexagonal_container_2_focus_ah
#define focus_xw mccHexagonal_container_2_focus_xw
#define focus_xh mccHexagonal_container_2_focus_xh
#define focus_r mccHexagonal_container_2_focus_r
#define p_interact mccHexagonal_container_2_p_interact
#define mask_string mccHexagonal_container_2_mask_string
#define mask_setting mccHexagonal_container_2_mask_setting
#define number_of_activations mccHexagonal_container_2_number_of_activations
#line 203 "/usr/share/mcstas/2.5/contrib/union/Union_box.comp"
// Needed for transport to the main component
struct global_geometry_element_struct global_geometry_element;

#ifndef ANY_GEOMETRY_DETECTOR_DECLARE
    #define ANY_GEOMETRY_DETECTOR_DECLARE dummy
	//struct pointer_to_global_geometry_list global_geometry_list = {0,NULL};
#endif

int loop_index,found_geometries;

double x_component;
double y_component;
double z_component;

struct Volume_struct this_box_volume;
struct box_storage this_box_storage;

#line 11840 "Union_manual_example.c"
#undef number_of_activations
#undef mask_setting
#undef mask_string
#undef p_interact
#undef focus_r
#undef focus_xh
#undef focus_xw
#undef focus_ah
#undef focus_aw
#undef target_z
#undef target_y
#undef target_x
#undef target_index
#undef visualize
#undef yheight2
#undef xwidth2
#undef zdepth
#undef yheight
#undef xwidth
#undef priority
#undef material_string
#undef this_box_storage
#undef global_geometry_element
#undef this_box_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Sample' [13]. */
#define mccompcurname  Sample
#define mccompcurtype  Union_cylinder
#define mccompcurindex 13
#define loop_index mccSample_loop_index
#define this_cylinder_volume mccSample_this_cylinder_volume
#define global_geometry_element mccSample_global_geometry_element
#define this_cylinder_storage mccSample_this_cylinder_storage
#define material_string mccSample_material_string
#define priority mccSample_priority
#define radius mccSample_radius
#define yheight mccSample_yheight
#define visualize mccSample_visualize
#define target_index mccSample_target_index
#define target_x mccSample_target_x
#define target_y mccSample_target_y
#define target_z mccSample_target_z
#define focus_aw mccSample_focus_aw
#define focus_ah mccSample_focus_ah
#define focus_xw mccSample_focus_xw
#define focus_xh mccSample_focus_xh
#define focus_r mccSample_focus_r
#define p_interact mccSample_p_interact
#define mask_string mccSample_mask_string
#define mask_setting mccSample_mask_setting
#define number_of_activations mccSample_number_of_activations
#line 166 "/usr/share/mcstas/2.5/contrib/union/Union_cylinder.comp"
// Needed for transport to the main component
struct global_geometry_element_struct global_geometry_element;

#ifndef ANY_GEOMETRY_DETECTOR_DECLARE
    #define ANY_GEOMETRY_DETECTOR_DECLARE dummy
	//struct pointer_to_global_geometry_list global_geometry_list = {0,NULL};
#endif

int dummy;
int loop_index,found_geometries;
int loop_2_index;
int material_index;

struct Volume_struct this_cylinder_volume;
struct cylinder_storage this_cylinder_storage;
#line 11912 "Union_manual_example.c"
#undef number_of_activations
#undef mask_setting
#undef mask_string
#undef p_interact
#undef focus_r
#undef focus_xh
#undef focus_xw
#undef focus_ah
#undef focus_aw
#undef target_z
#undef target_y
#undef target_x
#undef target_index
#undef visualize
#undef yheight
#undef radius
#undef priority
#undef material_string
#undef this_cylinder_storage
#undef global_geometry_element
#undef this_cylinder_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Master' [14]. */
#define mccompcurname  Master
#define mccompcurtype  Union_master
#define mccompcurindex 14
#define verbal mccMaster_verbal
#define list_verbal mccMaster_list_verbal
#define trace_verbal mccMaster_trace_verbal
#define finally_verbal mccMaster_finally_verbal
#define starting_volume_warning mccMaster_starting_volume_warning
#define global_master_element mccMaster_global_master_element
#define previous_master_index mccMaster_previous_master_index
#define geometry_list_index mccMaster_geometry_list_index
#define intersection_time_table mccMaster_intersection_time_table
#define Volumes mccMaster_Volumes
#define Geometries mccMaster_Geometries
#define starting_lists mccMaster_starting_lists
#define r mccMaster_r
#define r_start mccMaster_r_start
#define v mccMaster_v
#define error_msg mccMaster_error_msg
#define component_error_msg mccMaster_component_error_msg
#define string_output mccMaster_string_output
#define number_of_volumes mccMaster_number_of_volumes
#define volume_index mccMaster_volume_index
#define process_index mccMaster_process_index
#define solutions mccMaster_solutions
#define max_number_of_processes mccMaster_max_number_of_processes
#define limit mccMaster_limit
#define solution mccMaster_solution
#define min_solution mccMaster_min_solution
#define min_volume mccMaster_min_volume
#define time_found mccMaster_time_found
#define intersection_time mccMaster_intersection_time
#define min_intersection_time mccMaster_min_intersection_time
#define process mccMaster_process
#define process_start mccMaster_process_start
#define my_trace mccMaster_my_trace
#define p_my_trace mccMaster_p_my_trace
#define my_trace_fraction_control mccMaster_my_trace_fraction_control
#define k mccMaster_k
#define k_new mccMaster_k_new
#define v_length mccMaster_v_length
#define my_sum mccMaster_my_sum
#define my_sum_plus_abs mccMaster_my_sum_plus_abs
#define culmative_probability mccMaster_culmative_probability
#define mc_prop mccMaster_mc_prop
#define time_to_scattering mccMaster_time_to_scattering
#define length_to_scattering mccMaster_length_to_scattering
#define length_to_boundery mccMaster_length_to_boundery
#define time_to_boundery mccMaster_time_to_boundery
#define selected_process mccMaster_selected_process
#define scattering_event mccMaster_scattering_event
#define time_propagated_without_scattering mccMaster_time_propagated_without_scattering
#define a_next_volume_found mccMaster_a_next_volume_found
#define next_volume mccMaster_next_volume
#define next_volume_priority mccMaster_next_volume_priority
#define done mccMaster_done
#define current_volume mccMaster_current_volume
#define number_of_solutions mccMaster_number_of_solutions
#define number_of_solutions_static mccMaster_number_of_solutions_static
#define check mccMaster_check
#define start mccMaster_start
#define intersection_with_children mccMaster_intersection_with_children
#define geometry_output mccMaster_geometry_output
#define tree_next_volume mccMaster_tree_next_volume
#define pre_allocated1 mccMaster_pre_allocated1
#define pre_allocated2 mccMaster_pre_allocated2
#define pre_allocated3 mccMaster_pre_allocated3
#define ray_position mccMaster_ray_position
#define ray_velocity mccMaster_ray_velocity
#define ray_velocity_final mccMaster_ray_velocity_final
#define volume_0_found mccMaster_volume_0_found
#define scattered_flag mccMaster_scattered_flag
#define scattered_flag_VP mccMaster_scattered_flag_VP
#define master_transposed_rotation_matrix mccMaster_master_transposed_rotation_matrix
#define temp_rotation_matrix mccMaster_temp_rotation_matrix
#define non_rotated_position mccMaster_non_rotated_position
#define rotated_position mccMaster_rotated_position
#define enable_tagging mccMaster_enable_tagging
#define stop_tagging_ray mccMaster_stop_tagging_ray
#define stop_creating_nodes mccMaster_stop_creating_nodes
#define enable_tagging_check mccMaster_enable_tagging_check
#define master_tagging_node_list mccMaster_master_tagging_node_list
#define current_tagging_node mccMaster_current_tagging_node
#define tagging_leaf_counter mccMaster_tagging_leaf_counter
#define number_of_scattering_events mccMaster_number_of_scattering_events
#define real_transmission_probability mccMaster_real_transmission_probability
#define mc_transmission_probability mccMaster_mc_transmission_probability
#define number_of_masks mccMaster_number_of_masks
#define number_of_masked_volumes mccMaster_number_of_masked_volumes
#define need_to_run_within_which_volume mccMaster_need_to_run_within_which_volume
#define mask_index_main mccMaster_mask_index_main
#define mask_iterate mccMaster_mask_iterate
#define mask_status_list mccMaster_mask_status_list
#define current_mask_intersect_list_status mccMaster_current_mask_intersect_list_status
#define mask_volume_index_list mccMaster_mask_volume_index_list
#define geometry_component_index_list mccMaster_geometry_component_index_list
#define Volume_copies_allocated mccMaster_Volume_copies_allocated
#define loggers_with_data_array mccMaster_loggers_with_data_array
#define p_old mccMaster_p_old
#define this_logger mccMaster_this_logger
#define conditional_status mccMaster_conditional_status
#define tagging_conditional_list mccMaster_tagging_conditional_list
#define free_tagging_conditioanl_list mccMaster_free_tagging_conditioanl_list
#define logger_conditional_extend_array mccMaster_logger_conditional_extend_array
#define tagging_conditional_extend mccMaster_tagging_conditional_extend
#define max_conditional_extend_index mccMaster_max_conditional_extend_index
#define safty_distance mccMaster_safty_distance
#define safty_distance2 mccMaster_safty_distance2
#define number_of_processes_array mccMaster_number_of_processes_array
#define allow_inside_start mccMaster_allow_inside_start
#define history_limit mccMaster_history_limit
#line 64 "/usr/share/mcstas/2.5/contrib/union/Union_master.comp"
// Setting declare section, setting these at compile time is done to let the compiler optimize easier, 1 for enable, 0 for disable.
int verbal = 1;
int list_verbal = 1;
int trace_verbal = 0;
int enable_tagging = 1;
int finally_verbal = 0;

// It is possible to surpress warnings on starting volume by setting this to 1
int starting_volume_warning = 0;

// Declare the global variables (not to be in output parameters)
  struct global_master_element_struct global_master_element;
  
  #ifndef MASTER_DETECTOR
	//struct pointer_to_global_master_list global_master_list = {0,NULL};
	#define MASTER_DETECTOR dummy
  #endif
  
  // variables used for assigning global information to local variables
  int previous_master_index,geometry_list_index;

  // The main structures used in this component
  struct intersection_time_table_struct intersection_time_table;
  struct Volume_struct **Volumes;
  struct geometry_struct **Geometries;
  struct Volume_struct **Volume_copies;
  struct starting_lists_struct starting_lists;
  
  // garbage collection for volume_copies
  struct pointer_to_1d_int_list Volume_copies_allocated;

  // Vectors in old format (still used by intersect function, will go to Coords in future)
  double r[3],r_start[3],v[3];

  // Error handling
  int error_msg, component_error_msg=0;
  
  
  // For verbal output
  char string_output[128];
  
  int number_of_volumes, volume_index, process_index;
  int iterate,solutions,max_number_of_processes,limit;
  int solution,min_solution,min_volume,time_found;
  double intersection_time,min_intersection_time;
  
  struct scattering_process_struct *process,*process_start;
  double *my_trace,*p_my_trace,*my_trace_fraction_control;
  double k[3],k_new[3],k_rotated[3];
  double v_length,my_sum,my_sum_plus_abs,culmative_probability,mc_prop,time_to_scattering;
  double length_to_scattering,length_to_boundery,length_to_boundery_fp,time_to_boundery;
  int selected_process,scattering_event;
  double time_propagated_without_scattering;
  
  int a_next_volume_found,next_volume;
  double next_volume_priority;
  
  int done,current_volume,ray_sucseeded;
  int *number_of_solutions;
  int number_of_solutions_static;
  int *check,*start;
  int intersection_with_children,geometry_output;
  
  // For within_which_volume
  int tree_next_volume;
  int *pre_allocated1,*pre_allocated2,*pre_allocated3;
  Coords ray_position,ray_velocity,ray_velocity_rotated,ray_velocity_final,wavevector,wavevector_rotated;
  int volume_0_found=0;
  
  int *scattered_flag;
  int **scattered_flag_VP;
  
  // For coordinate transformations
  Rotation master_transposed_rotation_matrix;
  Rotation temp_rotation_matrix;
  Coords non_rotated_position;
  Coords rotated_position;
  int non_isotropic_found;
  
  // For tagging
  struct list_of_tagging_tree_node_pointers master_tagging_node_list;
  struct tagging_tree_node_struct *current_tagging_node;
  
  int tagging_leaf_counter=0,stop_tagging_ray,stop_creating_nodes;
  int number_of_scattering_events;
  
  // For geometry p interact
  double real_transmission_probability,mc_transmission_probability;
  
  // Process p interact
  int number_of_process_interacts_set,index_of_lacking_process;
  double total_process_interact;
  
  // Volume nr -> component index
  struct pointer_to_1d_int_list geometry_component_index_list;
  
  // Masks
  struct pointer_to_1d_int_list mask_volume_index_list;
  int number_of_masks=0;
  int number_of_masked_volumes=0;
  struct pointer_to_1d_int_list mask_status_list;
  struct pointer_to_1d_int_list current_mask_intersect_list_status;
  int mask_index_main,mask_iterate;
  int *mask_start,*mask_check;
  int need_to_run_within_which_volume;
  
  // Loggers
  struct logger_with_data_struct loggers_with_data_array;
  int *number_of_processes_array;
  double p_old;
  int log_index,conditional_status;
  struct logger_struct *this_logger;
  //  union detector_pointer_union detector_pointer;
  
  // Conditionals
  struct conditional_list_struct *tagging_conditional_list;
  int *logger_conditional_extend_array;
  int max_conditional_extend_index;
  int tagging_conditional_extend;
  int free_tagging_conditioanl_list;
  
  // Reliability control
  // Safty distance is needed to avoid having ray positions closer to a wall than the precision of intersection functions
  double safty_distance,safty_distance2;
  
#line 12177 "Union_manual_example.c"
#undef history_limit
#undef allow_inside_start
#undef number_of_processes_array
#undef safty_distance2
#undef safty_distance
#undef max_conditional_extend_index
#undef tagging_conditional_extend
#undef logger_conditional_extend_array
#undef free_tagging_conditioanl_list
#undef tagging_conditional_list
#undef conditional_status
#undef this_logger
#undef p_old
#undef loggers_with_data_array
#undef Volume_copies_allocated
#undef geometry_component_index_list
#undef mask_volume_index_list
#undef current_mask_intersect_list_status
#undef mask_status_list
#undef mask_iterate
#undef mask_index_main
#undef need_to_run_within_which_volume
#undef number_of_masked_volumes
#undef number_of_masks
#undef mc_transmission_probability
#undef real_transmission_probability
#undef number_of_scattering_events
#undef tagging_leaf_counter
#undef current_tagging_node
#undef master_tagging_node_list
#undef enable_tagging_check
#undef stop_creating_nodes
#undef stop_tagging_ray
#undef enable_tagging
#undef rotated_position
#undef non_rotated_position
#undef temp_rotation_matrix
#undef master_transposed_rotation_matrix
#undef scattered_flag_VP
#undef scattered_flag
#undef volume_0_found
#undef ray_velocity_final
#undef ray_velocity
#undef ray_position
#undef pre_allocated3
#undef pre_allocated2
#undef pre_allocated1
#undef tree_next_volume
#undef geometry_output
#undef intersection_with_children
#undef start
#undef check
#undef number_of_solutions_static
#undef number_of_solutions
#undef current_volume
#undef done
#undef next_volume_priority
#undef next_volume
#undef a_next_volume_found
#undef time_propagated_without_scattering
#undef scattering_event
#undef selected_process
#undef time_to_boundery
#undef length_to_boundery
#undef length_to_scattering
#undef time_to_scattering
#undef mc_prop
#undef culmative_probability
#undef my_sum_plus_abs
#undef my_sum
#undef v_length
#undef k_new
#undef k
#undef my_trace_fraction_control
#undef p_my_trace
#undef my_trace
#undef process_start
#undef process
#undef min_intersection_time
#undef intersection_time
#undef time_found
#undef min_volume
#undef min_solution
#undef solution
#undef limit
#undef max_number_of_processes
#undef solutions
#undef process_index
#undef volume_index
#undef number_of_volumes
#undef string_output
#undef component_error_msg
#undef error_msg
#undef v
#undef r_start
#undef r
#undef starting_lists
#undef Geometries
#undef Volumes
#undef intersection_time_table
#undef geometry_list_index
#undef previous_master_index
#undef global_master_element
#undef starting_volume_warning
#undef finally_verbal
#undef trace_verbal
#undef list_verbal
#undef verbal
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Banana_monitor' [15]. */
#define mccompcurname  Banana_monitor
#define mccompcurtype  Monitor_nD
#define mccompcurindex 15
#define user1 mccBanana_monitor_user1
#define user2 mccBanana_monitor_user2
#define user3 mccBanana_monitor_user3
#define DEFS mccBanana_monitor_DEFS
#define Vars mccBanana_monitor_Vars
#define detector mccBanana_monitor_detector
#define offdata mccBanana_monitor_offdata
#define xwidth mccBanana_monitor_xwidth
#define yheight mccBanana_monitor_yheight
#define zdepth mccBanana_monitor_zdepth
#define xmin mccBanana_monitor_xmin
#define xmax mccBanana_monitor_xmax
#define ymin mccBanana_monitor_ymin
#define ymax mccBanana_monitor_ymax
#define zmin mccBanana_monitor_zmin
#define zmax mccBanana_monitor_zmax
#define bins mccBanana_monitor_bins
#define min mccBanana_monitor_min
#define max mccBanana_monitor_max
#define restore_neutron mccBanana_monitor_restore_neutron
#define radius mccBanana_monitor_radius
#define options mccBanana_monitor_options
#define filename mccBanana_monitor_filename
#define geometry mccBanana_monitor_geometry
#define username1 mccBanana_monitor_username1
#define username2 mccBanana_monitor_username2
#define username3 mccBanana_monitor_username3
#define nowritefile mccBanana_monitor_nowritefile
#line 223 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
  MonitornD_Defines_type DEFS;
  MonitornD_Variables_type Vars;
  MCDETECTOR detector;
  off_struct offdata;
#line 12327 "Union_manual_example.c"
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

Coords mcposaAl_inc, mcposrAl_inc;
Rotation mcrotaAl_inc, mcrotrAl_inc;
Coords mcposaAl_powder, mcposrAl_powder;
Rotation mcrotaAl_powder, mcrotrAl_powder;
Coords mcposaAl, mcposrAl;
Rotation mcrotaAl, mcrotrAl;
Coords mcposaFe_inc, mcposrFe_inc;
Rotation mcrotaFe_inc, mcrotrFe_inc;
Coords mcposaFe_powder, mcposrFe_powder;
Rotation mcrotaFe_powder, mcrotrFe_powder;
Coords mcposaFe, mcposrFe;
Rotation mcrotaFe, mcrotrFe;
Coords mcposaProgress, mcposrProgress;
Rotation mcrotaProgress, mcrotrProgress;
Coords mcposaSource, mcposrSource;
Rotation mcrotaSource, mcrotrSource;
Coords mcposaGuide, mcposrGuide;
Rotation mcrotaGuide, mcrotrGuide;
Coords mcposasample_position, mcposrsample_position;
Rotation mcrotasample_position, mcrotrsample_position;
Coords mcposaHexagonal_container_1, mcposrHexagonal_container_1;
Rotation mcrotaHexagonal_container_1, mcrotrHexagonal_container_1;
Coords mcposaHexagonal_container_2, mcposrHexagonal_container_2;
Rotation mcrotaHexagonal_container_2, mcrotrHexagonal_container_2;
Coords mcposaSample, mcposrSample;
Rotation mcrotaSample, mcrotrSample;
Coords mcposaMaster, mcposrMaster;
Rotation mcrotaMaster, mcrotrMaster;
Coords mcposaBanana_monitor, mcposrBanana_monitor;
Rotation mcrotaBanana_monitor, mcrotrBanana_monitor;

MCNUM mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz, mcnt, mcnsx, mcnsy, mcnsz, mcnp;

/* end declare */

void mcinit(void) {
#define mccompcurname  Union_manual_example
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaUnion_manual_example coords_set(0,0,0)
#undef mcposaUnion_manual_example
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
    /* Component Al_inc. */
  /* Setting parameters for component Al_inc. */
  SIG_MESSAGE("Al_inc (Init:SetPar)");
#line 34 "Union_manual_example.instr"
  mccAl_inc_sigma = 4 * 0.0082;
#line 34 "Union_manual_example.instr"
  mccAl_inc_packing_factor = 1;
#line 34 "Union_manual_example.instr"
  mccAl_inc_unit_cell_volume = 66.4;
#line 52 "Union_manual_example.instr"
  mccAl_inc_interact_fraction = -1;
#line 12425 "Union_manual_example.c"

  SIG_MESSAGE("Al_inc (Init:Place/Rotate)");
  rot_set_rotation(mcrotaAl_inc,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12432 "Union_manual_example.c"
  rot_copy(mcrotrAl_inc, mcrotaAl_inc);
  mcposaAl_inc = coords_set(
#line 35 "Union_manual_example.instr"
    0,
#line 35 "Union_manual_example.instr"
    0,
#line 35 "Union_manual_example.instr"
    0);
#line 12441 "Union_manual_example.c"
  mctc1 = coords_neg(mcposaAl_inc);
  mcposrAl_inc = rot_apply(mcrotaAl_inc, mctc1);
  mcDEBUG_COMPONENT("Al_inc", mcposaAl_inc, mcrotaAl_inc)
  mccomp_posa[1] = mcposaAl_inc;
  mccomp_posr[1] = mcposrAl_inc;
  mcNCounter[1]  = mcPCounter[1] = mcP2Counter[1] = 0;
  mcAbsorbProp[1]= 0;
    /* Component Al_powder. */
  /* Setting parameters for component Al_powder. */
  SIG_MESSAGE("Al_powder (Init:SetPar)");
#line 37 "Union_manual_example.instr"
  if("Al.laz") strncpy(mccAl_powder_reflections, "Al.laz" ? "Al.laz" : "", 16384); else mccAl_powder_reflections[0]='\0';
#line 49 "Union_manual_example.instr"
  mccAl_powder_packing_factor = 1;
#line 49 "Union_manual_example.instr"
  mccAl_powder_Vc = 0;
#line 49 "Union_manual_example.instr"
  mccAl_powder_delta_d_d = 0;
#line 49 "Union_manual_example.instr"
  mccAl_powder_DW = 0;
#line 49 "Union_manual_example.instr"
  mccAl_powder_nb_atoms = 1;
#line 37 "Union_manual_example.instr"
  mccAl_powder_d_phi = 20;
#line 49 "Union_manual_example.instr"
  mccAl_powder_density = 0;
#line 49 "Union_manual_example.instr"
  mccAl_powder_weight = 0;
#line 49 "Union_manual_example.instr"
  mccAl_powder_barns = 1;
#line 49 "Union_manual_example.instr"
  mccAl_powder_Strain = 0;
#line 49 "Union_manual_example.instr"
  mccAl_powder_interact_fraction = -1;
#line 12476 "Union_manual_example.c"

  SIG_MESSAGE("Al_powder (Init:Place/Rotate)");
  rot_set_rotation(mcrotaAl_powder,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12483 "Union_manual_example.c"
  rot_transpose(mcrotaAl_inc, mctr1);
  rot_mul(mcrotaAl_powder, mctr1, mcrotrAl_powder);
  mcposaAl_powder = coords_set(
#line 38 "Union_manual_example.instr"
    0,
#line 38 "Union_manual_example.instr"
    0,
#line 38 "Union_manual_example.instr"
    0);
#line 12493 "Union_manual_example.c"
  mctc1 = coords_sub(mcposaAl_inc, mcposaAl_powder);
  mcposrAl_powder = rot_apply(mcrotaAl_powder, mctc1);
  mcDEBUG_COMPONENT("Al_powder", mcposaAl_powder, mcrotaAl_powder)
  mccomp_posa[2] = mcposaAl_powder;
  mccomp_posr[2] = mcposrAl_powder;
  mcNCounter[2]  = mcPCounter[2] = mcP2Counter[2] = 0;
  mcAbsorbProp[2]= 0;
    /* Component Al. */
  /* Setting parameters for component Al. */
  SIG_MESSAGE("Al (Init:SetPar)");
#line 40 "Union_manual_example.instr"
  if("Al_inc,Al_powder") strncpy(mccAl_process_string, "Al_inc,Al_powder" ? "Al_inc,Al_powder" : "", 16384); else mccAl_process_string[0]='\0';
#line 40 "Union_manual_example.instr"
  mccAl_my_absorption = 100 * 4 * 0.231 / 66.4;
#line 50 "Union_manual_example.instr"
  mccAl_absorber = 0;
#line 12510 "Union_manual_example.c"

  SIG_MESSAGE("Al (Init:Place/Rotate)");
  rot_set_rotation(mcrotaAl,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12517 "Union_manual_example.c"
  rot_transpose(mcrotaAl_powder, mctr1);
  rot_mul(mcrotaAl, mctr1, mcrotrAl);
  mcposaAl = coords_set(
#line 41 "Union_manual_example.instr"
    0,
#line 41 "Union_manual_example.instr"
    0,
#line 41 "Union_manual_example.instr"
    0);
#line 12527 "Union_manual_example.c"
  mctc1 = coords_sub(mcposaAl_powder, mcposaAl);
  mcposrAl = rot_apply(mcrotaAl, mctc1);
  mcDEBUG_COMPONENT("Al", mcposaAl, mcrotaAl)
  mccomp_posa[3] = mcposaAl;
  mccomp_posr[3] = mcposrAl;
  mcNCounter[3]  = mcPCounter[3] = mcP2Counter[3] = 0;
  mcAbsorbProp[3]= 0;
    /* Component Fe_inc. */
  /* Setting parameters for component Fe_inc. */
  SIG_MESSAGE("Fe_inc (Init:SetPar)");
#line 43 "Union_manual_example.instr"
  mccFe_inc_sigma = 2 * 0.4;
#line 43 "Union_manual_example.instr"
  mccFe_inc_packing_factor = 1;
#line 43 "Union_manual_example.instr"
  mccFe_inc_unit_cell_volume = 24.04;
#line 52 "Union_manual_example.instr"
  mccFe_inc_interact_fraction = -1;
#line 12546 "Union_manual_example.c"

  SIG_MESSAGE("Fe_inc (Init:Place/Rotate)");
  rot_set_rotation(mcrotaFe_inc,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12553 "Union_manual_example.c"
  rot_transpose(mcrotaAl, mctr1);
  rot_mul(mcrotaFe_inc, mctr1, mcrotrFe_inc);
  mcposaFe_inc = coords_set(
#line 44 "Union_manual_example.instr"
    0,
#line 44 "Union_manual_example.instr"
    0,
#line 44 "Union_manual_example.instr"
    0);
#line 12563 "Union_manual_example.c"
  mctc1 = coords_sub(mcposaAl, mcposaFe_inc);
  mcposrFe_inc = rot_apply(mcrotaFe_inc, mctc1);
  mcDEBUG_COMPONENT("Fe_inc", mcposaFe_inc, mcrotaFe_inc)
  mccomp_posa[4] = mcposaFe_inc;
  mccomp_posr[4] = mcposrFe_inc;
  mcNCounter[4]  = mcPCounter[4] = mcP2Counter[4] = 0;
  mcAbsorbProp[4]= 0;
    /* Component Fe_powder. */
  /* Setting parameters for component Fe_powder. */
  SIG_MESSAGE("Fe_powder (Init:SetPar)");
#line 46 "Union_manual_example.instr"
  if("Fe.laz") strncpy(mccFe_powder_reflections, "Fe.laz" ? "Fe.laz" : "", 16384); else mccFe_powder_reflections[0]='\0';
#line 49 "Union_manual_example.instr"
  mccFe_powder_packing_factor = 1;
#line 49 "Union_manual_example.instr"
  mccFe_powder_Vc = 0;
#line 49 "Union_manual_example.instr"
  mccFe_powder_delta_d_d = 0;
#line 49 "Union_manual_example.instr"
  mccFe_powder_DW = 0;
#line 49 "Union_manual_example.instr"
  mccFe_powder_nb_atoms = 1;
#line 46 "Union_manual_example.instr"
  mccFe_powder_d_phi = 20;
#line 49 "Union_manual_example.instr"
  mccFe_powder_density = 0;
#line 49 "Union_manual_example.instr"
  mccFe_powder_weight = 0;
#line 49 "Union_manual_example.instr"
  mccFe_powder_barns = 1;
#line 49 "Union_manual_example.instr"
  mccFe_powder_Strain = 0;
#line 49 "Union_manual_example.instr"
  mccFe_powder_interact_fraction = -1;
#line 12598 "Union_manual_example.c"

  SIG_MESSAGE("Fe_powder (Init:Place/Rotate)");
  rot_set_rotation(mcrotaFe_powder,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12605 "Union_manual_example.c"
  rot_transpose(mcrotaFe_inc, mctr1);
  rot_mul(mcrotaFe_powder, mctr1, mcrotrFe_powder);
  mcposaFe_powder = coords_set(
#line 47 "Union_manual_example.instr"
    0,
#line 47 "Union_manual_example.instr"
    0,
#line 47 "Union_manual_example.instr"
    0);
#line 12615 "Union_manual_example.c"
  mctc1 = coords_sub(mcposaFe_inc, mcposaFe_powder);
  mcposrFe_powder = rot_apply(mcrotaFe_powder, mctc1);
  mcDEBUG_COMPONENT("Fe_powder", mcposaFe_powder, mcrotaFe_powder)
  mccomp_posa[5] = mcposaFe_powder;
  mccomp_posr[5] = mcposrFe_powder;
  mcNCounter[5]  = mcPCounter[5] = mcP2Counter[5] = 0;
  mcAbsorbProp[5]= 0;
    /* Component Fe. */
  /* Setting parameters for component Fe. */
  SIG_MESSAGE("Fe (Init:SetPar)");
#line 49 "Union_manual_example.instr"
  if("Fe_inc,Fe_powder") strncpy(mccFe_process_string, "Fe_inc,Fe_powder" ? "Fe_inc,Fe_powder" : "", 16384); else mccFe_process_string[0]='\0';
#line 49 "Union_manual_example.instr"
  mccFe_my_absorption = 100 * 2 * 2.56 / 24.04;
#line 50 "Union_manual_example.instr"
  mccFe_absorber = 0;
#line 12632 "Union_manual_example.c"

  SIG_MESSAGE("Fe (Init:Place/Rotate)");
  rot_set_rotation(mcrotaFe,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12639 "Union_manual_example.c"
  rot_transpose(mcrotaFe_powder, mctr1);
  rot_mul(mcrotaFe, mctr1, mcrotrFe);
  mcposaFe = coords_set(
#line 50 "Union_manual_example.instr"
    0,
#line 50 "Union_manual_example.instr"
    0,
#line 50 "Union_manual_example.instr"
    0);
#line 12649 "Union_manual_example.c"
  mctc1 = coords_sub(mcposaFe_powder, mcposaFe);
  mcposrFe = rot_apply(mcrotaFe, mctc1);
  mcDEBUG_COMPONENT("Fe", mcposaFe, mcrotaFe)
  mccomp_posa[6] = mcposaFe;
  mccomp_posr[6] = mcposrFe;
  mcNCounter[6]  = mcPCounter[6] = mcP2Counter[6] = 0;
  mcAbsorbProp[6]= 0;
    /* Component Progress. */
  /* Setting parameters for component Progress. */
  SIG_MESSAGE("Progress (Init:SetPar)");
#line 39 "Union_manual_example.instr"
  if("NULL") strncpy(mccProgress_profile, "NULL" ? "NULL" : "", 16384); else mccProgress_profile[0]='\0';
#line 39 "Union_manual_example.instr"
  mccProgress_percent = 10;
#line 39 "Union_manual_example.instr"
  mccProgress_flag_save = 0;
#line 39 "Union_manual_example.instr"
  mccProgress_minutes = 0;
#line 12668 "Union_manual_example.c"

  SIG_MESSAGE("Progress (Init:Place/Rotate)");
  rot_set_rotation(mcrotaProgress,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12675 "Union_manual_example.c"
  rot_transpose(mcrotaFe, mctr1);
  rot_mul(mcrotaProgress, mctr1, mcrotrProgress);
  mcposaProgress = coords_set(
#line 53 "Union_manual_example.instr"
    0,
#line 53 "Union_manual_example.instr"
    0,
#line 53 "Union_manual_example.instr"
    0);
#line 12685 "Union_manual_example.c"
  mctc1 = coords_sub(mcposaFe, mcposaProgress);
  mcposrProgress = rot_apply(mcrotaProgress, mctc1);
  mcDEBUG_COMPONENT("Progress", mcposaProgress, mcrotaProgress)
  mccomp_posa[7] = mcposaProgress;
  mccomp_posr[7] = mcposrProgress;
  mcNCounter[7]  = mcPCounter[7] = mcP2Counter[7] = 0;
  mcAbsorbProp[7]= 0;
    /* Component Source. */
  /* Setting parameters for component Source. */
  SIG_MESSAGE("Source (Init:SetPar)");
#line 52 "Union_manual_example.instr"
  mccSource_radius = 0.1;
#line 55 "Union_manual_example.instr"
  mccSource_yheight = .05;
#line 55 "Union_manual_example.instr"
  mccSource_xwidth = .05;
#line 55 "Union_manual_example.instr"
  mccSource_dist = 2;
#line 55 "Union_manual_example.instr"
  mccSource_focus_xw = .05;
#line 55 "Union_manual_example.instr"
  mccSource_focus_yh = .05;
#line 55 "Union_manual_example.instr"
  mccSource_E0 = 15;
#line 55 "Union_manual_example.instr"
  mccSource_dE = 2;
#line 54 "Union_manual_example.instr"
  mccSource_lambda0 = 0;
#line 54 "Union_manual_example.instr"
  mccSource_dlambda = 0;
#line 55 "Union_manual_example.instr"
  mccSource_flux = 1;
#line 55 "Union_manual_example.instr"
  mccSource_gauss = 0;
#line 55 "Union_manual_example.instr"
  mccSource_target_index = + 1;
#line 12722 "Union_manual_example.c"

  SIG_MESSAGE("Source (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12729 "Union_manual_example.c"
  rot_mul(mctr1, mcrotaProgress, mcrotaSource);
  rot_transpose(mcrotaProgress, mctr1);
  rot_mul(mcrotaSource, mctr1, mcrotrSource);
  mctc1 = coords_set(
#line 56 "Union_manual_example.instr"
    0,
#line 56 "Union_manual_example.instr"
    0,
#line 56 "Union_manual_example.instr"
    0);
#line 12740 "Union_manual_example.c"
  rot_transpose(mcrotaProgress, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaSource = coords_add(mcposaProgress, mctc2);
  mctc1 = coords_sub(mcposaProgress, mcposaSource);
  mcposrSource = rot_apply(mcrotaSource, mctc1);
  mcDEBUG_COMPONENT("Source", mcposaSource, mcrotaSource)
  mccomp_posa[8] = mcposaSource;
  mccomp_posr[8] = mcposrSource;
  mcNCounter[8]  = mcPCounter[8] = mcP2Counter[8] = 0;
  mcAbsorbProp[8]= 0;
    /* Component Guide. */
  /* Setting parameters for component Guide. */
  SIG_MESSAGE("Guide (Init:SetPar)");
#line 58 "Union_manual_example.instr"
  mccGuide_w1 = 0.05;
#line 58 "Union_manual_example.instr"
  mccGuide_h1 = 0.05;
#line 58 "Union_manual_example.instr"
  mccGuide_w2 = 0.05;
#line 58 "Union_manual_example.instr"
  mccGuide_h2 = 0.05;
#line 58 "Union_manual_example.instr"
  mccGuide_l = 12;
#line 58 "Union_manual_example.instr"
  mccGuide_R0 = 0.99;
#line 58 "Union_manual_example.instr"
  mccGuide_Qc = 0.0219;
#line 58 "Union_manual_example.instr"
  mccGuide_alpha = 6.07;
#line 58 "Union_manual_example.instr"
  mccGuide_m = 3.0;
#line 58 "Union_manual_example.instr"
  mccGuide_W = 0.003;
#line 114 "Union_manual_example.instr"
  mccGuide_nslit = 1;
#line 114 "Union_manual_example.instr"
  mccGuide_d = 0.0005;
#line 115 "Union_manual_example.instr"
  mccGuide_mleft = -1;
#line 115 "Union_manual_example.instr"
  mccGuide_mright = -1;
#line 115 "Union_manual_example.instr"
  mccGuide_mtop = -1;
#line 115 "Union_manual_example.instr"
  mccGuide_mbottom = -1;
#line 115 "Union_manual_example.instr"
  mccGuide_nhslit = 1;
#line 115 "Union_manual_example.instr"
  mccGuide_G = 0;
#line 116 "Union_manual_example.instr"
  mccGuide_aleft = -1;
#line 116 "Union_manual_example.instr"
  mccGuide_aright = -1;
#line 116 "Union_manual_example.instr"
  mccGuide_atop = -1;
#line 116 "Union_manual_example.instr"
  mccGuide_abottom = -1;
#line 117 "Union_manual_example.instr"
  mccGuide_wavy = 0;
#line 117 "Union_manual_example.instr"
  mccGuide_wavy_z = 0;
#line 117 "Union_manual_example.instr"
  mccGuide_wavy_tb = 0;
#line 117 "Union_manual_example.instr"
  mccGuide_wavy_lr = 0;
#line 118 "Union_manual_example.instr"
  mccGuide_chamfers = 0;
#line 118 "Union_manual_example.instr"
  mccGuide_chamfers_z = 0;
#line 118 "Union_manual_example.instr"
  mccGuide_chamfers_lr = 0;
#line 118 "Union_manual_example.instr"
  mccGuide_chamfers_tb = 0;
#line 118 "Union_manual_example.instr"
  mccGuide_nelements = 1;
#line 119 "Union_manual_example.instr"
  mccGuide_nu = 0;
#line 119 "Union_manual_example.instr"
  mccGuide_phase = 0;
#line 119 "Union_manual_example.instr"
  if("NULL") strncpy(mccGuide_reflect, "NULL" ? "NULL" : "", 16384); else mccGuide_reflect[0]='\0';
#line 12822 "Union_manual_example.c"

  SIG_MESSAGE("Guide (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12829 "Union_manual_example.c"
  rot_mul(mctr1, mcrotaSource, mcrotaGuide);
  rot_transpose(mcrotaSource, mctr1);
  rot_mul(mcrotaGuide, mctr1, mcrotrGuide);
  mctc1 = coords_set(
#line 59 "Union_manual_example.instr"
    0,
#line 59 "Union_manual_example.instr"
    0,
#line 59 "Union_manual_example.instr"
    2);
#line 12840 "Union_manual_example.c"
  rot_transpose(mcrotaSource, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaGuide = coords_add(mcposaSource, mctc2);
  mctc1 = coords_sub(mcposaSource, mcposaGuide);
  mcposrGuide = rot_apply(mcrotaGuide, mctc1);
  mcDEBUG_COMPONENT("Guide", mcposaGuide, mcrotaGuide)
  mccomp_posa[9] = mcposaGuide;
  mccomp_posr[9] = mcposrGuide;
  mcNCounter[9]  = mcPCounter[9] = mcP2Counter[9] = 0;
  mcAbsorbProp[9]= 0;
    /* Component sample_position. */
  /* Setting parameters for component sample_position. */
  SIG_MESSAGE("sample_position (Init:SetPar)");

  SIG_MESSAGE("sample_position (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12860 "Union_manual_example.c"
  rot_mul(mctr1, mcrotaGuide, mcrotasample_position);
  rot_transpose(mcrotaGuide, mctr1);
  rot_mul(mcrotasample_position, mctr1, mcrotrsample_position);
  mctc1 = coords_set(
#line 62 "Union_manual_example.instr"
    0,
#line 62 "Union_manual_example.instr"
    0,
#line 62 "Union_manual_example.instr"
    12.5);
#line 12871 "Union_manual_example.c"
  rot_transpose(mcrotaGuide, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposasample_position = coords_add(mcposaGuide, mctc2);
  mctc1 = coords_sub(mcposaGuide, mcposasample_position);
  mcposrsample_position = rot_apply(mcrotasample_position, mctc1);
  mcDEBUG_COMPONENT("sample_position", mcposasample_position, mcrotasample_position)
  mccomp_posa[10] = mcposasample_position;
  mccomp_posr[10] = mcposrsample_position;
  mcNCounter[10]  = mcPCounter[10] = mcP2Counter[10] = 0;
  mcAbsorbProp[10]= 0;
    /* Component Hexagonal_container_1. */
  /* Setting parameters for component Hexagonal_container_1. */
  SIG_MESSAGE("Hexagonal_container_1 (Init:SetPar)");
#line 64 "Union_manual_example.instr"
  if("Al") strncpy(mccHexagonal_container_1_material_string, "Al" ? "Al" : "", 16384); else mccHexagonal_container_1_material_string[0]='\0';
#line 64 "Union_manual_example.instr"
  mccHexagonal_container_1_priority = 1;
#line 64 "Union_manual_example.instr"
  mccHexagonal_container_1_xwidth = 0.01;
#line 64 "Union_manual_example.instr"
  mccHexagonal_container_1_yheight = 0.04;
#line 64 "Union_manual_example.instr"
  mccHexagonal_container_1_zdepth = 0.01 * cos ( 30 * DEG2RAD );
#line 64 "Union_manual_example.instr"
  mccHexagonal_container_1_xwidth2 = 0.01 + 2 * 0.01 * sin ( 30 * DEG2RAD );
#line 70 "Union_manual_example.instr"
  mccHexagonal_container_1_yheight2 = -1;
#line 70 "Union_manual_example.instr"
  mccHexagonal_container_1_visualize = 1;
#line 70 "Union_manual_example.instr"
  mccHexagonal_container_1_target_index = 0;
#line 70 "Union_manual_example.instr"
  mccHexagonal_container_1_target_x = 0;
#line 70 "Union_manual_example.instr"
  mccHexagonal_container_1_target_y = 0;
#line 70 "Union_manual_example.instr"
  mccHexagonal_container_1_target_z = 0;
#line 70 "Union_manual_example.instr"
  mccHexagonal_container_1_focus_aw = 0;
#line 70 "Union_manual_example.instr"
  mccHexagonal_container_1_focus_ah = 0;
#line 70 "Union_manual_example.instr"
  mccHexagonal_container_1_focus_xw = 0;
#line 70 "Union_manual_example.instr"
  mccHexagonal_container_1_focus_xh = 0;
#line 70 "Union_manual_example.instr"
  mccHexagonal_container_1_focus_r = 0;
#line 64 "Union_manual_example.instr"
  mccHexagonal_container_1_p_interact = 0.3;
#line 70 "Union_manual_example.instr"
  if(0) strncpy(mccHexagonal_container_1_mask_string, 0 ? 0 : "", 16384); else mccHexagonal_container_1_mask_string[0]='\0';
#line 70 "Union_manual_example.instr"
  if(0) strncpy(mccHexagonal_container_1_mask_setting, 0 ? 0 : "", 16384); else mccHexagonal_container_1_mask_setting[0]='\0';
#line 70 "Union_manual_example.instr"
  mccHexagonal_container_1_number_of_activations = 1;
#line 12927 "Union_manual_example.c"

  SIG_MESSAGE("Hexagonal_container_1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 66 "Union_manual_example.instr"
    (0)*DEG2RAD,
#line 66 "Union_manual_example.instr"
    (0)*DEG2RAD,
#line 66 "Union_manual_example.instr"
    (0)*DEG2RAD);
#line 12937 "Union_manual_example.c"
  rot_mul(mctr1, mcrotasample_position, mcrotaHexagonal_container_1);
  rot_transpose(mcrotasample_position, mctr1);
  rot_mul(mcrotaHexagonal_container_1, mctr1, mcrotrHexagonal_container_1);
  mctc1 = coords_set(
#line 65 "Union_manual_example.instr"
    0,
#line 65 "Union_manual_example.instr"
    0,
#line 65 "Union_manual_example.instr"
    -0.5 * 0.01 * cos ( 30 * DEG2RAD ) -0.00001);
#line 12948 "Union_manual_example.c"
  rot_transpose(mcrotasample_position, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaHexagonal_container_1 = coords_add(mcposasample_position, mctc2);
  mctc1 = coords_sub(mcposasample_position, mcposaHexagonal_container_1);
  mcposrHexagonal_container_1 = rot_apply(mcrotaHexagonal_container_1, mctc1);
  mcDEBUG_COMPONENT("Hexagonal_container_1", mcposaHexagonal_container_1, mcrotaHexagonal_container_1)
  mccomp_posa[11] = mcposaHexagonal_container_1;
  mccomp_posr[11] = mcposrHexagonal_container_1;
  mcNCounter[11]  = mcPCounter[11] = mcP2Counter[11] = 0;
  mcAbsorbProp[11]= 0;
    /* Component Hexagonal_container_2. */
  /* Setting parameters for component Hexagonal_container_2. */
  SIG_MESSAGE("Hexagonal_container_2 (Init:SetPar)");
#line 68 "Union_manual_example.instr"
  if("Al") strncpy(mccHexagonal_container_2_material_string, "Al" ? "Al" : "", 16384); else mccHexagonal_container_2_material_string[0]='\0';
#line 68 "Union_manual_example.instr"
  mccHexagonal_container_2_priority = 2;
#line 68 "Union_manual_example.instr"
  mccHexagonal_container_2_xwidth = 0.01 + 2 * 0.01 * sin ( 30 * DEG2RAD );
#line 68 "Union_manual_example.instr"
  mccHexagonal_container_2_yheight = 0.04;
#line 68 "Union_manual_example.instr"
  mccHexagonal_container_2_zdepth = 0.01 * cos ( 30 * DEG2RAD );
#line 68 "Union_manual_example.instr"
  mccHexagonal_container_2_xwidth2 = 0.01;
#line 70 "Union_manual_example.instr"
  mccHexagonal_container_2_yheight2 = -1;
#line 70 "Union_manual_example.instr"
  mccHexagonal_container_2_visualize = 1;
#line 70 "Union_manual_example.instr"
  mccHexagonal_container_2_target_index = 0;
#line 70 "Union_manual_example.instr"
  mccHexagonal_container_2_target_x = 0;
#line 70 "Union_manual_example.instr"
  mccHexagonal_container_2_target_y = 0;
#line 70 "Union_manual_example.instr"
  mccHexagonal_container_2_target_z = 0;
#line 70 "Union_manual_example.instr"
  mccHexagonal_container_2_focus_aw = 0;
#line 70 "Union_manual_example.instr"
  mccHexagonal_container_2_focus_ah = 0;
#line 70 "Union_manual_example.instr"
  mccHexagonal_container_2_focus_xw = 0;
#line 70 "Union_manual_example.instr"
  mccHexagonal_container_2_focus_xh = 0;
#line 70 "Union_manual_example.instr"
  mccHexagonal_container_2_focus_r = 0;
#line 68 "Union_manual_example.instr"
  mccHexagonal_container_2_p_interact = 0.3;
#line 70 "Union_manual_example.instr"
  if(0) strncpy(mccHexagonal_container_2_mask_string, 0 ? 0 : "", 16384); else mccHexagonal_container_2_mask_string[0]='\0';
#line 70 "Union_manual_example.instr"
  if(0) strncpy(mccHexagonal_container_2_mask_setting, 0 ? 0 : "", 16384); else mccHexagonal_container_2_mask_setting[0]='\0';
#line 70 "Union_manual_example.instr"
  mccHexagonal_container_2_number_of_activations = 1;
#line 13004 "Union_manual_example.c"

  SIG_MESSAGE("Hexagonal_container_2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 70 "Union_manual_example.instr"
    (0)*DEG2RAD,
#line 70 "Union_manual_example.instr"
    (0)*DEG2RAD,
#line 70 "Union_manual_example.instr"
    (0)*DEG2RAD);
#line 13014 "Union_manual_example.c"
  rot_mul(mctr1, mcrotasample_position, mcrotaHexagonal_container_2);
  rot_transpose(mcrotaHexagonal_container_1, mctr1);
  rot_mul(mcrotaHexagonal_container_2, mctr1, mcrotrHexagonal_container_2);
  mctc1 = coords_set(
#line 69 "Union_manual_example.instr"
    0,
#line 69 "Union_manual_example.instr"
    0,
#line 69 "Union_manual_example.instr"
    0.5 * 0.01 * cos ( 30 * DEG2RAD ) + 0.00001);
#line 13025 "Union_manual_example.c"
  rot_transpose(mcrotasample_position, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaHexagonal_container_2 = coords_add(mcposasample_position, mctc2);
  mctc1 = coords_sub(mcposaHexagonal_container_1, mcposaHexagonal_container_2);
  mcposrHexagonal_container_2 = rot_apply(mcrotaHexagonal_container_2, mctc1);
  mcDEBUG_COMPONENT("Hexagonal_container_2", mcposaHexagonal_container_2, mcrotaHexagonal_container_2)
  mccomp_posa[12] = mcposaHexagonal_container_2;
  mccomp_posr[12] = mcposrHexagonal_container_2;
  mcNCounter[12]  = mcPCounter[12] = mcP2Counter[12] = 0;
  mcAbsorbProp[12]= 0;
    /* Component Sample. */
  /* Setting parameters for component Sample. */
  SIG_MESSAGE("Sample (Init:SetPar)");
#line 72 "Union_manual_example.instr"
  if("Fe") strncpy(mccSample_material_string, "Fe" ? "Fe" : "", 16384); else mccSample_material_string[0]='\0';
#line 72 "Union_manual_example.instr"
  mccSample_priority = 3;
#line 72 "Union_manual_example.instr"
  mccSample_radius = 0.008;
#line 72 "Union_manual_example.instr"
  mccSample_yheight = 0.36;
#line 66 "Union_manual_example.instr"
  mccSample_visualize = 1;
#line 66 "Union_manual_example.instr"
  mccSample_target_index = 0;
#line 66 "Union_manual_example.instr"
  mccSample_target_x = 0;
#line 66 "Union_manual_example.instr"
  mccSample_target_y = 0;
#line 66 "Union_manual_example.instr"
  mccSample_target_z = 0;
#line 66 "Union_manual_example.instr"
  mccSample_focus_aw = 0;
#line 66 "Union_manual_example.instr"
  mccSample_focus_ah = 0;
#line 66 "Union_manual_example.instr"
  mccSample_focus_xw = 0;
#line 66 "Union_manual_example.instr"
  mccSample_focus_xh = 0;
#line 66 "Union_manual_example.instr"
  mccSample_focus_r = 0;
#line 72 "Union_manual_example.instr"
  mccSample_p_interact = 0.5;
#line 66 "Union_manual_example.instr"
  if(0) strncpy(mccSample_mask_string, 0 ? 0 : "", 16384); else mccSample_mask_string[0]='\0';
#line 66 "Union_manual_example.instr"
  if(0) strncpy(mccSample_mask_setting, 0 ? 0 : "", 16384); else mccSample_mask_setting[0]='\0';
#line 66 "Union_manual_example.instr"
  mccSample_number_of_activations = 1;
#line 13075 "Union_manual_example.c"

  SIG_MESSAGE("Sample (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 74 "Union_manual_example.instr"
    (0)*DEG2RAD,
#line 74 "Union_manual_example.instr"
    (0)*DEG2RAD,
#line 74 "Union_manual_example.instr"
    (0)*DEG2RAD);
#line 13085 "Union_manual_example.c"
  rot_mul(mctr1, mcrotasample_position, mcrotaSample);
  rot_transpose(mcrotaHexagonal_container_2, mctr1);
  rot_mul(mcrotaSample, mctr1, mcrotrSample);
  mctc1 = coords_set(
#line 73 "Union_manual_example.instr"
    0,
#line 73 "Union_manual_example.instr"
    0,
#line 73 "Union_manual_example.instr"
    0);
#line 13096 "Union_manual_example.c"
  rot_transpose(mcrotasample_position, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaSample = coords_add(mcposasample_position, mctc2);
  mctc1 = coords_sub(mcposaHexagonal_container_2, mcposaSample);
  mcposrSample = rot_apply(mcrotaSample, mctc1);
  mcDEBUG_COMPONENT("Sample", mcposaSample, mcrotaSample)
  mccomp_posa[13] = mcposaSample;
  mccomp_posr[13] = mcposrSample;
  mcNCounter[13]  = mcPCounter[13] = mcP2Counter[13] = 0;
  mcAbsorbProp[13]= 0;
    /* Component Master. */
  /* Setting parameters for component Master. */
  SIG_MESSAGE("Master (Init:SetPar)");
#line 44 "Union_manual_example.instr"
  mccMaster_allow_inside_start = 0;
#line 44 "Union_manual_example.instr"
  mccMaster_history_limit = 300000;
#line 13114 "Union_manual_example.c"

  SIG_MESSAGE("Master (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 13121 "Union_manual_example.c"
  rot_mul(mctr1, mcrotasample_position, mcrotaMaster);
  rot_transpose(mcrotaSample, mctr1);
  rot_mul(mcrotaMaster, mctr1, mcrotrMaster);
  mctc1 = coords_set(
#line 77 "Union_manual_example.instr"
    0,
#line 77 "Union_manual_example.instr"
    0,
#line 77 "Union_manual_example.instr"
    0);
#line 13132 "Union_manual_example.c"
  rot_transpose(mcrotasample_position, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMaster = coords_add(mcposasample_position, mctc2);
  mctc1 = coords_sub(mcposaSample, mcposaMaster);
  mcposrMaster = rot_apply(mcrotaMaster, mctc1);
  mcDEBUG_COMPONENT("Master", mcposaMaster, mcrotaMaster)
  mccomp_posa[14] = mcposaMaster;
  mccomp_posr[14] = mcposrMaster;
  mcNCounter[14]  = mcPCounter[14] = mcP2Counter[14] = 0;
  mcAbsorbProp[14]= 0;
    /* Component Banana_monitor. */
  /* Setting parameters for component Banana_monitor. */
  SIG_MESSAGE("Banana_monitor (Init:SetPar)");
#line 201 "Union_manual_example.instr"
  mccBanana_monitor_xwidth = 0;
#line 79 "Union_manual_example.instr"
  mccBanana_monitor_yheight = 0.1;
#line 201 "Union_manual_example.instr"
  mccBanana_monitor_zdepth = 0;
#line 202 "Union_manual_example.instr"
  mccBanana_monitor_xmin = 0;
#line 202 "Union_manual_example.instr"
  mccBanana_monitor_xmax = 0;
#line 202 "Union_manual_example.instr"
  mccBanana_monitor_ymin = 0;
#line 202 "Union_manual_example.instr"
  mccBanana_monitor_ymax = 0;
#line 202 "Union_manual_example.instr"
  mccBanana_monitor_zmin = 0;
#line 202 "Union_manual_example.instr"
  mccBanana_monitor_zmax = 0;
#line 203 "Union_manual_example.instr"
  mccBanana_monitor_bins = 0;
#line 203 "Union_manual_example.instr"
  mccBanana_monitor_min = -1e40;
#line 203 "Union_manual_example.instr"
  mccBanana_monitor_max = 1e40;
#line 203 "Union_manual_example.instr"
  mccBanana_monitor_restore_neutron = 0;
#line 79 "Union_manual_example.instr"
  mccBanana_monitor_radius = 1;
#line 79 "Union_manual_example.instr"
  if("banana, theta limits=[20,170], bins=200") strncpy(mccBanana_monitor_options, "banana, theta limits=[20,170], bins=200" ? "banana, theta limits=[20,170], bins=200" : "", 16384); else mccBanana_monitor_options[0]='\0';
#line 79 "Union_manual_example.instr"
  if("banana.dat") strncpy(mccBanana_monitor_filename, "banana.dat" ? "banana.dat" : "", 16384); else mccBanana_monitor_filename[0]='\0';
#line 204 "Union_manual_example.instr"
  if("NULL") strncpy(mccBanana_monitor_geometry, "NULL" ? "NULL" : "", 16384); else mccBanana_monitor_geometry[0]='\0';
#line 205 "Union_manual_example.instr"
  if("NULL") strncpy(mccBanana_monitor_username1, "NULL" ? "NULL" : "", 16384); else mccBanana_monitor_username1[0]='\0';
#line 205 "Union_manual_example.instr"
  if("NULL") strncpy(mccBanana_monitor_username2, "NULL" ? "NULL" : "", 16384); else mccBanana_monitor_username2[0]='\0';
#line 205 "Union_manual_example.instr"
  if("NULL") strncpy(mccBanana_monitor_username3, "NULL" ? "NULL" : "", 16384); else mccBanana_monitor_username3[0]='\0';
#line 206 "Union_manual_example.instr"
  mccBanana_monitor_nowritefile = 0;
#line 13188 "Union_manual_example.c"

  SIG_MESSAGE("Banana_monitor (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 13195 "Union_manual_example.c"
  rot_mul(mctr1, mcrotasample_position, mcrotaBanana_monitor);
  rot_transpose(mcrotaMaster, mctr1);
  rot_mul(mcrotaBanana_monitor, mctr1, mcrotrBanana_monitor);
  mctc1 = coords_set(
#line 80 "Union_manual_example.instr"
    0,
#line 80 "Union_manual_example.instr"
    0,
#line 80 "Union_manual_example.instr"
    0);
#line 13206 "Union_manual_example.c"
  rot_transpose(mcrotasample_position, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaBanana_monitor = coords_add(mcposasample_position, mctc2);
  mctc1 = coords_sub(mcposaMaster, mcposaBanana_monitor);
  mcposrBanana_monitor = rot_apply(mcrotaBanana_monitor, mctc1);
  mcDEBUG_COMPONENT("Banana_monitor", mcposaBanana_monitor, mcrotaBanana_monitor)
  mccomp_posa[15] = mcposaBanana_monitor;
  mccomp_posr[15] = mcposrBanana_monitor;
  mcNCounter[15]  = mcPCounter[15] = mcP2Counter[15] = 0;
  mcAbsorbProp[15]= 0;
  /* Component initializations. */
  /* Initializations for component Al_inc. */
  SIG_MESSAGE("Al_inc (Init)");
#define mccompcurname  Al_inc
#define mccompcurtype  Incoherent_process
#define mccompcurindex 1
#define This_process mccAl_inc_This_process
#define Incoherent_storage mccAl_inc_Incoherent_storage
#define effective_my_scattering mccAl_inc_effective_my_scattering
#define sigma mccAl_inc_sigma
#define packing_factor mccAl_inc_packing_factor
#define unit_cell_volume mccAl_inc_unit_cell_volume
#define interact_fraction mccAl_inc_interact_fraction
#line 120 "/usr/share/mcstas/2.5/contrib/union/Incoherent_process.comp"
{
  // Initialize done in the component
  effective_my_scattering = ((packing_factor/unit_cell_volume) * 100 * sigma);

  // Need to specify if this process is isotropic
  This_process.non_isotropic_rot_index = -1; // Yes (powder)
  //This_process.non_isotropic_rot_index =  1;  // No (single crystal)

  // Packing the data into a structure that is transported to the main component
  sprintf(This_process.name,NAME_CURRENT_COMP);
  This_process.process_p_interact = interact_fraction;
  This_process.data_transfer.pointer_to_a_Incoherent_physics_storage_struct = &Incoherent_storage;
  This_process.data_transfer.pointer_to_a_Incoherent_physics_storage_struct->my_scattering = effective_my_scattering;
  This_process.probability_for_scattering_function = &Incoherent_physics_my;
  This_process.scattering_function = &Incoherent_physics_scattering;

  // This will be the same for all process's, and can thus be moved to an include.
  sprintf(global_process_element.name,NAME_CURRENT_COMP);
  global_process_element.component_index = INDEX_CURRENT_COMP;
  global_process_element.p_scattering_process = &This_process;
  add_element_to_process_list(&global_process_list,global_process_element);
}
#line 13253 "Union_manual_example.c"
#undef interact_fraction
#undef unit_cell_volume
#undef packing_factor
#undef sigma
#undef effective_my_scattering
#undef Incoherent_storage
#undef This_process
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Al_powder. */
  SIG_MESSAGE("Al_powder (Init)");
#define mccompcurname  Al_powder
#define mccompcurtype  Powder_process
#define mccompcurindex 2
#define format mccAl_powder_format
#define This_process mccAl_powder_This_process
#define Powder_storage mccAl_powder_Powder_storage
#define effective_my_scattering mccAl_powder_effective_my_scattering
#define d_phi mccAl_powder_d_phi
#define line_info mccAl_powder_line_info
#define columns mccAl_powder_columns
#define reflections mccAl_powder_reflections
#define packing_factor mccAl_powder_packing_factor
#define Vc mccAl_powder_Vc
#define delta_d_d mccAl_powder_delta_d_d
#define DW mccAl_powder_DW
#define nb_atoms mccAl_powder_nb_atoms
#define d_phi mccAl_powder_d_phi
#define density mccAl_powder_density
#define weight mccAl_powder_weight
#define barns mccAl_powder_barns
#define Strain mccAl_powder_Strain
#define interact_fraction mccAl_powder_interact_fraction
#line 637 "/usr/share/mcstas/2.5/contrib/union/Powder_process.comp"
{

  // Initialize done in the component
  
  // Copy from PowderN component
  int i=0;
  struct line_data_union *L;
  line_info.Dd       = delta_d_d;
  line_info.DWfactor = DW;
  line_info.V_0      = Vc;
  line_info.rho      = density;
  line_info.at_weight= weight;
  line_info.at_nb    = nb_atoms;
  line_info.sigma_a  = 0; // This inputs are not needed, as absorption is handled elsewhere
  line_info.sigma_i  = 0; // This input is not needed, as incoherent scattering is handled elsewhere
  line_info.flag_barns=barns;
  //line_info.shape    = 0;
  line_info.flag_warning=0;
  line_info.Epsilon  = Strain;
  line_info.radius_i =line_info.xwidth_i=line_info.yheight_i=line_info.zdepth_i=0;
  line_info.v  = 0;
  line_info.Nq = 0;
  //line_info.v_min = FLT_MAX; line_info.v_max = 0;
  line_info.v_min = 10000000000; line_info.v_max = 0;
  line_info.neutron_passed=0;
  line_info.nb_reuses = line_info.nb_refl = line_info.nb_refl_count = 0;
  line_info.xs_compute= line_info.xs_reuse= line_info.xs_calls =0;
  for (i=0; i< 9; i++) line_info.column_order[i] = columns[i];
  strncpy(line_info.compname, NAME_CURRENT_COMP, 256);

  // p_interact handled elsewhere
  //if (p_interact) {
  //  if (p_interact < p_inc) { double tmp=p_interact; p_interact=p_inc; p_inc=tmp; }
  //  p_transmit = 1-p_interact-p_inc;
  //}

  if (reflections && strlen(reflections) && strcmp(reflections, "NULL") && strcmp(reflections, "0")) {
    i = read_line_data_union(reflections, &line_info);
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

  if (line_info.V_0 <= 0)
    fprintf(stderr,"PowderN: %s: density/unit cell volume is NULL (Vc). Unactivating component.\n", NAME_CURRENT_COMP);


  if (line_info.flag_barns) { /* Factor 100 to convert from barns to fm^2 */
    line_info.XsectionFactor = 100;
  } else {
    line_info.XsectionFactor = 1;
  }

  if (line_info.V_0 && i) {
    L = line_info.list;

    line_info.q_v = malloc(line_info.count*sizeof(double));
    line_info.w_v = malloc(line_info.count*sizeof(double));
    line_info.my_s_v2 = malloc(line_info.count*sizeof(double));
    if (!line_info.q_v || !line_info.w_v || !line_info.my_s_v2)
      exit(fprintf(stderr,"PowderN: %s: ERROR allocating memory (init)\n", NAME_CURRENT_COMP));
    for(i=0; i<line_info.count; i++)
    {
      line_info.my_s_v2[i] = 4*PI*PI*PI*packing_factor*(L[i].DWfactor ? L[i].DWfactor : 1)
                 /(line_info.V_0*line_info.V_0*V2K*V2K)
                 *(L[i].j * L[i].F2 / L[i].q)*line_info.XsectionFactor;
      /* Is not yet divided by v^2 */
      /* Squires [3.103] */
      line_info.q_v[i] = L[i].q*K2V;
      line_info.w_v[i] = L[i].w;
    }
  }
  if (line_info.V_0) {
    /* Is not yet divided by v */
    line_info.my_a_v = packing_factor*line_info.sigma_a/line_info.V_0*2200*100;   // Factor 100 to convert from barns to fm^2
    line_info.my_inc = packing_factor*line_info.sigma_i/line_info.V_0*100;   // Factor 100 to convert from barns to fm^2

    printf("PowderN: %s: Vc=%g [Angs] sigma_abs=%g [barn] sigma_inc=%g [barn] reflections=%s\n",
      NAME_CURRENT_COMP, line_info.V_0, line_info.sigma_a, line_info.sigma_i, reflections && strlen(reflections) ? reflections : "NULL");
  }

  //printf("INTIALIZE line_info.v = %f, line_info.v_min = %f, line_info.v_max = %f, line_info.neutron_passed = %f\n",line_info.v,line_info.v_min,line_info.v_max,line_info.neutron_passed);

  Powder_storage.line_info_storage = &line_info;
  Powder_storage.vertical_angular_limit = d_phi;
  
  // Need to specify if this process is isotropic
  This_process.non_isotropic_rot_index = -1; // Yes (powder)
  //This_process.non_isotropic_rot_index =  1;  // No (single crystal)

  // Packing the data into a structure that is transported to the main component
  This_process.data_transfer.pointer_to_a_Powder_physics_storage_struct = &Powder_storage;
  This_process.data_transfer.pointer_to_a_Powder_physics_storage_struct->my_scattering = effective_my_scattering;
  This_process.probability_for_scattering_function = &Powder_physics_my;
  This_process.scattering_function = &Powder_physics_scattering;
    
  // This will be the same for all process's, and can thus be moved to an include.
  This_process.process_p_interact = interact_fraction;
  sprintf(This_process.name,NAME_CURRENT_COMP);
  rot_copy(This_process.rotation_matrix,ROT_A_CURRENT_COMP);
  sprintf(global_process_element.name,NAME_CURRENT_COMP);
  global_process_element.component_index = INDEX_CURRENT_COMP;
  global_process_element.p_scattering_process = &This_process;
  add_element_to_process_list(&global_process_list,global_process_element);
}
#line 13412 "Union_manual_example.c"
#undef interact_fraction
#undef Strain
#undef barns
#undef weight
#undef density
#undef d_phi
#undef nb_atoms
#undef DW
#undef delta_d_d
#undef Vc
#undef packing_factor
#undef reflections
#undef columns
#undef line_info
#undef d_phi
#undef effective_my_scattering
#undef Powder_storage
#undef This_process
#undef format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Al. */
  SIG_MESSAGE("Al (Init)");
#define mccompcurname  Al
#define mccompcurtype  Union_make_material
#define mccompcurindex 3
#define loop_index mccAl_loop_index
#define this_material mccAl_this_material
#define accepted_processes mccAl_accepted_processes
#define global_material_element mccAl_global_material_element
#define process_string mccAl_process_string
#define my_absorption mccAl_my_absorption
#define absorber mccAl_absorber
#line 193 "/usr/share/mcstas/2.5/contrib/union/Union_make_material.comp"
{

  /*
  // Comma test
  printf("Starting comma test on string: %s \n",process_string);
  printf("Number of commas in string: %d \n",count_commas(process_string));
  exit(1);
  */
  

  if (0 == strcmp(NAME_CURRENT_COMP,"vacuum") || 0 == strcmp(NAME_CURRENT_COMP,"Vacuum")) {
    printf("ERROR, a Union material may not be called Vacuum. A vacuum volume may be created by material=\"Vacuum\" in a geometry component.\n");
    exit(1);
  }
  if (0 == strcmp(NAME_CURRENT_COMP,"exit") || 0 == strcmp(NAME_CURRENT_COMP,"Exit")) {
    printf("ERROR, a Union material may not be called Exit. A exit volume may be created by material=\"Exit\" in a geometry component.\n");
    exit(1);
  }
  if (my_absorption < 0) {
    printf("ERROR, Union make material named %s have a negative absorption cross section!.\n",NAME_CURRENT_COMP);
    exit(1);
  }
  
  if (absorber == 0) {
    if (process_string && strlen(process_string) && strcmp(process_string,"NULL") && strcmp(process_string, "0")) {
        manual_linking_function_material(process_string, &global_process_list, &accepted_processes, NAME_CURRENT_COMP);
    } else {
      for (loop_index=0;loop_index<global_process_list.num_elements;loop_index++) {
        // printf("Automatic linking chosen [loop index = %d] with process_string = %s  \n",loop_index,process_string);
        // automatic linking
        // accept a process if index is between current and former index of make_material components
          if (1 == automatic_linking_materials_function(global_process_list.elements[loop_index],global_material_list,INDEX_CURRENT_COMP))
              add_element_to_int_list(&accepted_processes,loop_index);
      }
    }
  }

  this_material.number_of_processes = accepted_processes.num_elements; // Add number of processes
  this_material.is_vacuum = 0; // This material is not vacuum
  
  if (this_material.number_of_processes == 0 && my_absorption == 0) {
    printf("ERROR, the material named %s has no processes assigned and no absorption cross section, making it eqvialent to vacuum. Vacuums are assigned by setting material=\"Vacuum\" in a geometry component.\n",NAME_CURRENT_COMP);
    exit(1);
  }
  
  // add process element to this_material, building an array of processes called p_scattering_array
  if (this_material.number_of_processes > 0) this_material.p_scattering_array = malloc(this_material.number_of_processes * sizeof(struct scattering_process_struct));
  for (loop_index=0;loop_index<accepted_processes.num_elements;loop_index++) {
        this_material.p_scattering_array[loop_index]=*global_process_list.elements[accepted_processes.elements[loop_index]].p_scattering_process;
  }
  
  this_material.my_a = my_absorption;  // add the absorption to this material
  sprintf(this_material.name,NAME_CURRENT_COMP);
  
  // packing the information into the global_material_element, which is then included in the global_material_list.
  sprintf(global_material_element.name,NAME_CURRENT_COMP);
  global_material_element.component_index = INDEX_CURRENT_COMP;
  global_material_element.physics = &this_material; // Would be nicer if this material was a pointer, now we have the (small) data two places
  add_element_to_material_list(&global_material_list,global_material_element);
}
#line 13509 "Union_manual_example.c"
#undef absorber
#undef my_absorption
#undef process_string
#undef global_material_element
#undef accepted_processes
#undef this_material
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Fe_inc. */
  SIG_MESSAGE("Fe_inc (Init)");
#define mccompcurname  Fe_inc
#define mccompcurtype  Incoherent_process
#define mccompcurindex 4
#define This_process mccFe_inc_This_process
#define Incoherent_storage mccFe_inc_Incoherent_storage
#define effective_my_scattering mccFe_inc_effective_my_scattering
#define sigma mccFe_inc_sigma
#define packing_factor mccFe_inc_packing_factor
#define unit_cell_volume mccFe_inc_unit_cell_volume
#define interact_fraction mccFe_inc_interact_fraction
#line 120 "/usr/share/mcstas/2.5/contrib/union/Incoherent_process.comp"
{
  // Initialize done in the component
  effective_my_scattering = ((packing_factor/unit_cell_volume) * 100 * sigma);

  // Need to specify if this process is isotropic
  This_process.non_isotropic_rot_index = -1; // Yes (powder)
  //This_process.non_isotropic_rot_index =  1;  // No (single crystal)

  // Packing the data into a structure that is transported to the main component
  sprintf(This_process.name,NAME_CURRENT_COMP);
  This_process.process_p_interact = interact_fraction;
  This_process.data_transfer.pointer_to_a_Incoherent_physics_storage_struct = &Incoherent_storage;
  This_process.data_transfer.pointer_to_a_Incoherent_physics_storage_struct->my_scattering = effective_my_scattering;
  This_process.probability_for_scattering_function = &Incoherent_physics_my;
  This_process.scattering_function = &Incoherent_physics_scattering;

  // This will be the same for all process's, and can thus be moved to an include.
  sprintf(global_process_element.name,NAME_CURRENT_COMP);
  global_process_element.component_index = INDEX_CURRENT_COMP;
  global_process_element.p_scattering_process = &This_process;
  add_element_to_process_list(&global_process_list,global_process_element);
}
#line 13556 "Union_manual_example.c"
#undef interact_fraction
#undef unit_cell_volume
#undef packing_factor
#undef sigma
#undef effective_my_scattering
#undef Incoherent_storage
#undef This_process
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Fe_powder. */
  SIG_MESSAGE("Fe_powder (Init)");
#define mccompcurname  Fe_powder
#define mccompcurtype  Powder_process
#define mccompcurindex 5
#define format mccFe_powder_format
#define This_process mccFe_powder_This_process
#define Powder_storage mccFe_powder_Powder_storage
#define effective_my_scattering mccFe_powder_effective_my_scattering
#define d_phi mccFe_powder_d_phi
#define line_info mccFe_powder_line_info
#define columns mccFe_powder_columns
#define reflections mccFe_powder_reflections
#define packing_factor mccFe_powder_packing_factor
#define Vc mccFe_powder_Vc
#define delta_d_d mccFe_powder_delta_d_d
#define DW mccFe_powder_DW
#define nb_atoms mccFe_powder_nb_atoms
#define d_phi mccFe_powder_d_phi
#define density mccFe_powder_density
#define weight mccFe_powder_weight
#define barns mccFe_powder_barns
#define Strain mccFe_powder_Strain
#define interact_fraction mccFe_powder_interact_fraction
#line 637 "/usr/share/mcstas/2.5/contrib/union/Powder_process.comp"
{

  // Initialize done in the component
  
  // Copy from PowderN component
  int i=0;
  struct line_data_union *L;
  line_info.Dd       = delta_d_d;
  line_info.DWfactor = DW;
  line_info.V_0      = Vc;
  line_info.rho      = density;
  line_info.at_weight= weight;
  line_info.at_nb    = nb_atoms;
  line_info.sigma_a  = 0; // This inputs are not needed, as absorption is handled elsewhere
  line_info.sigma_i  = 0; // This input is not needed, as incoherent scattering is handled elsewhere
  line_info.flag_barns=barns;
  //line_info.shape    = 0;
  line_info.flag_warning=0;
  line_info.Epsilon  = Strain;
  line_info.radius_i =line_info.xwidth_i=line_info.yheight_i=line_info.zdepth_i=0;
  line_info.v  = 0;
  line_info.Nq = 0;
  //line_info.v_min = FLT_MAX; line_info.v_max = 0;
  line_info.v_min = 10000000000; line_info.v_max = 0;
  line_info.neutron_passed=0;
  line_info.nb_reuses = line_info.nb_refl = line_info.nb_refl_count = 0;
  line_info.xs_compute= line_info.xs_reuse= line_info.xs_calls =0;
  for (i=0; i< 9; i++) line_info.column_order[i] = columns[i];
  strncpy(line_info.compname, NAME_CURRENT_COMP, 256);

  // p_interact handled elsewhere
  //if (p_interact) {
  //  if (p_interact < p_inc) { double tmp=p_interact; p_interact=p_inc; p_inc=tmp; }
  //  p_transmit = 1-p_interact-p_inc;
  //}

  if (reflections && strlen(reflections) && strcmp(reflections, "NULL") && strcmp(reflections, "0")) {
    i = read_line_data_union(reflections, &line_info);
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

  if (line_info.V_0 <= 0)
    fprintf(stderr,"PowderN: %s: density/unit cell volume is NULL (Vc). Unactivating component.\n", NAME_CURRENT_COMP);


  if (line_info.flag_barns) { /* Factor 100 to convert from barns to fm^2 */
    line_info.XsectionFactor = 100;
  } else {
    line_info.XsectionFactor = 1;
  }

  if (line_info.V_0 && i) {
    L = line_info.list;

    line_info.q_v = malloc(line_info.count*sizeof(double));
    line_info.w_v = malloc(line_info.count*sizeof(double));
    line_info.my_s_v2 = malloc(line_info.count*sizeof(double));
    if (!line_info.q_v || !line_info.w_v || !line_info.my_s_v2)
      exit(fprintf(stderr,"PowderN: %s: ERROR allocating memory (init)\n", NAME_CURRENT_COMP));
    for(i=0; i<line_info.count; i++)
    {
      line_info.my_s_v2[i] = 4*PI*PI*PI*packing_factor*(L[i].DWfactor ? L[i].DWfactor : 1)
                 /(line_info.V_0*line_info.V_0*V2K*V2K)
                 *(L[i].j * L[i].F2 / L[i].q)*line_info.XsectionFactor;
      /* Is not yet divided by v^2 */
      /* Squires [3.103] */
      line_info.q_v[i] = L[i].q*K2V;
      line_info.w_v[i] = L[i].w;
    }
  }
  if (line_info.V_0) {
    /* Is not yet divided by v */
    line_info.my_a_v = packing_factor*line_info.sigma_a/line_info.V_0*2200*100;   // Factor 100 to convert from barns to fm^2
    line_info.my_inc = packing_factor*line_info.sigma_i/line_info.V_0*100;   // Factor 100 to convert from barns to fm^2

    printf("PowderN: %s: Vc=%g [Angs] sigma_abs=%g [barn] sigma_inc=%g [barn] reflections=%s\n",
      NAME_CURRENT_COMP, line_info.V_0, line_info.sigma_a, line_info.sigma_i, reflections && strlen(reflections) ? reflections : "NULL");
  }

  //printf("INTIALIZE line_info.v = %f, line_info.v_min = %f, line_info.v_max = %f, line_info.neutron_passed = %f\n",line_info.v,line_info.v_min,line_info.v_max,line_info.neutron_passed);

  Powder_storage.line_info_storage = &line_info;
  Powder_storage.vertical_angular_limit = d_phi;
  
  // Need to specify if this process is isotropic
  This_process.non_isotropic_rot_index = -1; // Yes (powder)
  //This_process.non_isotropic_rot_index =  1;  // No (single crystal)

  // Packing the data into a structure that is transported to the main component
  This_process.data_transfer.pointer_to_a_Powder_physics_storage_struct = &Powder_storage;
  This_process.data_transfer.pointer_to_a_Powder_physics_storage_struct->my_scattering = effective_my_scattering;
  This_process.probability_for_scattering_function = &Powder_physics_my;
  This_process.scattering_function = &Powder_physics_scattering;
    
  // This will be the same for all process's, and can thus be moved to an include.
  This_process.process_p_interact = interact_fraction;
  sprintf(This_process.name,NAME_CURRENT_COMP);
  rot_copy(This_process.rotation_matrix,ROT_A_CURRENT_COMP);
  sprintf(global_process_element.name,NAME_CURRENT_COMP);
  global_process_element.component_index = INDEX_CURRENT_COMP;
  global_process_element.p_scattering_process = &This_process;
  add_element_to_process_list(&global_process_list,global_process_element);
}
#line 13715 "Union_manual_example.c"
#undef interact_fraction
#undef Strain
#undef barns
#undef weight
#undef density
#undef d_phi
#undef nb_atoms
#undef DW
#undef delta_d_d
#undef Vc
#undef packing_factor
#undef reflections
#undef columns
#undef line_info
#undef d_phi
#undef effective_my_scattering
#undef Powder_storage
#undef This_process
#undef format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Fe. */
  SIG_MESSAGE("Fe (Init)");
#define mccompcurname  Fe
#define mccompcurtype  Union_make_material
#define mccompcurindex 6
#define loop_index mccFe_loop_index
#define this_material mccFe_this_material
#define accepted_processes mccFe_accepted_processes
#define global_material_element mccFe_global_material_element
#define process_string mccFe_process_string
#define my_absorption mccFe_my_absorption
#define absorber mccFe_absorber
#line 193 "/usr/share/mcstas/2.5/contrib/union/Union_make_material.comp"
{

  /*
  // Comma test
  printf("Starting comma test on string: %s \n",process_string);
  printf("Number of commas in string: %d \n",count_commas(process_string));
  exit(1);
  */
  

  if (0 == strcmp(NAME_CURRENT_COMP,"vacuum") || 0 == strcmp(NAME_CURRENT_COMP,"Vacuum")) {
    printf("ERROR, a Union material may not be called Vacuum. A vacuum volume may be created by material=\"Vacuum\" in a geometry component.\n");
    exit(1);
  }
  if (0 == strcmp(NAME_CURRENT_COMP,"exit") || 0 == strcmp(NAME_CURRENT_COMP,"Exit")) {
    printf("ERROR, a Union material may not be called Exit. A exit volume may be created by material=\"Exit\" in a geometry component.\n");
    exit(1);
  }
  if (my_absorption < 0) {
    printf("ERROR, Union make material named %s have a negative absorption cross section!.\n",NAME_CURRENT_COMP);
    exit(1);
  }
  
  if (absorber == 0) {
    if (process_string && strlen(process_string) && strcmp(process_string,"NULL") && strcmp(process_string, "0")) {
        manual_linking_function_material(process_string, &global_process_list, &accepted_processes, NAME_CURRENT_COMP);
    } else {
      for (loop_index=0;loop_index<global_process_list.num_elements;loop_index++) {
        // printf("Automatic linking chosen [loop index = %d] with process_string = %s  \n",loop_index,process_string);
        // automatic linking
        // accept a process if index is between current and former index of make_material components
          if (1 == automatic_linking_materials_function(global_process_list.elements[loop_index],global_material_list,INDEX_CURRENT_COMP))
              add_element_to_int_list(&accepted_processes,loop_index);
      }
    }
  }

  this_material.number_of_processes = accepted_processes.num_elements; // Add number of processes
  this_material.is_vacuum = 0; // This material is not vacuum
  
  if (this_material.number_of_processes == 0 && my_absorption == 0) {
    printf("ERROR, the material named %s has no processes assigned and no absorption cross section, making it eqvialent to vacuum. Vacuums are assigned by setting material=\"Vacuum\" in a geometry component.\n",NAME_CURRENT_COMP);
    exit(1);
  }
  
  // add process element to this_material, building an array of processes called p_scattering_array
  if (this_material.number_of_processes > 0) this_material.p_scattering_array = malloc(this_material.number_of_processes * sizeof(struct scattering_process_struct));
  for (loop_index=0;loop_index<accepted_processes.num_elements;loop_index++) {
        this_material.p_scattering_array[loop_index]=*global_process_list.elements[accepted_processes.elements[loop_index]].p_scattering_process;
  }
  
  this_material.my_a = my_absorption;  // add the absorption to this material
  sprintf(this_material.name,NAME_CURRENT_COMP);
  
  // packing the information into the global_material_element, which is then included in the global_material_list.
  sprintf(global_material_element.name,NAME_CURRENT_COMP);
  global_material_element.component_index = INDEX_CURRENT_COMP;
  global_material_element.physics = &this_material; // Would be nicer if this material was a pointer, now we have the (small) data two places
  add_element_to_material_list(&global_material_list,global_material_element);
}
#line 13812 "Union_manual_example.c"
#undef absorber
#undef my_absorption
#undef process_string
#undef global_material_element
#undef accepted_processes
#undef this_material
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Progress. */
  SIG_MESSAGE("Progress (Init)");
#define mccompcurname  Progress
#define mccompcurtype  Progress_bar
#define mccompcurindex 7
#define IntermediateCnts mccProgress_IntermediateCnts
#define StartTime mccProgress_StartTime
#define EndTime mccProgress_EndTime
#define CurrentTime mccProgress_CurrentTime
#define profile mccProgress_profile
#define percent mccProgress_percent
#define flag_save mccProgress_flag_save
#define minutes mccProgress_minutes
#line 57 "/usr/share/mcstas/2.5/misc/Progress_bar.comp"
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
#line 13849 "Union_manual_example.c"
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
#define mccompcurindex 8
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
#line 65 "/usr/share/mcstas/2.5/sources/Source_simple.comp"
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
#line 13943 "Union_manual_example.c"
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

  /* Initializations for component Guide. */
  SIG_MESSAGE("Guide (Init)");
#define mccompcurname  Guide
#define mccompcurtype  Guide_gravity
#define mccompcurindex 9
#define GVars mccGuide_GVars
#define pTable mccGuide_pTable
#define w1 mccGuide_w1
#define h1 mccGuide_h1
#define w2 mccGuide_w2
#define h2 mccGuide_h2
#define l mccGuide_l
#define R0 mccGuide_R0
#define Qc mccGuide_Qc
#define alpha mccGuide_alpha
#define m mccGuide_m
#define W mccGuide_W
#define nslit mccGuide_nslit
#define d mccGuide_d
#define mleft mccGuide_mleft
#define mright mccGuide_mright
#define mtop mccGuide_mtop
#define mbottom mccGuide_mbottom
#define nhslit mccGuide_nhslit
#define G mccGuide_G
#define aleft mccGuide_aleft
#define aright mccGuide_aright
#define atop mccGuide_atop
#define abottom mccGuide_abottom
#define wavy mccGuide_wavy
#define wavy_z mccGuide_wavy_z
#define wavy_tb mccGuide_wavy_tb
#define wavy_lr mccGuide_wavy_lr
#define chamfers mccGuide_chamfers
#define chamfers_z mccGuide_chamfers_z
#define chamfers_lr mccGuide_chamfers_lr
#define chamfers_tb mccGuide_chamfers_tb
#define nelements mccGuide_nelements
#define nu mccGuide_nu
#define phase mccGuide_phase
#define reflect mccGuide_reflect
#line 339 "/usr/share/mcstas/2.5/optics/Guide_gravity.comp"
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
#line 14057 "Union_manual_example.c"
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

  /* Initializations for component sample_position. */
  SIG_MESSAGE("sample_position (Init)");

  /* Initializations for component Hexagonal_container_1. */
  SIG_MESSAGE("Hexagonal_container_1 (Init)");
#define mccompcurname  Hexagonal_container_1
#define mccompcurtype  Union_box
#define mccompcurindex 11
#define loop_index mccHexagonal_container_1_loop_index
#define this_box_volume mccHexagonal_container_1_this_box_volume
#define global_geometry_element mccHexagonal_container_1_global_geometry_element
#define this_box_storage mccHexagonal_container_1_this_box_storage
#define material_string mccHexagonal_container_1_material_string
#define priority mccHexagonal_container_1_priority
#define xwidth mccHexagonal_container_1_xwidth
#define yheight mccHexagonal_container_1_yheight
#define zdepth mccHexagonal_container_1_zdepth
#define xwidth2 mccHexagonal_container_1_xwidth2
#define yheight2 mccHexagonal_container_1_yheight2
#define visualize mccHexagonal_container_1_visualize
#define target_index mccHexagonal_container_1_target_index
#define target_x mccHexagonal_container_1_target_x
#define target_y mccHexagonal_container_1_target_y
#define target_z mccHexagonal_container_1_target_z
#define focus_aw mccHexagonal_container_1_focus_aw
#define focus_ah mccHexagonal_container_1_focus_ah
#define focus_xw mccHexagonal_container_1_focus_xw
#define focus_xh mccHexagonal_container_1_focus_xh
#define focus_r mccHexagonal_container_1_focus_r
#define p_interact mccHexagonal_container_1_p_interact
#define mask_string mccHexagonal_container_1_mask_string
#define mask_setting mccHexagonal_container_1_mask_setting
#define number_of_activations mccHexagonal_container_1_number_of_activations
#line 223 "/usr/share/mcstas/2.5/contrib/union/Union_box.comp"
{
// Initializes the focusing system for this volume including input sanitation.
focus_initialize(&this_box_volume.geometry, POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index), POS_A_CURRENT_COMP, ROT_A_CURRENT_COMP, target_index, target_x, target_y, target_z, focus_aw, focus_ah, focus_xw, focus_xh, focus_r, NAME_CURRENT_COMP);

// Input sanitation for this geometry
if (xwidth <= 0) {
  printf("\nERROR in Union_box named %s, the xwidth is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (yheight <= 0) {
  printf("\nERROR in Union_box named %s, yheight is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (zdepth <= 0) {
  printf("\nERROR in Union_box named %s, zdepth is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (xwidth2 <= 0 && xwidth2 != -1) {
  printf("\nERROR in Union_box named %s, the xwidth2 is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (yheight2 <= 0 && yheight2 != -1) {
  printf("\nERROR in Union_box named %s, yheight2 is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}

// Use sanitation
if (global_material_list.num_elements == 0) {
  printf("\nERROR: Need to define a material using Union_make_material before using a Union geometry component. \n");
  printf("       %s was defined before first use of Union_make_material.\n",NAME_CURRENT_COMP);
  exit(1);
}

this_box_volume.geometry.is_masked_volume = 0;
this_box_volume.geometry.is_exit_volume = 0;
this_box_volume.geometry.is_mask_volume = 0;
// check if the volume is a mask, if it is the material string is irelevant.
if (mask_string && strlen(mask_string) && strcmp(mask_string, "NULL") && strcmp(mask_string, "0")) {
    // A mask volume is used to limit the extend of other volumes, called the masked volumes. These are specified in the mask_string.
    // In order for a ray to enter a masked volume, it needs to be both in the region covered by that volume AND the mask volume.
    // When more than
    this_box_volume.geometry.mask_mode = 1; // Default is mask mode is ALL
    if (mask_setting && strlen(mask_setting) && strcmp(mask_setting, "NULL") && strcmp(mask_setting, "0")) {
        if (strcmp(mask_setting,"ALL") == 0 || strcmp(mask_setting,"All") == 0) this_box_volume.geometry.mask_mode = 1;
        else if (strcmp(mask_setting,"ANY") == 0 || strcmp(mask_setting,"Any") == 0) this_box_volume.geometry.mask_mode = 2;
        else {
            printf("The mask_mode of component %s is set to %s, but must be either ALL or ANY.\n",NAME_CURRENT_COMP,mask_setting);
            exit(1);
        }
    }
    
    for (loop_index=0;loop_index<global_geometry_list.num_elements;loop_index++) {
        // Add mask list
        if (1 == manual_linking_function(global_geometry_list.elements[loop_index].name,mask_string)) {
            add_element_to_int_list(&this_box_volume.geometry.mask_list,global_geometry_list.elements[loop_index].component_index);
            add_element_to_int_list(&global_geometry_list.elements[loop_index].Volume->geometry.masked_by_list,INDEX_CURRENT_COMP);
            global_geometry_list.elements[loop_index].Volume->geometry.is_masked_volume = 1;
            if (this_box_volume.geometry.mask_mode == 2)
                global_geometry_list.elements[loop_index].Volume->geometry.mask_mode = 2;
            if (this_box_volume.geometry.mask_mode == 1) {
                if (global_geometry_list.elements[loop_index].Volume->geometry.is_masked_volume == 1 && global_geometry_list.elements[loop_index].Volume->geometry.mask_mode != 2)
                    // If more than one mask is added to one volume, the ANY mode overwrites the (default) ALL mode.
                    global_geometry_list.elements[loop_index].Volume->geometry.mask_mode = 1;
            }
            
            found_geometries = 1;
        }
    }
    if (found_geometries == 0) {
        printf("The mask_string in geometry: %s did not find any of the specified volumes in the mask_string %s \n",NAME_CURRENT_COMP,mask_string);
        exit(1);
    }
    this_box_volume.p_physics = malloc(sizeof(struct physics_struct));
    this_box_volume.p_physics->is_vacuum = 0; // Makes this volume a vacuum
    this_box_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
    this_box_volume.p_physics->my_a = 0; // Should not be used.
    sprintf(this_box_volume.p_physics->name,"Mask");
    this_box_volume.geometry.is_mask_volume = 1;
    
    
// Read the material input, or if it lacks, use automatic linking.
} else if (material_string && strlen(material_string) && strcmp(material_string, "NULL") && strcmp(material_string, "0")) {
    // A geometry string was given, use it to determine which material
    if (0 == strcmp(material_string,"vacuum") || 0 == strcmp(material_string,"Vacuum")) {
        // One could have a global physics struct for vacuum instead of creating one for each
        this_box_volume.p_physics = malloc(sizeof(struct physics_struct));
        this_box_volume.p_physics->is_vacuum = 1; // Makes this volume a vacuum
        this_box_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
        this_box_volume.p_physics->my_a = 0; // Should not be used.
        sprintf(this_box_volume.p_physics->name,"Vacuum");
    } else if (0 == strcmp(material_string,"exit") || 0 == strcmp(material_string,"Exit")) {
        // One could have a global physics struct for vacuum instead of creating one for each
        this_box_volume.p_physics = malloc(sizeof(struct physics_struct));
        this_box_volume.p_physics->is_vacuum = 1; // Makes this volume a vacuum
        this_box_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
        this_box_volume.p_physics->my_a = 0; // Should not be used.
        this_box_volume.geometry.is_exit_volume = 1;
        sprintf(this_box_volume.p_physics->name,"Exit");
    } else {
        #ifndef MATERIAL_DETECTOR
            printf("Need to define a material before refering to it in a geometry %s.\n",NAME_CURRENT_COMP);
            exit(1);
        #endif
        for (loop_index=0;loop_index<global_material_list.num_elements;loop_index++) {
            if (0 == strcmp(material_string,global_material_list.elements[loop_index].name)) {
                this_box_volume.p_physics = global_material_list.elements[loop_index].physics;
                break;
            }
            if (loop_index == global_material_list.num_elements-1) {
            printf("\n");
            printf("ERROR: The material string \"%s\" in Union geometry \"%s\" did not match a specified material. \n",material_string,NAME_CURRENT_COMP);
            printf("       The materials available at this point (need to be defined before the geometry): \n");
            for (loop_index=0;loop_index<global_material_list.num_elements;loop_index++)
              printf("         %s\n",global_material_list.elements[loop_index].name);
            printf("\n");
            printf("       It is also possible to use one of the defualt materials avaiable: \n");
            printf("           Vacuum (for a Volume without scattering or absorption)\n");
            printf("           Exit (for a Volume where the ray exits the component if it enters)\n");
            printf("           Mask (for a Volume that masks existing volumes specified in the mask_string\n");
            exit(1);
          }
        }
    }
} else {
    // Automatic linking, simply using the last defined material.
    #ifndef MATERIAL_DETECTOR
        printf("Need to define a material before the geometry to use automatic linking %s.\n",NAME_CURRENT_COMP);
        exit(1);
    #endif
    this_box_volume.p_physics = global_material_list.elements[global_material_list.num_elements-1].physics;
}

sprintf(this_box_volume.name,NAME_CURRENT_COMP);
sprintf(this_box_volume.geometry.shape,"box");
this_box_volume.geometry.priority_value = priority;
this_box_volume.geometry.geometry_p_interact = p_interact;
// Currently the coordinates will be in absolute space.
this_box_volume.geometry.center = POS_A_CURRENT_COMP;

this_box_storage.z_depth = zdepth;
this_box_storage.x_width1 = xwidth;
this_box_storage.y_height1 = yheight;

this_box_storage.is_rectangle = 0;
if (xwidth2 < 0 && yheight2 < 0) this_box_storage.is_rectangle = 1;
if (xwidth == xwidth2 && yheight == yheight2) this_box_storage.is_rectangle = 1;

if (xwidth2 < 0) {
    this_box_storage.x_width2 = xwidth;
    xwidth2 = xwidth;
} else this_box_storage.x_width2 = xwidth2;

if (yheight2 < 0) {
    this_box_storage.y_height2 = yheight;
    yheight2 = yheight;
} else this_box_storage.y_height2 = yheight2;


this_box_storage.normal_vectors[0] = coords_set(0,0,1);
this_box_storage.normal_vectors[1] = coords_set(0,0,1);

// for sides with y component = 0
x_component = 2*zdepth/sqrt((xwidth-xwidth2)*(xwidth-xwidth2)+4*zdepth*zdepth);
z_component = (xwidth-xwidth2)/sqrt(4*zdepth*zdepth+(xwidth-xwidth2)*(xwidth-xwidth2));

this_box_storage.normal_vectors[2] = coords_set(x_component,0,z_component);
this_box_storage.normal_vectors[3] = coords_set(-x_component,0,z_component);

// for sides with x component = 0
y_component = 2*zdepth/sqrt((yheight-yheight2)*(yheight-yheight2)+4*zdepth*zdepth);
z_component = (yheight-yheight2)/sqrt(4*zdepth*zdepth+(yheight-yheight2)*(yheight-yheight2));

this_box_storage.normal_vectors[4] = coords_set(0,y_component,z_component);
this_box_storage.normal_vectors[5] = coords_set(0,-y_component,z_component);

this_box_volume.geometry.visualization_on = visualize;

this_box_volume.geometry.geometry_parameters.p_box_storage = &this_box_storage;

// Assign pointers to functions for intersection with the shape, checking if a point is inside the shape
if (this_box_storage.is_rectangle == 1) {
this_box_volume.geometry.intersect_function = &sample_box_intersect_simple;
this_box_volume.geometry.within_function = &r_within_box_simple;
}
else {
this_box_volume.geometry.intersect_function = &sample_box_intersect_advanced;
this_box_volume.geometry.within_function = &r_within_box_advanced;
}

this_box_volume.geometry.shell_points = &box_shell_points;
this_box_volume.geometry.mcdisplay_function = &mcdisplay_box_function;
this_box_volume.geometry.initialize_from_main_function = &initialize_box_geometry_from_main_component;
this_box_volume.geometry.process_rot_allocated = 0;

this_box_volume.geometry.copy_geometry_parameters = &allocate_box_storage_copy;

rot_copy(this_box_volume.geometry.rotation_matrix,ROT_A_CURRENT_COMP); //check how ROT_R_CURRENT_COMP would work
rot_transpose(ROT_A_CURRENT_COMP,this_box_volume.geometry.transpose_rotation_matrix);

// Initialize loggers
this_box_volume.loggers.num_elements = 0;

// packing the information into the global_geometry_element, which is then included in the global_geometry_list.
sprintf(global_geometry_element.name,NAME_CURRENT_COMP);
global_geometry_element.activation_counter = number_of_activations;
global_geometry_element.component_index = INDEX_CURRENT_COMP;
global_geometry_element.Volume = &this_box_volume; // Would be nicer if this m was a pointer, now we have the (small) data two places
add_element_to_geometry_list(&global_geometry_list,global_geometry_element);
}
#line 14341 "Union_manual_example.c"
#undef number_of_activations
#undef mask_setting
#undef mask_string
#undef p_interact
#undef focus_r
#undef focus_xh
#undef focus_xw
#undef focus_ah
#undef focus_aw
#undef target_z
#undef target_y
#undef target_x
#undef target_index
#undef visualize
#undef yheight2
#undef xwidth2
#undef zdepth
#undef yheight
#undef xwidth
#undef priority
#undef material_string
#undef this_box_storage
#undef global_geometry_element
#undef this_box_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Hexagonal_container_2. */
  SIG_MESSAGE("Hexagonal_container_2 (Init)");
#define mccompcurname  Hexagonal_container_2
#define mccompcurtype  Union_box
#define mccompcurindex 12
#define loop_index mccHexagonal_container_2_loop_index
#define this_box_volume mccHexagonal_container_2_this_box_volume
#define global_geometry_element mccHexagonal_container_2_global_geometry_element
#define this_box_storage mccHexagonal_container_2_this_box_storage
#define material_string mccHexagonal_container_2_material_string
#define priority mccHexagonal_container_2_priority
#define xwidth mccHexagonal_container_2_xwidth
#define yheight mccHexagonal_container_2_yheight
#define zdepth mccHexagonal_container_2_zdepth
#define xwidth2 mccHexagonal_container_2_xwidth2
#define yheight2 mccHexagonal_container_2_yheight2
#define visualize mccHexagonal_container_2_visualize
#define target_index mccHexagonal_container_2_target_index
#define target_x mccHexagonal_container_2_target_x
#define target_y mccHexagonal_container_2_target_y
#define target_z mccHexagonal_container_2_target_z
#define focus_aw mccHexagonal_container_2_focus_aw
#define focus_ah mccHexagonal_container_2_focus_ah
#define focus_xw mccHexagonal_container_2_focus_xw
#define focus_xh mccHexagonal_container_2_focus_xh
#define focus_r mccHexagonal_container_2_focus_r
#define p_interact mccHexagonal_container_2_p_interact
#define mask_string mccHexagonal_container_2_mask_string
#define mask_setting mccHexagonal_container_2_mask_setting
#define number_of_activations mccHexagonal_container_2_number_of_activations
#line 223 "/usr/share/mcstas/2.5/contrib/union/Union_box.comp"
{
// Initializes the focusing system for this volume including input sanitation.
focus_initialize(&this_box_volume.geometry, POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index), POS_A_CURRENT_COMP, ROT_A_CURRENT_COMP, target_index, target_x, target_y, target_z, focus_aw, focus_ah, focus_xw, focus_xh, focus_r, NAME_CURRENT_COMP);

// Input sanitation for this geometry
if (xwidth <= 0) {
  printf("\nERROR in Union_box named %s, the xwidth is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (yheight <= 0) {
  printf("\nERROR in Union_box named %s, yheight is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (zdepth <= 0) {
  printf("\nERROR in Union_box named %s, zdepth is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (xwidth2 <= 0 && xwidth2 != -1) {
  printf("\nERROR in Union_box named %s, the xwidth2 is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (yheight2 <= 0 && yheight2 != -1) {
  printf("\nERROR in Union_box named %s, yheight2 is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}

// Use sanitation
if (global_material_list.num_elements == 0) {
  printf("\nERROR: Need to define a material using Union_make_material before using a Union geometry component. \n");
  printf("       %s was defined before first use of Union_make_material.\n",NAME_CURRENT_COMP);
  exit(1);
}

this_box_volume.geometry.is_masked_volume = 0;
this_box_volume.geometry.is_exit_volume = 0;
this_box_volume.geometry.is_mask_volume = 0;
// check if the volume is a mask, if it is the material string is irelevant.
if (mask_string && strlen(mask_string) && strcmp(mask_string, "NULL") && strcmp(mask_string, "0")) {
    // A mask volume is used to limit the extend of other volumes, called the masked volumes. These are specified in the mask_string.
    // In order for a ray to enter a masked volume, it needs to be both in the region covered by that volume AND the mask volume.
    // When more than
    this_box_volume.geometry.mask_mode = 1; // Default is mask mode is ALL
    if (mask_setting && strlen(mask_setting) && strcmp(mask_setting, "NULL") && strcmp(mask_setting, "0")) {
        if (strcmp(mask_setting,"ALL") == 0 || strcmp(mask_setting,"All") == 0) this_box_volume.geometry.mask_mode = 1;
        else if (strcmp(mask_setting,"ANY") == 0 || strcmp(mask_setting,"Any") == 0) this_box_volume.geometry.mask_mode = 2;
        else {
            printf("The mask_mode of component %s is set to %s, but must be either ALL or ANY.\n",NAME_CURRENT_COMP,mask_setting);
            exit(1);
        }
    }
    
    for (loop_index=0;loop_index<global_geometry_list.num_elements;loop_index++) {
        // Add mask list
        if (1 == manual_linking_function(global_geometry_list.elements[loop_index].name,mask_string)) {
            add_element_to_int_list(&this_box_volume.geometry.mask_list,global_geometry_list.elements[loop_index].component_index);
            add_element_to_int_list(&global_geometry_list.elements[loop_index].Volume->geometry.masked_by_list,INDEX_CURRENT_COMP);
            global_geometry_list.elements[loop_index].Volume->geometry.is_masked_volume = 1;
            if (this_box_volume.geometry.mask_mode == 2)
                global_geometry_list.elements[loop_index].Volume->geometry.mask_mode = 2;
            if (this_box_volume.geometry.mask_mode == 1) {
                if (global_geometry_list.elements[loop_index].Volume->geometry.is_masked_volume == 1 && global_geometry_list.elements[loop_index].Volume->geometry.mask_mode != 2)
                    // If more than one mask is added to one volume, the ANY mode overwrites the (default) ALL mode.
                    global_geometry_list.elements[loop_index].Volume->geometry.mask_mode = 1;
            }
            
            found_geometries = 1;
        }
    }
    if (found_geometries == 0) {
        printf("The mask_string in geometry: %s did not find any of the specified volumes in the mask_string %s \n",NAME_CURRENT_COMP,mask_string);
        exit(1);
    }
    this_box_volume.p_physics = malloc(sizeof(struct physics_struct));
    this_box_volume.p_physics->is_vacuum = 0; // Makes this volume a vacuum
    this_box_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
    this_box_volume.p_physics->my_a = 0; // Should not be used.
    sprintf(this_box_volume.p_physics->name,"Mask");
    this_box_volume.geometry.is_mask_volume = 1;
    
    
// Read the material input, or if it lacks, use automatic linking.
} else if (material_string && strlen(material_string) && strcmp(material_string, "NULL") && strcmp(material_string, "0")) {
    // A geometry string was given, use it to determine which material
    if (0 == strcmp(material_string,"vacuum") || 0 == strcmp(material_string,"Vacuum")) {
        // One could have a global physics struct for vacuum instead of creating one for each
        this_box_volume.p_physics = malloc(sizeof(struct physics_struct));
        this_box_volume.p_physics->is_vacuum = 1; // Makes this volume a vacuum
        this_box_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
        this_box_volume.p_physics->my_a = 0; // Should not be used.
        sprintf(this_box_volume.p_physics->name,"Vacuum");
    } else if (0 == strcmp(material_string,"exit") || 0 == strcmp(material_string,"Exit")) {
        // One could have a global physics struct for vacuum instead of creating one for each
        this_box_volume.p_physics = malloc(sizeof(struct physics_struct));
        this_box_volume.p_physics->is_vacuum = 1; // Makes this volume a vacuum
        this_box_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
        this_box_volume.p_physics->my_a = 0; // Should not be used.
        this_box_volume.geometry.is_exit_volume = 1;
        sprintf(this_box_volume.p_physics->name,"Exit");
    } else {
        #ifndef MATERIAL_DETECTOR
            printf("Need to define a material before refering to it in a geometry %s.\n",NAME_CURRENT_COMP);
            exit(1);
        #endif
        for (loop_index=0;loop_index<global_material_list.num_elements;loop_index++) {
            if (0 == strcmp(material_string,global_material_list.elements[loop_index].name)) {
                this_box_volume.p_physics = global_material_list.elements[loop_index].physics;
                break;
            }
            if (loop_index == global_material_list.num_elements-1) {
            printf("\n");
            printf("ERROR: The material string \"%s\" in Union geometry \"%s\" did not match a specified material. \n",material_string,NAME_CURRENT_COMP);
            printf("       The materials available at this point (need to be defined before the geometry): \n");
            for (loop_index=0;loop_index<global_material_list.num_elements;loop_index++)
              printf("         %s\n",global_material_list.elements[loop_index].name);
            printf("\n");
            printf("       It is also possible to use one of the defualt materials avaiable: \n");
            printf("           Vacuum (for a Volume without scattering or absorption)\n");
            printf("           Exit (for a Volume where the ray exits the component if it enters)\n");
            printf("           Mask (for a Volume that masks existing volumes specified in the mask_string\n");
            exit(1);
          }
        }
    }
} else {
    // Automatic linking, simply using the last defined material.
    #ifndef MATERIAL_DETECTOR
        printf("Need to define a material before the geometry to use automatic linking %s.\n",NAME_CURRENT_COMP);
        exit(1);
    #endif
    this_box_volume.p_physics = global_material_list.elements[global_material_list.num_elements-1].physics;
}

sprintf(this_box_volume.name,NAME_CURRENT_COMP);
sprintf(this_box_volume.geometry.shape,"box");
this_box_volume.geometry.priority_value = priority;
this_box_volume.geometry.geometry_p_interact = p_interact;
// Currently the coordinates will be in absolute space.
this_box_volume.geometry.center = POS_A_CURRENT_COMP;

this_box_storage.z_depth = zdepth;
this_box_storage.x_width1 = xwidth;
this_box_storage.y_height1 = yheight;

this_box_storage.is_rectangle = 0;
if (xwidth2 < 0 && yheight2 < 0) this_box_storage.is_rectangle = 1;
if (xwidth == xwidth2 && yheight == yheight2) this_box_storage.is_rectangle = 1;

if (xwidth2 < 0) {
    this_box_storage.x_width2 = xwidth;
    xwidth2 = xwidth;
} else this_box_storage.x_width2 = xwidth2;

if (yheight2 < 0) {
    this_box_storage.y_height2 = yheight;
    yheight2 = yheight;
} else this_box_storage.y_height2 = yheight2;


this_box_storage.normal_vectors[0] = coords_set(0,0,1);
this_box_storage.normal_vectors[1] = coords_set(0,0,1);

// for sides with y component = 0
x_component = 2*zdepth/sqrt((xwidth-xwidth2)*(xwidth-xwidth2)+4*zdepth*zdepth);
z_component = (xwidth-xwidth2)/sqrt(4*zdepth*zdepth+(xwidth-xwidth2)*(xwidth-xwidth2));

this_box_storage.normal_vectors[2] = coords_set(x_component,0,z_component);
this_box_storage.normal_vectors[3] = coords_set(-x_component,0,z_component);

// for sides with x component = 0
y_component = 2*zdepth/sqrt((yheight-yheight2)*(yheight-yheight2)+4*zdepth*zdepth);
z_component = (yheight-yheight2)/sqrt(4*zdepth*zdepth+(yheight-yheight2)*(yheight-yheight2));

this_box_storage.normal_vectors[4] = coords_set(0,y_component,z_component);
this_box_storage.normal_vectors[5] = coords_set(0,-y_component,z_component);

this_box_volume.geometry.visualization_on = visualize;

this_box_volume.geometry.geometry_parameters.p_box_storage = &this_box_storage;

// Assign pointers to functions for intersection with the shape, checking if a point is inside the shape
if (this_box_storage.is_rectangle == 1) {
this_box_volume.geometry.intersect_function = &sample_box_intersect_simple;
this_box_volume.geometry.within_function = &r_within_box_simple;
}
else {
this_box_volume.geometry.intersect_function = &sample_box_intersect_advanced;
this_box_volume.geometry.within_function = &r_within_box_advanced;
}

this_box_volume.geometry.shell_points = &box_shell_points;
this_box_volume.geometry.mcdisplay_function = &mcdisplay_box_function;
this_box_volume.geometry.initialize_from_main_function = &initialize_box_geometry_from_main_component;
this_box_volume.geometry.process_rot_allocated = 0;

this_box_volume.geometry.copy_geometry_parameters = &allocate_box_storage_copy;

rot_copy(this_box_volume.geometry.rotation_matrix,ROT_A_CURRENT_COMP); //check how ROT_R_CURRENT_COMP would work
rot_transpose(ROT_A_CURRENT_COMP,this_box_volume.geometry.transpose_rotation_matrix);

// Initialize loggers
this_box_volume.loggers.num_elements = 0;

// packing the information into the global_geometry_element, which is then included in the global_geometry_list.
sprintf(global_geometry_element.name,NAME_CURRENT_COMP);
global_geometry_element.activation_counter = number_of_activations;
global_geometry_element.component_index = INDEX_CURRENT_COMP;
global_geometry_element.Volume = &this_box_volume; // Would be nicer if this m was a pointer, now we have the (small) data two places
add_element_to_geometry_list(&global_geometry_list,global_geometry_element);
}
#line 14611 "Union_manual_example.c"
#undef number_of_activations
#undef mask_setting
#undef mask_string
#undef p_interact
#undef focus_r
#undef focus_xh
#undef focus_xw
#undef focus_ah
#undef focus_aw
#undef target_z
#undef target_y
#undef target_x
#undef target_index
#undef visualize
#undef yheight2
#undef xwidth2
#undef zdepth
#undef yheight
#undef xwidth
#undef priority
#undef material_string
#undef this_box_storage
#undef global_geometry_element
#undef this_box_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Sample. */
  SIG_MESSAGE("Sample (Init)");
#define mccompcurname  Sample
#define mccompcurtype  Union_cylinder
#define mccompcurindex 13
#define loop_index mccSample_loop_index
#define this_cylinder_volume mccSample_this_cylinder_volume
#define global_geometry_element mccSample_global_geometry_element
#define this_cylinder_storage mccSample_this_cylinder_storage
#define material_string mccSample_material_string
#define priority mccSample_priority
#define radius mccSample_radius
#define yheight mccSample_yheight
#define visualize mccSample_visualize
#define target_index mccSample_target_index
#define target_x mccSample_target_x
#define target_y mccSample_target_y
#define target_z mccSample_target_z
#define focus_aw mccSample_focus_aw
#define focus_ah mccSample_focus_ah
#define focus_xw mccSample_focus_xw
#define focus_xh mccSample_focus_xh
#define focus_r mccSample_focus_r
#define p_interact mccSample_p_interact
#define mask_string mccSample_mask_string
#define mask_setting mccSample_mask_setting
#define number_of_activations mccSample_number_of_activations
#line 184 "/usr/share/mcstas/2.5/contrib/union/Union_cylinder.comp"
{
// Initializes the focusing system for this volume including input sanitation.
focus_initialize(&this_cylinder_volume.geometry, POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index), POS_A_CURRENT_COMP, ROT_A_CURRENT_COMP, target_index, target_x, target_y, target_z, focus_aw, focus_ah, focus_xw, focus_xh, focus_r, NAME_CURRENT_COMP);

// Input sanitation for this geometry
if (radius <= 0) {
  printf("\nERROR in Union_cylinder named %s, the radius is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}

if (yheight <= 0) {
  printf("\nERROR in Union_cylinder named %s, yheight is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}

// Use sanitation
#ifdef MATERIAL_DETECTOR
if (global_material_list.num_elements == 0) {
  // Here if the user have defined a material, but only after this material
  printf("\nERROR: Need to define a material using Union_make_material before using a Union geometry component. \n");
  printf("       %s was defined before first use of Union_make_material.\n",NAME_CURRENT_COMP);
  exit(1);
}
#endif
#ifndef MATERIAL_DETECTOR
  printf("\nERROR: Need to define a material using Union_make_material before using a Union geometry component. \n");
  exit(1);
#endif


this_cylinder_volume.geometry.is_masked_volume = 0;
this_cylinder_volume.geometry.is_exit_volume = 0;
this_cylinder_volume.geometry.is_mask_volume = 0;

// Read the material input, or if it lacks, use automatic linking.
if (mask_string && strlen(mask_string) && strcmp(mask_string, "NULL") && strcmp(mask_string, "0")) {
    // A mask volume is used to limit the extend of other volumes, called the masked volumes. These are specified in the mask_string.
    // In order for a ray to enter a masked volume, it needs to be both in the region covered by that volume AND the mask volume.
    // When more than
    this_cylinder_volume.geometry.mask_mode = 1; // Default is mask mode is ALL
    if (mask_setting && strlen(mask_setting) && strcmp(mask_setting, "NULL") && strcmp(mask_setting, "0")) {
        if (strcmp(mask_setting,"ALL") == 0 || strcmp(mask_setting,"All") == 0) this_cylinder_volume.geometry.mask_mode = 1;
        else if (strcmp(mask_setting,"ANY") == 0 || strcmp(mask_setting,"Any") == 0) this_cylinder_volume.geometry.mask_mode = 2;
        else {
            printf("The mask_mode of component %s is set to %s, but must be either ALL or ANY.\n",NAME_CURRENT_COMP,mask_setting);
            exit(1);
        }
    }
    
    for (loop_index=0;loop_index<global_geometry_list.num_elements;loop_index++) {
        // Add mask list
        if (1 == manual_linking_function(global_geometry_list.elements[loop_index].name,mask_string)) {
            add_element_to_int_list(&this_cylinder_volume.geometry.mask_list,global_geometry_list.elements[loop_index].component_index);
            add_element_to_int_list(&global_geometry_list.elements[loop_index].Volume->geometry.masked_by_list,INDEX_CURRENT_COMP);
            global_geometry_list.elements[loop_index].Volume->geometry.is_masked_volume = 1;
            if (this_cylinder_volume.geometry.mask_mode == 2)
                global_geometry_list.elements[loop_index].Volume->geometry.mask_mode = 2;
            if (this_cylinder_volume.geometry.mask_mode == 1) {
                if (global_geometry_list.elements[loop_index].Volume->geometry.is_masked_volume == 1 && global_geometry_list.elements[loop_index].Volume->geometry.mask_mode != 2)
                    // If more than one mask is added to one volume, the ANY mode overwrites the (default) ALL mode.
                    global_geometry_list.elements[loop_index].Volume->geometry.mask_mode = 1;
            }
            
            found_geometries = 1;
        }
    }
    if (found_geometries == 0) {
        printf("The mask_string in geometry: %s did not find any of the specified volumes in the mask_string %s \n",NAME_CURRENT_COMP,mask_string);
        exit(1);
    }
    this_cylinder_volume.p_physics = malloc(sizeof(struct physics_struct));
    this_cylinder_volume.p_physics->is_vacuum = 0; // Makes this volume a vacuum
    this_cylinder_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
    this_cylinder_volume.p_physics->my_a = 0; // Should not be used.
    sprintf(this_cylinder_volume.p_physics->name,"Mask");
    this_cylinder_volume.geometry.is_mask_volume = 1;
    
    
// Read the material input, or if it lacks, use automatic linking.
} else if (material_string && strlen(material_string) && strcmp(material_string, "NULL") && strcmp(material_string, "0")) {
    // A geometry string was given, use it to determine which material
    if (0 == strcmp(material_string,"vacuum") || 0 == strcmp(material_string,"Vacuum")) {
        // One could have a global physics struct for vacuum instead of creating one for each
        this_cylinder_volume.p_physics = malloc(sizeof(struct physics_struct));
        this_cylinder_volume.p_physics->is_vacuum = 1; // Makes this volume a vacuum
        this_cylinder_volume.p_physics->number_of_processes = (int) 0;
        this_cylinder_volume.p_physics->my_a = 0; // Should not be used.
        sprintf(this_cylinder_volume.p_physics->name,"Vacuum");
    } else if (0 == strcmp(material_string,"exit") || 0 == strcmp(material_string,"Exit")) {
        // One could have a global physics struct for exit instead of creating one for each
        this_cylinder_volume.p_physics = malloc(sizeof(struct physics_struct));
        this_cylinder_volume.p_physics->is_vacuum = 1; // Makes this volume a vacuum
        this_cylinder_volume.p_physics->number_of_processes = (int) 0;
        this_cylinder_volume.p_physics->my_a = 0; // Should not be used.
        this_cylinder_volume.geometry.is_exit_volume = 1;
        sprintf(this_cylinder_volume.p_physics->name,"Exit");
    } else {
        for (loop_index=0;loop_index<global_material_list.num_elements;loop_index++) {
          if (0 == strcmp(material_string,global_material_list.elements[loop_index].name)) {
            this_cylinder_volume.p_physics = global_material_list.elements[loop_index].physics;
            break;
          }
          if (loop_index == global_material_list.num_elements-1) {
            printf("\n");
            printf("ERROR: The material string \"%s\" in Union geometry \"%s\" did not match a specified material. \n",material_string,NAME_CURRENT_COMP);
            printf("       The materials available at this point (need to be defined before the geometry): \n");
            for (loop_index=0;loop_index<global_material_list.num_elements;loop_index++)
              printf("         %s\n",global_material_list.elements[loop_index].name);
            printf("\n");
            printf("       It is also possible to use one of the defualt materials avaiable: \n");
            printf("           Vacuum (for a Volume without scattering or absorption)\n");
            printf("           Exit (for a Volume where the ray exits the component if it enters)\n");
            printf("           Mask (for a Volume that masks existing volumes specified in the mask_string\n");
            exit(1);
          }
        }
    }
} else {
    // Automatic linking, simply using the last defined material.
    #ifndef MATERIAL_DETECTOR
        printf("Need to define a material before the geometry to use automatic linking %s.\n",NAME_CURRENT_COMP);
        exit(1);
    #endif
    this_cylinder_volume.p_physics = global_material_list.elements[global_material_list.num_elements-1].physics;
}

sprintf(this_cylinder_volume.name,NAME_CURRENT_COMP);
sprintf(this_cylinder_volume.geometry.shape,"cylinder");
this_cylinder_volume.geometry.priority_value = priority;
// Currently the coordinates will be in absolute space.
this_cylinder_volume.geometry.center = POS_A_CURRENT_COMP;

this_cylinder_volume.geometry.geometry_p_interact = p_interact;
this_cylinder_storage.cyl_radius = radius;
this_cylinder_storage.height = yheight;
this_cylinder_volume.geometry.visualization_on = visualize;
this_cylinder_volume.geometry.geometry_parameters.p_cylinder_storage = &this_cylinder_storage;
this_cylinder_volume.geometry.within_function = &r_within_cylinder;
this_cylinder_volume.geometry.intersect_function = &sample_cylinder_intersect;
this_cylinder_volume.geometry.mcdisplay_function = &mcdisplay_cylinder_function;
this_cylinder_volume.geometry.shell_points = &cylinder_shell_points;
this_cylinder_volume.geometry.initialize_from_main_function = &initialize_cylinder_geometry_from_main_component;
this_cylinder_volume.geometry.process_rot_allocated = 0;
this_cylinder_volume.geometry.copy_geometry_parameters = &allocate_cylinder_storage_copy;

rot_copy(this_cylinder_volume.geometry.rotation_matrix,ROT_A_CURRENT_COMP);
rot_transpose(ROT_A_CURRENT_COMP,this_cylinder_volume.geometry.transpose_rotation_matrix);

// Initialize loggers
this_cylinder_volume.loggers.num_elements = 0;

// packing the information into the global_geometry_element, which is then included in the global_geometry_list.
sprintf(global_geometry_element.name,NAME_CURRENT_COMP);
global_geometry_element.activation_counter = number_of_activations;
global_geometry_element.component_index = INDEX_CURRENT_COMP;
global_geometry_element.Volume = &this_cylinder_volume; // Would be nicer if this m was a pointer, now we have the (small) data two places
add_element_to_geometry_list(&global_geometry_list,global_geometry_element);
}
#line 14827 "Union_manual_example.c"
#undef number_of_activations
#undef mask_setting
#undef mask_string
#undef p_interact
#undef focus_r
#undef focus_xh
#undef focus_xw
#undef focus_ah
#undef focus_aw
#undef target_z
#undef target_y
#undef target_x
#undef target_index
#undef visualize
#undef yheight
#undef radius
#undef priority
#undef material_string
#undef this_cylinder_storage
#undef global_geometry_element
#undef this_cylinder_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Master. */
  SIG_MESSAGE("Master (Init)");
#define mccompcurname  Master
#define mccompcurtype  Union_master
#define mccompcurindex 14
#define verbal mccMaster_verbal
#define list_verbal mccMaster_list_verbal
#define trace_verbal mccMaster_trace_verbal
#define finally_verbal mccMaster_finally_verbal
#define starting_volume_warning mccMaster_starting_volume_warning
#define global_master_element mccMaster_global_master_element
#define previous_master_index mccMaster_previous_master_index
#define geometry_list_index mccMaster_geometry_list_index
#define intersection_time_table mccMaster_intersection_time_table
#define Volumes mccMaster_Volumes
#define Geometries mccMaster_Geometries
#define starting_lists mccMaster_starting_lists
#define r mccMaster_r
#define r_start mccMaster_r_start
#define v mccMaster_v
#define error_msg mccMaster_error_msg
#define component_error_msg mccMaster_component_error_msg
#define string_output mccMaster_string_output
#define number_of_volumes mccMaster_number_of_volumes
#define volume_index mccMaster_volume_index
#define process_index mccMaster_process_index
#define solutions mccMaster_solutions
#define max_number_of_processes mccMaster_max_number_of_processes
#define limit mccMaster_limit
#define solution mccMaster_solution
#define min_solution mccMaster_min_solution
#define min_volume mccMaster_min_volume
#define time_found mccMaster_time_found
#define intersection_time mccMaster_intersection_time
#define min_intersection_time mccMaster_min_intersection_time
#define process mccMaster_process
#define process_start mccMaster_process_start
#define my_trace mccMaster_my_trace
#define p_my_trace mccMaster_p_my_trace
#define my_trace_fraction_control mccMaster_my_trace_fraction_control
#define k mccMaster_k
#define k_new mccMaster_k_new
#define v_length mccMaster_v_length
#define my_sum mccMaster_my_sum
#define my_sum_plus_abs mccMaster_my_sum_plus_abs
#define culmative_probability mccMaster_culmative_probability
#define mc_prop mccMaster_mc_prop
#define time_to_scattering mccMaster_time_to_scattering
#define length_to_scattering mccMaster_length_to_scattering
#define length_to_boundery mccMaster_length_to_boundery
#define time_to_boundery mccMaster_time_to_boundery
#define selected_process mccMaster_selected_process
#define scattering_event mccMaster_scattering_event
#define time_propagated_without_scattering mccMaster_time_propagated_without_scattering
#define a_next_volume_found mccMaster_a_next_volume_found
#define next_volume mccMaster_next_volume
#define next_volume_priority mccMaster_next_volume_priority
#define done mccMaster_done
#define current_volume mccMaster_current_volume
#define number_of_solutions mccMaster_number_of_solutions
#define number_of_solutions_static mccMaster_number_of_solutions_static
#define check mccMaster_check
#define start mccMaster_start
#define intersection_with_children mccMaster_intersection_with_children
#define geometry_output mccMaster_geometry_output
#define tree_next_volume mccMaster_tree_next_volume
#define pre_allocated1 mccMaster_pre_allocated1
#define pre_allocated2 mccMaster_pre_allocated2
#define pre_allocated3 mccMaster_pre_allocated3
#define ray_position mccMaster_ray_position
#define ray_velocity mccMaster_ray_velocity
#define ray_velocity_final mccMaster_ray_velocity_final
#define volume_0_found mccMaster_volume_0_found
#define scattered_flag mccMaster_scattered_flag
#define scattered_flag_VP mccMaster_scattered_flag_VP
#define master_transposed_rotation_matrix mccMaster_master_transposed_rotation_matrix
#define temp_rotation_matrix mccMaster_temp_rotation_matrix
#define non_rotated_position mccMaster_non_rotated_position
#define rotated_position mccMaster_rotated_position
#define enable_tagging mccMaster_enable_tagging
#define stop_tagging_ray mccMaster_stop_tagging_ray
#define stop_creating_nodes mccMaster_stop_creating_nodes
#define enable_tagging_check mccMaster_enable_tagging_check
#define master_tagging_node_list mccMaster_master_tagging_node_list
#define current_tagging_node mccMaster_current_tagging_node
#define tagging_leaf_counter mccMaster_tagging_leaf_counter
#define number_of_scattering_events mccMaster_number_of_scattering_events
#define real_transmission_probability mccMaster_real_transmission_probability
#define mc_transmission_probability mccMaster_mc_transmission_probability
#define number_of_masks mccMaster_number_of_masks
#define number_of_masked_volumes mccMaster_number_of_masked_volumes
#define need_to_run_within_which_volume mccMaster_need_to_run_within_which_volume
#define mask_index_main mccMaster_mask_index_main
#define mask_iterate mccMaster_mask_iterate
#define mask_status_list mccMaster_mask_status_list
#define current_mask_intersect_list_status mccMaster_current_mask_intersect_list_status
#define mask_volume_index_list mccMaster_mask_volume_index_list
#define geometry_component_index_list mccMaster_geometry_component_index_list
#define Volume_copies_allocated mccMaster_Volume_copies_allocated
#define loggers_with_data_array mccMaster_loggers_with_data_array
#define p_old mccMaster_p_old
#define this_logger mccMaster_this_logger
#define conditional_status mccMaster_conditional_status
#define tagging_conditional_list mccMaster_tagging_conditional_list
#define free_tagging_conditioanl_list mccMaster_free_tagging_conditioanl_list
#define logger_conditional_extend_array mccMaster_logger_conditional_extend_array
#define tagging_conditional_extend mccMaster_tagging_conditional_extend
#define max_conditional_extend_index mccMaster_max_conditional_extend_index
#define safty_distance mccMaster_safty_distance
#define safty_distance2 mccMaster_safty_distance2
#define number_of_processes_array mccMaster_number_of_processes_array
#define allow_inside_start mccMaster_allow_inside_start
#define history_limit mccMaster_history_limit
#line 192 "/usr/share/mcstas/2.5/contrib/union/Union_master.comp"
{
  // Use sanitation
  #ifndef ANY_GEOMETRY_DETECTOR_DECLARE
    printf("\nERROR: Need to define at least one Volume using Union_cylinder or Union_box before using the Union_master component. \n");
    exit(1);
  #endif
  #ifdef ANY_GEOMETRY_DETECTOR_DECLARE
    if (global_geometry_list.num_elements == 0) {
      printf("\nERROR: Need to define at least one Volume using Union_cylinder or Union_box before using the Union_master component. \n");
      printf("       Union_master component named \"%s\" is before any Volumes in the instrument file. At least one Volume need to be defined before\n",NAME_CURRENT_COMP);
    
      exit(1);
    }
  #endif
  
  safty_distance = 1E-11;
  safty_distance2 = safty_distance*2;

  sprintf(global_master_element.name,NAME_CURRENT_COMP);
  global_master_element.component_index = INDEX_CURRENT_COMP;
  add_element_to_master_list(&global_master_list,global_master_element);
  
  
  if (global_master_list.num_elements == 1) previous_master_index = 0; // no previous index
  else previous_master_index = global_master_list.elements[global_master_list.num_elements-2].component_index; // -2 because of zero indexing and needing the previous index.
  //printf("Assigned previous_master_index = %d \n",previous_master_index);
  
  number_of_volumes = 1; // Starting with 1 as the surrounding vacuum is considered a volume
  number_of_masks = 0;   // Starting with 0 mask volumes
  number_of_masked_volumes = 0; // Starting with 0 masked volumes
  for (iterate=0;iterate<global_geometry_list.num_elements;iterate++) {
    if (global_geometry_list.elements[iterate].component_index < INDEX_CURRENT_COMP && global_geometry_list.elements[iterate].activation_counter > 0) {
        global_geometry_list.elements[iterate].active = 1;
        global_geometry_list.elements[iterate].activation_counter--;
        number_of_volumes++;
        if (global_geometry_list.elements[iterate].Volume->geometry.is_mask_volume == 1) number_of_masks++;
        if (global_geometry_list.elements[iterate].Volume->geometry.is_masked_volume == 1) number_of_masked_volumes++;
    } else global_geometry_list.elements[iterate].active = 0;
  }
  //printf("Found number of volumes to be %d \n",number_of_volumes);
  
  // Allocation of global lists
  geometry_component_index_list.num_elements = number_of_volumes;
  geometry_component_index_list.elements = malloc( geometry_component_index_list.num_elements * sizeof(int));
  mask_volume_index_list.num_elements = number_of_masks;
  if (number_of_masks >0) mask_volume_index_list.elements = malloc( number_of_masks * sizeof(int));
  mask_status_list.num_elements = number_of_masks;
  if (number_of_masks >0) mask_status_list.elements = malloc( number_of_masks * sizeof(int));
  current_mask_intersect_list_status.num_elements = number_of_masked_volumes;
  if (number_of_masked_volumes >0) current_mask_intersect_list_status.elements = malloc( number_of_masked_volumes * sizeof(int));
  
  volume_index = 0;
  for (iterate=0;iterate<global_geometry_list.num_elements;iterate++) {
    if (global_geometry_list.elements[iterate].active == 1)
        geometry_component_index_list.elements[++volume_index] = global_geometry_list.elements[iterate].component_index;
      
  }
  geometry_component_index_list.elements[0] = 0; // Volume 0 is never set in the above code, but should never be used.
  

  MPI_MASTER(
  // The input for this component is done through a series of input components
  // All information needed is stored in global lists, which are printed here.
  if (verbal == 1) {
      printf("---------------------------------------------------------------------\n");
      printf("global_process_list.num_elements: %d\n",global_process_list.num_elements);
      for (iterate=0;iterate<global_process_list.num_elements;iterate++) {
      printf("name of process [%d]: %s \n",iterate,global_process_list.elements[iterate].name);
      printf("component index [%d]: %d \n",iterate,global_process_list.elements[iterate].component_index);
      }
      
      printf("---------------------------------------------------------------------\n");
      printf("global_material_list.num_elements: %d\n",global_material_list.num_elements);
      for (iterate=0;iterate<global_material_list.num_elements;iterate++) {
      printf("name of material    [%d]: %s \n",iterate,global_material_list.elements[iterate].name);
      printf("component index     [%d]: %d \n",iterate,global_material_list.elements[iterate].component_index);
      printf("my_absoprtion       [%d]: %f \n",iterate,global_material_list.elements[iterate].physics->my_a);
      printf("number of processes [%d]: %d \n",iterate,global_material_list.elements[iterate].physics->number_of_processes);
      }

      printf("---------------------------------------------------------------------\n");
      printf("global_geometry_list.num_elements: %d\n",global_material_list.num_elements);
      for (iterate=0;iterate<global_geometry_list.num_elements;iterate++) {
        if (global_geometry_list.elements[iterate].active == 1) {
          printf("\n");
          printf("name of geometry    [%d]: %s \n",iterate,global_geometry_list.elements[iterate].name);
          printf("component index     [%d]: %d \n",iterate,global_geometry_list.elements[iterate].component_index);
          printf("Volume.name         [%d]: %s \n",iterate,global_geometry_list.elements[iterate].Volume->name);
          if (global_geometry_list.elements[iterate].Volume->geometry.is_mask_volume == 0) {
          printf("Volume.p_physics.is_vacuum           [%d]: %d \n",iterate,global_geometry_list.elements[iterate].Volume->p_physics->is_vacuum);
          printf("Volume.p_physics.my_absoprtion       [%d]: %f \n",iterate,global_geometry_list.elements[iterate].Volume->p_physics->my_a);
          printf("Volume.p_physics.number of processes [%d]: %d \n",iterate,global_geometry_list.elements[iterate].Volume->p_physics->number_of_processes);
          }
          printf("Volume.geometry.shape                [%d]: %s \n",iterate,global_geometry_list.elements[iterate].Volume->geometry.shape);
          printf("Volume.geometry.center.x             [%d]: %f \n",iterate,global_geometry_list.elements[iterate].Volume->geometry.center.x);
          printf("Volume.geometry.center.y             [%d]: %f \n",iterate,global_geometry_list.elements[iterate].Volume->geometry.center.y);
          printf("Volume.geometry.center.z             [%d]: %f \n",iterate,global_geometry_list.elements[iterate].Volume->geometry.center.z);
          printf("Volume.geometry.rotation_matrix[0]           [%d]: [%f %f %f] \n",iterate,global_geometry_list.elements[iterate].Volume->geometry.rotation_matrix[0][0],global_geometry_list.elements[iterate].Volume->geometry.rotation_matrix[0][1],global_geometry_list.elements[iterate].Volume->geometry.rotation_matrix[0][2]);
          printf("Volume.geometry.rotation_matrix[1]           [%d]: [%f %f %f] \n",iterate,global_geometry_list.elements[iterate].Volume->geometry.rotation_matrix[1][0],global_geometry_list.elements[iterate].Volume->geometry.rotation_matrix[1][1],global_geometry_list.elements[iterate].Volume->geometry.rotation_matrix[1][2]);
          printf("Volume.geometry.rotation_matrix[2]           [%d]: [%f %f %f] \n",iterate,global_geometry_list.elements[iterate].Volume->geometry.rotation_matrix[2][0],global_geometry_list.elements[iterate].Volume->geometry.rotation_matrix[2][1],global_geometry_list.elements[iterate].Volume->geometry.rotation_matrix[2][2]);
          if (strcmp(global_geometry_list.elements[iterate].Volume->geometry.shape,"cylinder") == 0) {
          printf("Volume.geometry.geometry_parameters.cyl_radius [%d]: %f \n",iterate,global_geometry_list.elements[iterate].Volume->geometry.geometry_parameters.p_cylinder_storage->cyl_radius);
          printf("Volume.geometry.geometry_parameters.height [%d]: %f \n",iterate,global_geometry_list.elements[iterate].Volume->geometry.geometry_parameters.p_cylinder_storage->height);
          }
          printf("Volume.geometry.focus_data.Aim             [%d]: [%f %f %f] \n",iterate,global_geometry_list.elements[iterate].Volume->geometry.focus_data.Aim.x,global_geometry_list.elements[iterate].Volume->geometry.focus_data.Aim.y,global_geometry_list.elements[iterate].Volume->geometry.focus_data.Aim.z);
        }
      }
      printf("---------------------------------------------------------------------\n");
      printf("number_of_volumes = %d\n",number_of_volumes);
      printf("number_of_masks = %d\n",number_of_masks);
      printf("number_of_masked_volumes = %d\n",number_of_masked_volumes);
  }
  ); // End MPI_MASTER

  
  // This information is to be stored in one array of structures that is allocated here
  Volumes = malloc(number_of_volumes * sizeof(struct Volume_struct*)); // BUG, should be Volume_struct*, 23/9/2016
  scattered_flag = malloc(number_of_volumes*sizeof(int)); // BUG, was sizeof(struct Volume_struct) 23/9/2016
  scattered_flag_VP = (int**) malloc(number_of_volumes * sizeof(int*));
  number_of_processes_array = malloc(number_of_volumes*sizeof(int));
  
  // The plotting functions need access to the other geomtries, but can not use the Volumes struct because of order of definition.
  Geometries = malloc(number_of_volumes * sizeof(struct geometry_struct *));
  
  // When activation counter is used to have several copies of one volume, it can become necessary to have soft copies of volumes
  // Not necessarily all of these will be allocated or used.
  Volume_copies = malloc(number_of_volumes * sizeof(struct Volume_struct *));
  Volume_copies_allocated.num_elements = 0;
  
  // The central structure is called a "Volume", tt describes a region in space with certain scattering processes and absorption cross section

  // ---  Volume 0 ---------------------------------------------------------------------------------
  // Volume zero is the vacuum surrounding the experiment (infinite, everywhere)
  Volumes[0] = malloc(sizeof(struct Volume_struct));
  strcpy(Volumes[0]->name,"Surrounding vacuum");
  // Assign geometry
  
  // This information is meaningless for vacuum, and should never be acsessed in the logic.
  Volumes[0]->geometry.priority_value = 0.0;
  Volumes[0]->geometry.center.x = 0;
  Volumes[0]->geometry.center.y = 0;
  Volumes[0]->geometry.center.z = 0;
  strcpy(Volumes[0]->geometry.shape,"vacuum");
  Volumes[0]->geometry.within_function = &r_within_surroundings;
  // No physics struct allocated
  Volumes[0]->p_physics = NULL;
  number_of_processes_array[volume_index] = 0;
  
  // These are never used for volume 0, but by setting the length to 0 it is automatically skipped in many forloops without the need for a if statement
  Volumes[0]->geometry.children.num_elements=0;
  Volumes[0]->geometry.direct_children.num_elements=0;
  Volumes[0]->geometry.destinations_list.num_elements=0;
  Volumes[0]->geometry.reduced_destinations_list.num_elements=0;
  
  Volumes[0]->geometry.masked_by_list.num_elements = 0;
  Volumes[0]->geometry.mask_list.num_elements = 0;
  Volumes[0]->geometry.masked_by_mask_index_list.num_elements = 0;
  Volumes[0]->geometry.mask_mode=0;
  Volumes[0]->geometry.is_mask_volume=0;
  Volumes[0]->geometry.is_masked_volume=0;
  
  
  // Store a pointer to the conditional list and update the current index;
  if (global_tagging_conditional_list.num_elements == global_tagging_conditional_list.current_index + 1) {
    tagging_conditional_list = &global_tagging_conditional_list.elements[global_tagging_conditional_list.current_index++].conditional_list;
    free_tagging_conditioanl_list = 0;
  } else {
    tagging_conditional_list = malloc(sizeof(struct conditional_list_struct));
    tagging_conditional_list->num_elements = 0;
    free_tagging_conditioanl_list = 1;
  }
  
    // The absolute rotation of this component is saved for use in initialization
  rot_transpose(ROT_A_CURRENT_COMP,master_transposed_rotation_matrix);
  
  // Update positions and rotations in transform lists
  for (iterate=0;iterate<global_positions_to_transform_list.num_elements;iterate++) {
    //print_position(*global_positions_to_transform_list.positions[iterate],"Position to transform before");
    non_rotated_position = coords_sub(*(global_positions_to_transform_list.positions[iterate]),POS_A_CURRENT_COMP);
    *(global_positions_to_transform_list.positions[iterate]) = rot_apply(ROT_A_CURRENT_COMP,non_rotated_position);
    //print_position(*global_positions_to_transform_list.positions[iterate],"Position to transform after");
  }
  // Free the list so that the Union components after this master get their own list
  if (global_positions_to_transform_list.num_elements > 0) {
    global_positions_to_transform_list.num_elements = 0;
    free(global_positions_to_transform_list.positions);
  }
  for (iterate=0;iterate<global_rotations_to_transform_list.num_elements;iterate++) {
    //print_rotation(*(global_rotations_to_transform_list.rotations[iterate]),"rotation matrix to be updated");
    
    //rot_mul(ROT_A_CURRENT_COMP,*(global_rotations_to_transform_list.rotations[iterate]),temp_rotation_matrix);
    rot_mul(master_transposed_rotation_matrix,*(global_rotations_to_transform_list.rotations[iterate]),temp_rotation_matrix);
    rot_copy(*(global_rotations_to_transform_list.rotations[iterate]),temp_rotation_matrix);
  }
  if (global_rotations_to_transform_list.num_elements > 0) {
    global_rotations_to_transform_list.num_elements = 0;
    free(global_rotations_to_transform_list.rotations);
  }
  
  // A pointer to the geometry structure
  Geometries[0] = &Volumes[0]->geometry;
  
  // Logging initialization
  Volumes[0]->loggers.num_elements = 0;
  
  // Initialize all Volumes that originate from input components
  max_number_of_processes = 0;
  max_conditional_extend_index = -1;
  for (iterate=0;iterate<global_all_volume_logger_list.num_elements;iterate++) {
    if (global_all_volume_logger_list.elements[iterate].logger->logger_extend_index > max_conditional_extend_index) {
      max_conditional_extend_index = global_all_volume_logger_list.elements[iterate].logger->logger_extend_index;
    }
  }
    
    
  volume_index = 0;
  mask_index_main = 0;
  for (geometry_list_index=0;geometry_list_index<global_geometry_list.num_elements;geometry_list_index++) {
    if (global_geometry_list.elements[geometry_list_index].active == 1) {
      volume_index++;
      // Connect a volume for each of the geometry.comp instances in the McStas instrument files
      if (global_geometry_list.elements[geometry_list_index].activation_counter == 0) {
        // This is the last time this volume is used, use the hard copy from the geometry component
        Volumes[volume_index] = global_geometry_list.elements[geometry_list_index].Volume;
        //printf("used hard copy of volume %d \n",volume_index);
      } else {
        //printf("\n");
        //printf("making local copy of volume %d \n",volume_index);
        // Since this volume is still needed more than this once, we need to make a shallow copy and use instead
        Volume_copies[volume_index] = malloc(sizeof(struct Volume_struct));
        *(Volume_copies[volume_index]) = *global_geometry_list.elements[geometry_list_index].Volume; // Makes shallow copy
        Volumes[volume_index] = Volume_copies[volume_index];
        add_element_to_int_list(&Volume_copies_allocated,volume_index);
        // Here I need to test how deep this copy is.
        //  - Are matrix elements transfered, or are they stored as pointers? YES
        //  - Are Coords transfered, or are they stored as pointers? YES
        // Needed hard copy: center (Coords center, Rotation rotation_matrix (+transpose)) ALL WORKS
        
        // The geometry storage needs a shallow copy as well (hard copy not necessary for any current geometries)
        // A simple copy_geometry_parameters function is added to the geometry in each geometry component
        Volumes[volume_index]->geometry.geometry_parameters = Volumes[volume_index]->geometry.copy_geometry_parameters(&global_geometry_list.elements[geometry_list_index].Volume->geometry.geometry_parameters);
        
        /*
        // These tests show that the shallow copy copes matrices and arrays, but not two layer deep structes, meaning p_physics is the same for the copy and original, which is perfect for our purposes.
        // Lists are also copied
        // The main problem is the geometry storage structure, which should be copied, but isn't
        
        print_position(Volumes[volume_index]->geometry.center,"Volume copy center");
        print_position(global_geometry_list.elements[geometry_list_index].Volume->geometry.center,"Volume original center");
        
        printf("changed center of original to 0 0 0 \n");
        global_geometry_list.elements[geometry_list_index].Volume->geometry.center = coords_set(0,0,0);
        
        print_position(Volumes[volume_index]->geometry.center,"Volume copy center");
        print_position(global_geometry_list.elements[geometry_list_index].Volume->geometry.center,"Volume original center");
        
        printf("\n");
        print_rotation(Volumes[volume_index]->geometry.rotation_matrix,"Volume copy rotation_matrix");
        print_rotation(global_geometry_list.elements[geometry_list_index].Volume->geometry.rotation_matrix,"Volume original rotation_matrix");
        
        printf("changed rotation matrix [0][1] original 15 \n");
        global_geometry_list.elements[geometry_list_index].Volume->geometry.rotation_matrix[0][1] = 15.0;
        
        print_rotation(Volumes[volume_index]->geometry.rotation_matrix,"Volume copy rotation_matrix");
        print_rotation(global_geometry_list.elements[geometry_list_index].Volume->geometry.rotation_matrix,"Volume original rotation_matrix");
        
        print_1d_int_list(Volumes[volume_index]->geometry.mask_list,"Volume copy mask_list");
        print_1d_int_list(global_geometry_list.elements[geometry_list_index].Volume->geometry.mask_list,"Volume original mask_list");
        
        printf("Addding element 15 to mask list of volume \n");
        add_element_to_int_list(&global_geometry_list.elements[geometry_list_index].Volume->geometry.mask_list,15);
    
        print_1d_int_list(Volumes[volume_index]->geometry.mask_list,"Volume copy mask_list");
        print_1d_int_list(global_geometry_list.elements[geometry_list_index].Volume->geometry.mask_list,"Volume original mask_list");
        
        printf("my_a copy physical material: %f \n",Volumes[volume_index]->p_physics->my_a);
        printf("my_a of original physical material: 0: %f \n",global_geometry_list.elements[geometry_list_index].Volume->p_physics->my_a);
        
        printf("changeing originals my_a for physical material \n");
        global_geometry_list.elements[geometry_list_index].Volume->p_physics->my_a = 15.0;
        
        printf("my_a copy physical material: %f \n",Volumes[volume_index]->p_physics->my_a);
        printf("my_a of original physical material: 0: %f \n",global_geometry_list.elements[geometry_list_index].Volume->p_physics->my_a);
        
        
        
        
        
        //printf("name of copy physical process 0: %s \n",Volumes[volume_index]->p_physics->name);
        //printf("name of original physical process 0: %s \n",(global_geometry_list.elements[geometry_list_index].Volume->p_physics->name);
        
        printf("is_exit_volume = %d : is_mask_volume = %d \n",Volumes[volume_index]->geometry.is_exit_volume,Volumes[volume_index]->geometry.is_mask_volume);
        
        printf("\n");
        */
        
        
        /*
        // Test showing the issue, the geometry_parameters point to the same file, and need to be copied instead
        
        // Fixed by adding this line
        Volumes[volume_index]->geometry.geometry_parameters = Volumes[volume_index]->geometry.copy_geometry_parameters(&global_geometry_list.elements[geometry_list_index].Volume->geometry.geometry_parameters);
        
        print_position(Volumes[volume_index]->geometry.geometry_parameters.p_cylinder_storage->direction_vector,"Volume copy geometry_storage.direction_vector");
        print_position(global_geometry_list.elements[geometry_list_index].Volume->geometry.geometry_parameters.p_cylinder_storage->direction_vector,"Volume original geometry_storage.direction_vector");
        
        printf("changed geometry_storage.direction_vector of original to 15 16 17 \n");
        global_geometry_list.elements[geometry_list_index].Volume->geometry.geometry_parameters.p_cylinder_storage->direction_vector = coords_set(15,16,17);
        
        print_position(Volumes[volume_index]->geometry.geometry_parameters.p_cylinder_storage->direction_vector,"Volume copy geometry_storage.direction_vector");
        print_position(global_geometry_list.elements[geometry_list_index].Volume->geometry.geometry_parameters.p_cylinder_storage->direction_vector,"Volume original geometry_storage.direction_vector");
        */
        
          
      }
      
      // This section identifies the different non isotropic processes in the current volume and give them appropriate transformation matrices
      // Identify the number of non isotropic processes in a material (this code can be safely executed for the same material many times)
      non_isotropic_found = 0;
      for (iterate=0;iterate<Volumes[volume_index]->p_physics->number_of_processes;iterate++) {
        if (Volumes[volume_index]->p_physics->p_scattering_array[iterate].non_isotropic_rot_index != -1) {
            Volumes[volume_index]->p_physics->p_scattering_array[iterate].non_isotropic_rot_index = non_isotropic_found;
            non_isotropic_found++;
        }
      }
      // If there were any, allocate memory to hold all of them in the volume struct
      if (non_isotropic_found > 0) {
        if (Volumes[volume_index]->geometry.process_rot_allocated == 0) {
          Volumes[volume_index]->geometry.process_rot_matrix_array = malloc(non_isotropic_found * sizeof(Rotation));
          Volumes[volume_index]->geometry.transpose_process_rot_matrix_array = malloc(non_isotropic_found * sizeof(Rotation));
          Volumes[volume_index]->geometry.process_rot_allocated = 1;
        }
      
        non_isotropic_found = 0;
        for (iterate=0;iterate<Volumes[volume_index]->p_physics->number_of_processes;iterate++) {
          if (Volumes[volume_index]->p_physics->p_scattering_array[iterate].non_isotropic_rot_index != -1) {
            
            // Create rotation matrix for this specific volume / process combination to transform from master coordinate system to the non-isotropics process coordinate system
            // This is done by multipling the transpose master component roration matrix, the volume rotation, and then the process rotation matrix onto the velocity / wavevector
            rot_mul(Volumes[volume_index]->geometry.rotation_matrix,master_transposed_rotation_matrix,temp_rotation_matrix);
            rot_mul(Volumes[volume_index]->p_physics->p_scattering_array[iterate].rotation_matrix,temp_rotation_matrix,Volumes[volume_index]->geometry.process_rot_matrix_array[non_isotropic_found]);
            
            // Need to transpose as well to transform back to the master coordinate system
            rot_transpose(Volumes[volume_index]->geometry.process_rot_matrix_array[non_isotropic_found],Volumes[volume_index]->geometry.transpose_process_rot_matrix_array[non_isotropic_found]);

            // Debug print
            //print_rotation(Volumes[volume_index]->geometry.process_rot_matrix_array[non_isotropic_found],"Process rotation matrix");
            //print_rotation(Volumes[volume_index]->geometry.transpose_process_rot_matrix_array[non_isotropic_found],"Transpose process rotation matrix");
            
            non_isotropic_found++;
          }
        }
      }
      
      // This component works in its local coordinate system, and thus all information from the input components should be transformed to its coordinate system.
      // All the input components saved their absolute rotation/position into their Volume structure, and the absolute rotation of the current component is known.
      // The next section finds the relative rotation and translation of all the volumes and the master component.
        
      // transpose(R_V)*R_C
      rot_mul(ROT_A_CURRENT_COMP,Volumes[volume_index]->geometry.transpose_rotation_matrix,temp_rotation_matrix);
      // Copy the result back to the volumes structure
      rot_copy(Volumes[volume_index]->geometry.rotation_matrix,temp_rotation_matrix);
      // Now update the transpose as well
      rot_transpose(Volumes[volume_index]->geometry.rotation_matrix,temp_rotation_matrix);
      rot_copy(Volumes[volume_index]->geometry.transpose_rotation_matrix,temp_rotation_matrix);
        
      non_rotated_position.x = Volumes[volume_index]->geometry.center.x - POS_A_CURRENT_COMP.x;
      non_rotated_position.y = Volumes[volume_index]->geometry.center.y - POS_A_CURRENT_COMP.y;
      non_rotated_position.z = Volumes[volume_index]->geometry.center.z - POS_A_CURRENT_COMP.z;

      rot_transpose(ROT_A_CURRENT_COMP,temp_rotation_matrix);
      rotated_position = rot_apply(ROT_A_CURRENT_COMP,non_rotated_position);
      //rotated_position = rot_apply(Volumes[volume_index]->geometry.rotation_matrix,non_rotated_position);

      Volumes[volume_index]->geometry.center.x = rotated_position.x;
      Volumes[volume_index]->geometry.center.y = rotated_position.y;
      Volumes[volume_index]->geometry.center.z = rotated_position.z;
      
      // The focus_data information need to be updated as well
      rot_mul(ROT_A_CURRENT_COMP,Volumes[volume_index]->geometry.focus_data.absolute_rotation,temp_rotation_matrix);
      // Copy the result back to the volumes structure
      rot_copy(Volumes[volume_index]->geometry.focus_data.absolute_rotation,temp_rotation_matrix);
      // Now update the transpose as well
      //rot_transpose(Volumes[volume_index]->geometry.focus_data.absolute_rotation,temp_rotation_matrix);
      //rot_copy(Volumes[volume_index]->geometry.transpose_rotation_matrix,temp_rotation_matrix);
      
      // Using the same algorithm, the Aim vector of the focus_data struct is transformed as well
      // There may be an issue with focusing in a isotropic and non isotropic scattering process
      non_rotated_position.x = Volumes[volume_index]->geometry.focus_data.Aim.x;
      non_rotated_position.y = Volumes[volume_index]->geometry.focus_data.Aim.y;
      non_rotated_position.z = Volumes[volume_index]->geometry.focus_data.Aim.z;

      rotated_position = rot_apply(ROT_A_CURRENT_COMP,non_rotated_position);
      //rotated_position = rot_apply(master_transposed_rotation_matrix,non_rotated_position);

      /*
      Volumes[volume_index]->geometry.focus_data.Aim.x = rotated_position.x;
      Volumes[volume_index]->geometry.focus_data.Aim.y = rotated_position.y;
      Volumes[volume_index]->geometry.focus_data.Aim.z = rotated_position.z;
      */

      //sprintf(string_output,"Rotated Aim vector for volume [%d]",volume_index);
      //print_position(Volumes[volume_index]->geometry.focus_data.Aim,string_output);
      
        
      // To allocate enough memory to hold information on all processes, the maximum of these is found
      if (Volumes[volume_index]->p_physics->number_of_processes > max_number_of_processes)
          max_number_of_processes = Volumes[volume_index]->p_physics->number_of_processes;
        
      // Allocate memory to scattered_flag_VP
      scattered_flag_VP[volume_index] = malloc(Volumes[volume_index]->p_physics->number_of_processes * sizeof(int));
      number_of_processes_array[volume_index] = Volumes[volume_index]->p_physics->number_of_processes;
      
      // Normalizing and error checking process interact fraction
      number_of_process_interacts_set = 0; total_process_interact=0;
      for (process_index=0;process_index<Volumes[volume_index]->p_physics->number_of_processes;process_index++) {
        if (Volumes[volume_index]->p_physics->p_scattering_array[process_index].process_p_interact != -1) {
          number_of_process_interacts_set++;
          total_process_interact += Volumes[volume_index]->p_physics->p_scattering_array[process_index].process_p_interact;
        } else {
          index_of_lacking_process = process_index;
        }
      }
      
      if (number_of_process_interacts_set == 0) Volumes[volume_index]->p_physics->interact_control = 0;
      else Volumes[volume_index]->p_physics->interact_control = 1;
      
      // If all are set, check if they need renormalization so that the sum is one.
      if (number_of_process_interacts_set == Volumes[volume_index]->p_physics->number_of_processes) {
        if (total_process_interact > 1.001 || total_process_interact < 0.999) {
          for (process_index=0;process_index<Volumes[volume_index]->p_physics->number_of_processes;process_index++) {
            Volumes[volume_index]->p_physics->p_scattering_array[process_index].process_p_interact = Volumes[volume_index]->p_physics->p_scattering_array[process_index].process_p_interact/total_process_interact;
          }
        }
      } else if ( number_of_process_interacts_set != 0) {
        if (number_of_process_interacts_set == Volumes[volume_index]->p_physics->number_of_processes - 1) // If all but one is set, it is an easy fix
          Volumes[volume_index]->p_physics->p_scattering_array[index_of_lacking_process].process_p_interact = 1 - total_process_interact;
        else {
            printf("ERROR, material %s needs to have all, all minus one or none of its processes with an interact_fraction \n",Volumes[volume_index]->p_physics->name);
            exit(1);
        }
      }
      
      //if (verbal) {
      //  printf("Process interact fraction for volume_index = %d\n",volume_index);
      //  for (process_index=0;process_index<Volumes[volume_index]->p_physics->number_of_processes;process_index++) {
      //      printf("process interact for %d out of %d = %f \n",process_index ,Volumes[volume_index]->p_physics->number_of_processes ,Volumes[volume_index]->p_physics->p_scattering_array[process_index].process_p_interact);
      //  }
      //}
        
      // Some initialization can only happen after the rotation matrix relative to the master is known
      // Such initialization is placed in the geometry component, and executed here through a function pointer
      Volumes[volume_index]->geometry.initialize_from_main_function(&Volumes[volume_index]->geometry);
      
      /*
      Coords simple_vector;
      Coords cyl_vector;
      if (strcmp("cylinder",Volumes[volume_index]->geometry.shape) == 0) {
          
          // This section should be moved to a initialize_from_main function like for the box volume
          // Start with vector that points along the cylinder in the simple frame, and rotate to global
          simple_vector = coords_set(0,1,0);

          // Rotate the position of the neutron around the center of the cylinder
          cyl_vector = rot_apply(Volumes[volume_index]->geometry.rotation_matrix,simple_vector);
          NORM(cyl_vector.x,cyl_vector.y,cyl_vector.z);
          Volumes[volume_index]->geometry.geometry_parameters.p_cylinder_storage->direction_vector.x = cyl_vector.x;
          Volumes[volume_index]->geometry.geometry_parameters.p_cylinder_storage->direction_vector.y = cyl_vector.y;
          Volumes[volume_index]->geometry.geometry_parameters.p_cylinder_storage->direction_vector.z = cyl_vector.z;
          // if (verbal == 1) printf("Cords vector1 = (%f,%f,%f)\n",cyl_vector.x,cyl_vector.y,cyl_vector.z);
          
          Volumes[volume_index]->geometry.initialize_from_main_function(&Volumes[volume_index]->geometry);
      }
      if (strcmp("box",Volumes[volume_index]->geometry.shape) == 0) {
          // Sets the vectors in the coordinate system of the main component (this component)
          Volumes[volume_index]->geometry.initialize_from_main_function(&Volumes[volume_index]->geometry);
      }
      */
      
      // Add pointer to geometry to Geometries
      Geometries[volume_index] = &Volumes[volume_index]->geometry;
      
      // Initialize mask intersect list
      Volumes[volume_index]->geometry.mask_intersect_list.num_elements = 0;
      
      // Here the mask_list and masked_by_list for the volume is updated from being component index values to being volume indexes
      for (iterate=0;iterate<Volumes[volume_index]->geometry.mask_list.num_elements;iterate++)
        Volumes[volume_index]->geometry.mask_list.elements[iterate] = find_on_int_list(geometry_component_index_list,Volumes[volume_index]->geometry.mask_list.elements[iterate]);
      
      for (iterate=0;iterate<Volumes[volume_index]->geometry.masked_by_list.num_elements;iterate++)
        Volumes[volume_index]->geometry.masked_by_list.elements[iterate] = find_on_int_list(geometry_component_index_list,Volumes[volume_index]->geometry.masked_by_list.elements[iterate]);
        

      // These should maybe be in mask_indices instead, but not before after list generation
      
      // If the volume is a mask, it's volume number is added to the mask_volume_index list so volume index can be converted to mask_index.
      if (Volumes[volume_index]->geometry.is_mask_volume == 1) Volumes[volume_index]->geometry.mask_index = mask_index_main;
      if (Volumes[volume_index]->geometry.is_mask_volume == 1) mask_volume_index_list.elements[mask_index_main++] = volume_index;
     
      // Check all loggers assosiated with this volume
      //printf("reached max_test for volume %d \n",volume_index);
      for (iterate=0;iterate<Volumes[volume_index]->loggers.num_elements;iterate++) {
        //printf("iterate = %d \n",iterate);
        for (process_index=0;process_index<Volumes[volume_index]->loggers.p_logger_volume[iterate].num_elements;process_index++) {
        //printf("process_index = %d \n",process_index);
          if (Volumes[volume_index]->loggers.p_logger_volume[iterate].p_logger_process[process_index] != NULL) {
            if (Volumes[volume_index]->loggers.p_logger_volume[iterate].p_logger_process[process_index]->logger_extend_index > max_conditional_extend_index)
              max_conditional_extend_index = Volumes[volume_index]->loggers.p_logger_volume[iterate].p_logger_process[process_index]->logger_extend_index;
          }
        }
      }
      //printf("did max_test for volume %d\n",volume_index);
      
    
    }
  } // Initialization for each volume done
  
  
  
  my_trace = malloc(max_number_of_processes*sizeof(double));
  my_trace_fraction_control = malloc(max_number_of_processes*sizeof(double));
  number_of_solutions = &number_of_solutions_static;
  component_error_msg = 0;
  
  // Pre allocated memory for destination list search
  pre_allocated1 = malloc(number_of_volumes * sizeof(int));
  pre_allocated2 = malloc(number_of_volumes * sizeof(int));
  pre_allocated3 = malloc(number_of_volumes * sizeof(int));
  
  
  // Allocate memory for logger_conditional_extend_array used in the extend section of the master component, if it is needed.
  if (max_conditional_extend_index > -1) {
    logger_conditional_extend_array = malloc((max_conditional_extend_index + 1)*sizeof(int));
  }
  
  
  // In this function different lists of volume indecies are generated. They are the key to the speed of the component and central for the logic.
  // They use simple set algebra to generate these lists for each volume:
  // Children list for volume n: Indicies of volumes that are entirely within the set described by volume n
  // Overlap list for volume n: Indicies of volume that contains some of the set described by volume n (excluding volume n)
  // Intersect check list for volume n: Indicies of volumes to check for intersection if a ray originates from volume n (is generated from the children and overlap lists)
  // Parents list for volume n: Indicies of volumes that contain the entire set of volume n
  // Grandparents lists for volume n: Indicies of volumes that contain the entire set of at least one parent of volume n
  // Destination list for volume n: Indicies of volumes that could be the destination volume when a ray leaves volume n
  // The overlap, parents and grandparents lists are local variables in the function, and not in the main scope.
  
  generate_lists(Volumes, &starting_lists, number_of_volumes, list_verbal);
  
  // Generate "safe starting list", which contains all volumes that the ray may enter from other components
  // These are all volumes without scattering or absorption
  
  for (volume_index=0;volume_index<number_of_volumes;volume_index++) {
    // Updating mask lists from volume index to global_mask_indices
    // Filling out the masked by list that uses mask indices
    Volumes[volume_index]->geometry.masked_by_mask_index_list.num_elements = Volumes[volume_index]->geometry.masked_by_list.num_elements;
    Volumes[volume_index]->geometry.masked_by_mask_index_list.elements = malloc(Volumes[volume_index]->geometry.masked_by_mask_index_list.num_elements * sizeof(int));
    for (iterate=0;iterate<Volumes[volume_index]->geometry.masked_by_list.num_elements;iterate++)
        Volumes[volume_index]->geometry.masked_by_mask_index_list.elements[iterate] = find_on_int_list(mask_volume_index_list,Volumes[volume_index]->geometry.masked_by_list.elements[iterate]);
  }
  
  // Optimizing speed of the within_which_volume search algorithm
  int volume_index_main;
  /*
  for (volume_index_main=0;volume_index_main<number_of_volumes;volume_index_main++) {
    Volumes[volume_index_main]->geometry.destinations_logic_list.elements[0] = 0;
  }
  */
  
  // Checking for equal priorities
  for (volume_index_main=0;volume_index_main<number_of_volumes;volume_index_main++) {
    for (volume_index=0;volume_index<number_of_volumes;volume_index++)
        if (Volumes[volume_index_main]->geometry.priority_value == Volumes[volume_index]->geometry.priority_value && volume_index_main != volume_index) {
            if (Volumes[volume_index_main]->geometry.is_mask_volume == 0 && Volumes[volume_index]->geometry.is_mask_volume == 0) {
              // Priority of masks do not matter
              printf("ERROR in Union_master with name %s. The volumes named %s and %s have the same priority. Change the priorities so the one present in case of overlap has highest priority.\n",NAME_CURRENT_COMP,Volumes[volume_index_main]->name,Volumes[volume_index]->name);
              exit(1);
            }
        }
  }
  
  
  
  MPI_MASTER(
  if (verbal) printf("\n ---- Overview of the lists generated for each volume ---- \n");
  
  for (volume_index_main=0;volume_index_main<number_of_volumes;volume_index_main++) {
      if (verbal) sprintf(string_output,"Children for Volume                  %d",volume_index_main);
      if (verbal) print_1d_int_list(Volumes[volume_index_main]->geometry.children,string_output);
      
      if (verbal) sprintf(string_output,"Direct_children for Volume           %d",volume_index_main);
      if (verbal) print_1d_int_list(Volumes[volume_index_main]->geometry.direct_children,string_output);
      
      if (verbal) sprintf(string_output,"Intersect_check_list for Volume      %d",volume_index_main);
      if (verbal) print_1d_int_list(Volumes[volume_index_main]->geometry.intersect_check_list,string_output);
      
      if (verbal) sprintf(string_output,"Mask_intersect_list for Volume       %d",volume_index_main);
      if (verbal) print_1d_int_list(Volumes[volume_index_main]->geometry.mask_intersect_list,string_output);
      
      if (verbal) sprintf(string_output,"Destinations_list for Volume         %d",volume_index_main);
      if (verbal) print_1d_int_list(Volumes[volume_index_main]->geometry.destinations_list,string_output);
      
      //if (verbal) sprintf(string_output,"Destinations_logic_list for Volume   %d",volume_index_main);
      //if (verbal) print_1d_int_list(Volumes[volume_index_main]->geometry.destinations_logic_list,string_output);
      
      if (verbal) sprintf(string_output,"Reduced_destinations_list for Volume %d",volume_index_main);
      if (verbal) print_1d_int_list(Volumes[volume_index_main]->geometry.reduced_destinations_list,string_output);
      
      if (verbal) sprintf(string_output,"Next_volume_list for Volume          %d",volume_index_main);
      if (verbal) print_1d_int_list(Volumes[volume_index_main]->geometry.next_volume_list,string_output);
      
      if (verbal) {
        if (volume_index_main != 0)
                           printf("      Is_vacuum for Volume                 %d = %d\n",volume_index_main,Volumes[volume_index_main]->p_physics->is_vacuum);
      }
      if (verbal) {
        if (volume_index_main != 0)
                           printf("      is_mask_volume for Volume            %d = %d\n",volume_index_main,Volumes[volume_index_main]->geometry.is_mask_volume);
      }
      if (verbal) {
        if (volume_index_main != 0)
                           printf("      is_masked_volume for Volume          %d = %d\n",volume_index_main,Volumes[volume_index_main]->geometry.is_masked_volume);
      }
      if (verbal) {
        if (volume_index_main != 0)
                           printf("      is_exit_volume for Volume            %d = %d\n",volume_index_main,Volumes[volume_index_main]->geometry.is_exit_volume);
      }
      
      if (verbal) sprintf(string_output,"mask_list for Volume                 %d",volume_index_main);
      if (verbal) print_1d_int_list(Volumes[volume_index_main]->geometry.mask_list,string_output);
      
      if (verbal) sprintf(string_output,"masked_by_list for Volume            %d",volume_index_main);
      if (verbal) print_1d_int_list(Volumes[volume_index_main]->geometry.masked_by_list,string_output);
      
      if (verbal) sprintf(string_output,"masked_by_mask_index_list for Volume %d",volume_index_main);
      if (verbal) print_1d_int_list(Volumes[volume_index_main]->geometry.masked_by_mask_index_list,string_output);
      
      if (verbal)          printf("      mask_mode for Volume                 %d = %d\n",volume_index_main,Volumes[volume_index_main]->geometry.mask_mode);
      if (verbal) printf("\n");
  }
  )
  


  // Initializing intersection_time_table
  // The intersection time table contains all information on intersection times for the current position/direction, and is cleared everytime a ray changes direction.
  // Not all entries needs to be calculated, so there is a variable that keeps track of which intersection times have been calculated in order to avoid redoing that.
  // When the intersections times are calculated for a volume, all future intersections are kept in the time table.
  // Thus the memory allocation have to take into account how many intersections there can be with each volume, but it is currently set to 2, but can easily be changed.
  
  intersection_time_table.num_volumes = number_of_volumes;
  
  intersection_time_table.n_elements = (int*) malloc(intersection_time_table.num_volumes * sizeof(int));
  intersection_time_table.calculated = (int*) malloc(intersection_time_table.num_volumes * sizeof(int));
  intersection_time_table.intersection_times = (double**) malloc(intersection_time_table.num_volumes * sizeof(double));
  for (iterate = 0;iterate < intersection_time_table.num_volumes;iterate++){
      intersection_time_table.n_elements[iterate] = (int) 2; // number of intersection solutions
      if (iterate == 0) intersection_time_table.n_elements[iterate] = (int) 0; // number of intersection solutions
      intersection_time_table.calculated[iterate] = (int) 0; // Initializing calculated logic
      
      if (iterate == 0) {
           intersection_time_table.intersection_times[0] = NULL;
      }
      else {
        intersection_time_table.intersection_times[iterate] = (double*) malloc(intersection_time_table.n_elements[iterate]*sizeof(double));
        for (solutions = 0;solutions < intersection_time_table.n_elements[iterate];solutions++) {
              intersection_time_table.intersection_times[iterate][solutions] = -1.0;
        }
      }
  }
  
  
  
  
  // If enabled, the tagging system tracks all different histories sampled by the program.

  // Initialize the tagging tree
  // Allocate a list of host nodes with the same length as the number of volumes
  
  stop_creating_nodes = 0; stop_tagging_ray = 0; tagging_leaf_counter = 0;
  if (enable_tagging) {
    master_tagging_node_list.num_elements = number_of_volumes;
    master_tagging_node_list.elements = malloc(master_tagging_node_list.num_elements * sizeof(struct tagging_tree_node_struct*));
  
    // Initialize
    for (volume_index=0;volume_index<number_of_volumes;volume_index++) {
      //if (verbal) printf("Allocating master tagging node for volume number %d \n",volume_index);
      master_tagging_node_list.elements[volume_index] = initialize_tagging_tree_node(master_tagging_node_list.elements[volume_index], NULL, Volumes[volume_index]);
      //if (verbal) printf("Allocated master tagging node for volume number %d \n",volume_index);
    }
  }
  
  // Initialize loggers
  loggers_with_data_array.allocated_elements = 0;
  loggers_with_data_array.used_elements = 0;
  
  
  // Signal initialization complete
  MPI_MASTER(
  printf("Union_master component %s initialized sucessfully\n",NAME_CURRENT_COMP);
  )

}
#line 15670 "Union_manual_example.c"
#undef history_limit
#undef allow_inside_start
#undef number_of_processes_array
#undef safty_distance2
#undef safty_distance
#undef max_conditional_extend_index
#undef tagging_conditional_extend
#undef logger_conditional_extend_array
#undef free_tagging_conditioanl_list
#undef tagging_conditional_list
#undef conditional_status
#undef this_logger
#undef p_old
#undef loggers_with_data_array
#undef Volume_copies_allocated
#undef geometry_component_index_list
#undef mask_volume_index_list
#undef current_mask_intersect_list_status
#undef mask_status_list
#undef mask_iterate
#undef mask_index_main
#undef need_to_run_within_which_volume
#undef number_of_masked_volumes
#undef number_of_masks
#undef mc_transmission_probability
#undef real_transmission_probability
#undef number_of_scattering_events
#undef tagging_leaf_counter
#undef current_tagging_node
#undef master_tagging_node_list
#undef enable_tagging_check
#undef stop_creating_nodes
#undef stop_tagging_ray
#undef enable_tagging
#undef rotated_position
#undef non_rotated_position
#undef temp_rotation_matrix
#undef master_transposed_rotation_matrix
#undef scattered_flag_VP
#undef scattered_flag
#undef volume_0_found
#undef ray_velocity_final
#undef ray_velocity
#undef ray_position
#undef pre_allocated3
#undef pre_allocated2
#undef pre_allocated1
#undef tree_next_volume
#undef geometry_output
#undef intersection_with_children
#undef start
#undef check
#undef number_of_solutions_static
#undef number_of_solutions
#undef current_volume
#undef done
#undef next_volume_priority
#undef next_volume
#undef a_next_volume_found
#undef time_propagated_without_scattering
#undef scattering_event
#undef selected_process
#undef time_to_boundery
#undef length_to_boundery
#undef length_to_scattering
#undef time_to_scattering
#undef mc_prop
#undef culmative_probability
#undef my_sum_plus_abs
#undef my_sum
#undef v_length
#undef k_new
#undef k
#undef my_trace_fraction_control
#undef p_my_trace
#undef my_trace
#undef process_start
#undef process
#undef min_intersection_time
#undef intersection_time
#undef time_found
#undef min_volume
#undef min_solution
#undef solution
#undef limit
#undef max_number_of_processes
#undef solutions
#undef process_index
#undef volume_index
#undef number_of_volumes
#undef string_output
#undef component_error_msg
#undef error_msg
#undef v
#undef r_start
#undef r
#undef starting_lists
#undef Geometries
#undef Volumes
#undef intersection_time_table
#undef geometry_list_index
#undef previous_master_index
#undef global_master_element
#undef starting_volume_warning
#undef finally_verbal
#undef trace_verbal
#undef list_verbal
#undef verbal
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Banana_monitor. */
  SIG_MESSAGE("Banana_monitor (Init)");
#define mccompcurname  Banana_monitor
#define mccompcurtype  Monitor_nD
#define mccompcurindex 15
#define user1 mccBanana_monitor_user1
#define user2 mccBanana_monitor_user2
#define user3 mccBanana_monitor_user3
#define DEFS mccBanana_monitor_DEFS
#define Vars mccBanana_monitor_Vars
#define detector mccBanana_monitor_detector
#define offdata mccBanana_monitor_offdata
#define xwidth mccBanana_monitor_xwidth
#define yheight mccBanana_monitor_yheight
#define zdepth mccBanana_monitor_zdepth
#define xmin mccBanana_monitor_xmin
#define xmax mccBanana_monitor_xmax
#define ymin mccBanana_monitor_ymin
#define ymax mccBanana_monitor_ymax
#define zmin mccBanana_monitor_zmin
#define zmax mccBanana_monitor_zmax
#define bins mccBanana_monitor_bins
#define min mccBanana_monitor_min
#define max mccBanana_monitor_max
#define restore_neutron mccBanana_monitor_restore_neutron
#define radius mccBanana_monitor_radius
#define options mccBanana_monitor_options
#define filename mccBanana_monitor_filename
#define geometry mccBanana_monitor_geometry
#define username1 mccBanana_monitor_username1
#define username2 mccBanana_monitor_username2
#define username3 mccBanana_monitor_username3
#define nowritefile mccBanana_monitor_nowritefile
#line 230 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
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
  if (geometry && strlen(geometry) && strcmp(geometry,"0") && strcmp(geometry, "NULL"))
    if (!off_init(  geometry, xwidth, yheight, zdepth, 1, &offdata )) {
      printf("Monitor_nD: %s could not initiate the OFF geometry %s. \n"
             "            Defaulting to normal Monitor dimensions.\n",
             NAME_CURRENT_COMP, geometry);
      strcpy(geometry, "");
    } else {
      offflag=1;
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
#line 15895 "Union_manual_example.c"
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
  /* TRACE Component Al_inc [1] */
  mccoordschange(mcposrAl_inc, mcrotrAl_inc,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Al_inc (without coords transformations) */
  mcJumpTrace_Al_inc:
  SIG_MESSAGE("Al_inc (Trace)");
  mcDEBUG_COMP("Al_inc")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompAl_inc
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
#define mccompcurname  Al_inc
#define mccompcurtype  Incoherent_process
#define mccompcurindex 1
#define This_process mccAl_inc_This_process
#define Incoherent_storage mccAl_inc_Incoherent_storage
#define effective_my_scattering mccAl_inc_effective_my_scattering
{   /* Declarations of Al_inc=Incoherent_process() SETTING parameters. */
MCNUM sigma = mccAl_inc_sigma;
MCNUM packing_factor = mccAl_inc_packing_factor;
MCNUM unit_cell_volume = mccAl_inc_unit_cell_volume;
MCNUM interact_fraction = mccAl_inc_interact_fraction;
}   /* End of Al_inc=Incoherent_process() SETTING parameter declarations. */
#undef effective_my_scattering
#undef Incoherent_storage
#undef This_process
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompAl_inc:
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

  /* TRACE Component Al_powder [2] */
  mccoordschange(mcposrAl_powder, mcrotrAl_powder,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Al_powder (without coords transformations) */
  mcJumpTrace_Al_powder:
  SIG_MESSAGE("Al_powder (Trace)");
  mcDEBUG_COMP("Al_powder")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompAl_powder
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
#define mccompcurname  Al_powder
#define mccompcurtype  Powder_process
#define mccompcurindex 2
#define format mccAl_powder_format
#define This_process mccAl_powder_This_process
#define Powder_storage mccAl_powder_Powder_storage
#define effective_my_scattering mccAl_powder_effective_my_scattering
#define d_phi mccAl_powder_d_phi
#define line_info mccAl_powder_line_info
#define columns mccAl_powder_columns
{   /* Declarations of Al_powder=Powder_process() SETTING parameters. */
char* reflections = mccAl_powder_reflections;
MCNUM packing_factor = mccAl_powder_packing_factor;
MCNUM Vc = mccAl_powder_Vc;
MCNUM delta_d_d = mccAl_powder_delta_d_d;
MCNUM DW = mccAl_powder_DW;
MCNUM nb_atoms = mccAl_powder_nb_atoms;
MCNUM d_phi = mccAl_powder_d_phi;
MCNUM density = mccAl_powder_density;
MCNUM weight = mccAl_powder_weight;
MCNUM barns = mccAl_powder_barns;
MCNUM Strain = mccAl_powder_Strain;
MCNUM interact_fraction = mccAl_powder_interact_fraction;
}   /* End of Al_powder=Powder_process() SETTING parameter declarations. */
#undef columns
#undef line_info
#undef d_phi
#undef effective_my_scattering
#undef Powder_storage
#undef This_process
#undef format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompAl_powder:
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

  /* TRACE Component Al [3] */
  mccoordschange(mcposrAl, mcrotrAl,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Al (without coords transformations) */
  mcJumpTrace_Al:
  SIG_MESSAGE("Al (Trace)");
  mcDEBUG_COMP("Al")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompAl
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
#define mccompcurname  Al
#define mccompcurtype  Union_make_material
#define mccompcurindex 3
#define loop_index mccAl_loop_index
#define this_material mccAl_this_material
#define accepted_processes mccAl_accepted_processes
#define global_material_element mccAl_global_material_element
{   /* Declarations of Al=Union_make_material() SETTING parameters. */
char* process_string = mccAl_process_string;
MCNUM my_absorption = mccAl_my_absorption;
MCNUM absorber = mccAl_absorber;
}   /* End of Al=Union_make_material() SETTING parameter declarations. */
#undef global_material_element
#undef accepted_processes
#undef this_material
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompAl:
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

  /* TRACE Component Fe_inc [4] */
  mccoordschange(mcposrFe_inc, mcrotrFe_inc,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Fe_inc (without coords transformations) */
  mcJumpTrace_Fe_inc:
  SIG_MESSAGE("Fe_inc (Trace)");
  mcDEBUG_COMP("Fe_inc")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompFe_inc
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
#define mccompcurname  Fe_inc
#define mccompcurtype  Incoherent_process
#define mccompcurindex 4
#define This_process mccFe_inc_This_process
#define Incoherent_storage mccFe_inc_Incoherent_storage
#define effective_my_scattering mccFe_inc_effective_my_scattering
{   /* Declarations of Fe_inc=Incoherent_process() SETTING parameters. */
MCNUM sigma = mccFe_inc_sigma;
MCNUM packing_factor = mccFe_inc_packing_factor;
MCNUM unit_cell_volume = mccFe_inc_unit_cell_volume;
MCNUM interact_fraction = mccFe_inc_interact_fraction;
}   /* End of Fe_inc=Incoherent_process() SETTING parameter declarations. */
#undef effective_my_scattering
#undef Incoherent_storage
#undef This_process
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompFe_inc:
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

  /* TRACE Component Fe_powder [5] */
  mccoordschange(mcposrFe_powder, mcrotrFe_powder,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Fe_powder (without coords transformations) */
  mcJumpTrace_Fe_powder:
  SIG_MESSAGE("Fe_powder (Trace)");
  mcDEBUG_COMP("Fe_powder")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompFe_powder
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
#define mccompcurname  Fe_powder
#define mccompcurtype  Powder_process
#define mccompcurindex 5
#define format mccFe_powder_format
#define This_process mccFe_powder_This_process
#define Powder_storage mccFe_powder_Powder_storage
#define effective_my_scattering mccFe_powder_effective_my_scattering
#define d_phi mccFe_powder_d_phi
#define line_info mccFe_powder_line_info
#define columns mccFe_powder_columns
{   /* Declarations of Fe_powder=Powder_process() SETTING parameters. */
char* reflections = mccFe_powder_reflections;
MCNUM packing_factor = mccFe_powder_packing_factor;
MCNUM Vc = mccFe_powder_Vc;
MCNUM delta_d_d = mccFe_powder_delta_d_d;
MCNUM DW = mccFe_powder_DW;
MCNUM nb_atoms = mccFe_powder_nb_atoms;
MCNUM d_phi = mccFe_powder_d_phi;
MCNUM density = mccFe_powder_density;
MCNUM weight = mccFe_powder_weight;
MCNUM barns = mccFe_powder_barns;
MCNUM Strain = mccFe_powder_Strain;
MCNUM interact_fraction = mccFe_powder_interact_fraction;
}   /* End of Fe_powder=Powder_process() SETTING parameter declarations. */
#undef columns
#undef line_info
#undef d_phi
#undef effective_my_scattering
#undef Powder_storage
#undef This_process
#undef format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompFe_powder:
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

  /* TRACE Component Fe [6] */
  mccoordschange(mcposrFe, mcrotrFe,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Fe (without coords transformations) */
  mcJumpTrace_Fe:
  SIG_MESSAGE("Fe (Trace)");
  mcDEBUG_COMP("Fe")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompFe
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
#define mccompcurname  Fe
#define mccompcurtype  Union_make_material
#define mccompcurindex 6
#define loop_index mccFe_loop_index
#define this_material mccFe_this_material
#define accepted_processes mccFe_accepted_processes
#define global_material_element mccFe_global_material_element
{   /* Declarations of Fe=Union_make_material() SETTING parameters. */
char* process_string = mccFe_process_string;
MCNUM my_absorption = mccFe_my_absorption;
MCNUM absorber = mccFe_absorber;
}   /* End of Fe=Union_make_material() SETTING parameter declarations. */
#undef global_material_element
#undef accepted_processes
#undef this_material
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompFe:
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

  /* TRACE Component Progress [7] */
  mccoordschange(mcposrProgress, mcrotrProgress,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Progress (without coords transformations) */
  mcJumpTrace_Progress:
  SIG_MESSAGE("Progress (Trace)");
  mcDEBUG_COMP("Progress")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompProgress
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
#define mccompcurname  Progress
#define mccompcurtype  Progress_bar
#define mccompcurindex 7
#define IntermediateCnts mccProgress_IntermediateCnts
#define StartTime mccProgress_StartTime
#define EndTime mccProgress_EndTime
#define CurrentTime mccProgress_CurrentTime
{   /* Declarations of Progress=Progress_bar() SETTING parameters. */
char* profile = mccProgress_profile;
MCNUM percent = mccProgress_percent;
MCNUM flag_save = mccProgress_flag_save;
MCNUM minutes = mccProgress_minutes;
#line 70 "/usr/share/mcstas/2.5/misc/Progress_bar.comp"
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
#line 16805 "Union_manual_example.c"
}   /* End of Progress=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompProgress:
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

  /* TRACE Component Source [8] */
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
#define mccompcurname  Source
#define mccompcurtype  Source_simple
#define mccompcurindex 8
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
#line 125 "/usr/share/mcstas/2.5/sources/Source_simple.comp"
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
#line 16976 "Union_manual_example.c"
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

  /* TRACE Component Guide [9] */
  mccoordschange(mcposrGuide, mcrotrGuide,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Guide (without coords transformations) */
  mcJumpTrace_Guide:
  SIG_MESSAGE("Guide (Trace)");
  mcDEBUG_COMP("Guide")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompGuide
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
#define mccompcurname  Guide
#define mccompcurtype  Guide_gravity
#define mccompcurindex 9
#define GVars mccGuide_GVars
#define pTable mccGuide_pTable
{   /* Declarations of Guide=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccGuide_w1;
MCNUM h1 = mccGuide_h1;
MCNUM w2 = mccGuide_w2;
MCNUM h2 = mccGuide_h2;
MCNUM l = mccGuide_l;
MCNUM R0 = mccGuide_R0;
MCNUM Qc = mccGuide_Qc;
MCNUM alpha = mccGuide_alpha;
MCNUM m = mccGuide_m;
MCNUM W = mccGuide_W;
MCNUM nslit = mccGuide_nslit;
MCNUM d = mccGuide_d;
MCNUM mleft = mccGuide_mleft;
MCNUM mright = mccGuide_mright;
MCNUM mtop = mccGuide_mtop;
MCNUM mbottom = mccGuide_mbottom;
MCNUM nhslit = mccGuide_nhslit;
MCNUM G = mccGuide_G;
MCNUM aleft = mccGuide_aleft;
MCNUM aright = mccGuide_aright;
MCNUM atop = mccGuide_atop;
MCNUM abottom = mccGuide_abottom;
MCNUM wavy = mccGuide_wavy;
MCNUM wavy_z = mccGuide_wavy_z;
MCNUM wavy_tb = mccGuide_wavy_tb;
MCNUM wavy_lr = mccGuide_wavy_lr;
MCNUM chamfers = mccGuide_chamfers;
MCNUM chamfers_z = mccGuide_chamfers_z;
MCNUM chamfers_lr = mccGuide_chamfers_lr;
MCNUM chamfers_tb = mccGuide_chamfers_tb;
MCNUM nelements = mccGuide_nelements;
MCNUM nu = mccGuide_nu;
MCNUM phase = mccGuide_phase;
char* reflect = mccGuide_reflect;
#line 392 "/usr/share/mcstas/2.5/optics/Guide_gravity.comp"
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
#line 17290 "Union_manual_example.c"
}   /* End of Guide=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompGuide:
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

  /* TRACE Component sample_position [10] */
  mccoordschange(mcposrsample_position, mcrotrsample_position,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component sample_position (without coords transformations) */
  mcJumpTrace_sample_position:
  SIG_MESSAGE("sample_position (Trace)");
  mcDEBUG_COMP("sample_position")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompsample_position
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
#define mccompcurname  sample_position
#define mccompcurtype  Arm
#define mccompcurindex 10
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompsample_position:
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

  /* TRACE Component Hexagonal_container_1 [11] */
  mccoordschange(mcposrHexagonal_container_1, mcrotrHexagonal_container_1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Hexagonal_container_1 (without coords transformations) */
  mcJumpTrace_Hexagonal_container_1:
  SIG_MESSAGE("Hexagonal_container_1 (Trace)");
  mcDEBUG_COMP("Hexagonal_container_1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompHexagonal_container_1
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
#define mccompcurname  Hexagonal_container_1
#define mccompcurtype  Union_box
#define mccompcurindex 11
#define loop_index mccHexagonal_container_1_loop_index
#define this_box_volume mccHexagonal_container_1_this_box_volume
#define global_geometry_element mccHexagonal_container_1_global_geometry_element
#define this_box_storage mccHexagonal_container_1_this_box_storage
{   /* Declarations of Hexagonal_container_1=Union_box() SETTING parameters. */
char* material_string = mccHexagonal_container_1_material_string;
MCNUM priority = mccHexagonal_container_1_priority;
MCNUM xwidth = mccHexagonal_container_1_xwidth;
MCNUM yheight = mccHexagonal_container_1_yheight;
MCNUM zdepth = mccHexagonal_container_1_zdepth;
MCNUM xwidth2 = mccHexagonal_container_1_xwidth2;
MCNUM yheight2 = mccHexagonal_container_1_yheight2;
MCNUM visualize = mccHexagonal_container_1_visualize;
int target_index = mccHexagonal_container_1_target_index;
MCNUM target_x = mccHexagonal_container_1_target_x;
MCNUM target_y = mccHexagonal_container_1_target_y;
MCNUM target_z = mccHexagonal_container_1_target_z;
MCNUM focus_aw = mccHexagonal_container_1_focus_aw;
MCNUM focus_ah = mccHexagonal_container_1_focus_ah;
MCNUM focus_xw = mccHexagonal_container_1_focus_xw;
MCNUM focus_xh = mccHexagonal_container_1_focus_xh;
MCNUM focus_r = mccHexagonal_container_1_focus_r;
MCNUM p_interact = mccHexagonal_container_1_p_interact;
char* mask_string = mccHexagonal_container_1_mask_string;
char* mask_setting = mccHexagonal_container_1_mask_setting;
MCNUM number_of_activations = mccHexagonal_container_1_number_of_activations;
}   /* End of Hexagonal_container_1=Union_box() SETTING parameter declarations. */
#undef this_box_storage
#undef global_geometry_element
#undef this_box_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompHexagonal_container_1:
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

  /* TRACE Component Hexagonal_container_2 [12] */
  mccoordschange(mcposrHexagonal_container_2, mcrotrHexagonal_container_2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Hexagonal_container_2 (without coords transformations) */
  mcJumpTrace_Hexagonal_container_2:
  SIG_MESSAGE("Hexagonal_container_2 (Trace)");
  mcDEBUG_COMP("Hexagonal_container_2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompHexagonal_container_2
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
#define mccompcurname  Hexagonal_container_2
#define mccompcurtype  Union_box
#define mccompcurindex 12
#define loop_index mccHexagonal_container_2_loop_index
#define this_box_volume mccHexagonal_container_2_this_box_volume
#define global_geometry_element mccHexagonal_container_2_global_geometry_element
#define this_box_storage mccHexagonal_container_2_this_box_storage
{   /* Declarations of Hexagonal_container_2=Union_box() SETTING parameters. */
char* material_string = mccHexagonal_container_2_material_string;
MCNUM priority = mccHexagonal_container_2_priority;
MCNUM xwidth = mccHexagonal_container_2_xwidth;
MCNUM yheight = mccHexagonal_container_2_yheight;
MCNUM zdepth = mccHexagonal_container_2_zdepth;
MCNUM xwidth2 = mccHexagonal_container_2_xwidth2;
MCNUM yheight2 = mccHexagonal_container_2_yheight2;
MCNUM visualize = mccHexagonal_container_2_visualize;
int target_index = mccHexagonal_container_2_target_index;
MCNUM target_x = mccHexagonal_container_2_target_x;
MCNUM target_y = mccHexagonal_container_2_target_y;
MCNUM target_z = mccHexagonal_container_2_target_z;
MCNUM focus_aw = mccHexagonal_container_2_focus_aw;
MCNUM focus_ah = mccHexagonal_container_2_focus_ah;
MCNUM focus_xw = mccHexagonal_container_2_focus_xw;
MCNUM focus_xh = mccHexagonal_container_2_focus_xh;
MCNUM focus_r = mccHexagonal_container_2_focus_r;
MCNUM p_interact = mccHexagonal_container_2_p_interact;
char* mask_string = mccHexagonal_container_2_mask_string;
char* mask_setting = mccHexagonal_container_2_mask_setting;
MCNUM number_of_activations = mccHexagonal_container_2_number_of_activations;
}   /* End of Hexagonal_container_2=Union_box() SETTING parameter declarations. */
#undef this_box_storage
#undef global_geometry_element
#undef this_box_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompHexagonal_container_2:
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

  /* TRACE Component Sample [13] */
  mccoordschange(mcposrSample, mcrotrSample,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Sample (without coords transformations) */
  mcJumpTrace_Sample:
  SIG_MESSAGE("Sample (Trace)");
  mcDEBUG_COMP("Sample")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompSample
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
#define mccompcurname  Sample
#define mccompcurtype  Union_cylinder
#define mccompcurindex 13
#define loop_index mccSample_loop_index
#define this_cylinder_volume mccSample_this_cylinder_volume
#define global_geometry_element mccSample_global_geometry_element
#define this_cylinder_storage mccSample_this_cylinder_storage
{   /* Declarations of Sample=Union_cylinder() SETTING parameters. */
char* material_string = mccSample_material_string;
MCNUM priority = mccSample_priority;
MCNUM radius = mccSample_radius;
MCNUM yheight = mccSample_yheight;
MCNUM visualize = mccSample_visualize;
int target_index = mccSample_target_index;
MCNUM target_x = mccSample_target_x;
MCNUM target_y = mccSample_target_y;
MCNUM target_z = mccSample_target_z;
MCNUM focus_aw = mccSample_focus_aw;
MCNUM focus_ah = mccSample_focus_ah;
MCNUM focus_xw = mccSample_focus_xw;
MCNUM focus_xh = mccSample_focus_xh;
MCNUM focus_r = mccSample_focus_r;
MCNUM p_interact = mccSample_p_interact;
char* mask_string = mccSample_mask_string;
char* mask_setting = mccSample_mask_setting;
MCNUM number_of_activations = mccSample_number_of_activations;
#line 344 "/usr/share/mcstas/2.5/contrib/union/Union_cylinder.comp"
{
dummy = 1;
}
#line 17795 "Union_manual_example.c"
}   /* End of Sample=Union_cylinder() SETTING parameter declarations. */
#undef this_cylinder_storage
#undef global_geometry_element
#undef this_cylinder_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompSample:
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

  /* TRACE Component Master [14] */
  mccoordschange(mcposrMaster, mcrotrMaster,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Master (without coords transformations) */
  mcJumpTrace_Master:
  SIG_MESSAGE("Master (Trace)");
  mcDEBUG_COMP("Master")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompMaster
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
#define mccompcurname  Master
#define mccompcurtype  Union_master
#define mccompcurindex 14
#define verbal mccMaster_verbal
#define list_verbal mccMaster_list_verbal
#define trace_verbal mccMaster_trace_verbal
#define finally_verbal mccMaster_finally_verbal
#define starting_volume_warning mccMaster_starting_volume_warning
#define global_master_element mccMaster_global_master_element
#define previous_master_index mccMaster_previous_master_index
#define geometry_list_index mccMaster_geometry_list_index
#define intersection_time_table mccMaster_intersection_time_table
#define Volumes mccMaster_Volumes
#define Geometries mccMaster_Geometries
#define starting_lists mccMaster_starting_lists
#define r mccMaster_r
#define r_start mccMaster_r_start
#define v mccMaster_v
#define error_msg mccMaster_error_msg
#define component_error_msg mccMaster_component_error_msg
#define string_output mccMaster_string_output
#define number_of_volumes mccMaster_number_of_volumes
#define volume_index mccMaster_volume_index
#define process_index mccMaster_process_index
#define solutions mccMaster_solutions
#define max_number_of_processes mccMaster_max_number_of_processes
#define limit mccMaster_limit
#define solution mccMaster_solution
#define min_solution mccMaster_min_solution
#define min_volume mccMaster_min_volume
#define time_found mccMaster_time_found
#define intersection_time mccMaster_intersection_time
#define min_intersection_time mccMaster_min_intersection_time
#define process mccMaster_process
#define process_start mccMaster_process_start
#define my_trace mccMaster_my_trace
#define p_my_trace mccMaster_p_my_trace
#define my_trace_fraction_control mccMaster_my_trace_fraction_control
#define k mccMaster_k
#define k_new mccMaster_k_new
#define v_length mccMaster_v_length
#define my_sum mccMaster_my_sum
#define my_sum_plus_abs mccMaster_my_sum_plus_abs
#define culmative_probability mccMaster_culmative_probability
#define mc_prop mccMaster_mc_prop
#define time_to_scattering mccMaster_time_to_scattering
#define length_to_scattering mccMaster_length_to_scattering
#define length_to_boundery mccMaster_length_to_boundery
#define time_to_boundery mccMaster_time_to_boundery
#define selected_process mccMaster_selected_process
#define scattering_event mccMaster_scattering_event
#define time_propagated_without_scattering mccMaster_time_propagated_without_scattering
#define a_next_volume_found mccMaster_a_next_volume_found
#define next_volume mccMaster_next_volume
#define next_volume_priority mccMaster_next_volume_priority
#define done mccMaster_done
#define current_volume mccMaster_current_volume
#define number_of_solutions mccMaster_number_of_solutions
#define number_of_solutions_static mccMaster_number_of_solutions_static
#define check mccMaster_check
#define start mccMaster_start
#define intersection_with_children mccMaster_intersection_with_children
#define geometry_output mccMaster_geometry_output
#define tree_next_volume mccMaster_tree_next_volume
#define pre_allocated1 mccMaster_pre_allocated1
#define pre_allocated2 mccMaster_pre_allocated2
#define pre_allocated3 mccMaster_pre_allocated3
#define ray_position mccMaster_ray_position
#define ray_velocity mccMaster_ray_velocity
#define ray_velocity_final mccMaster_ray_velocity_final
#define volume_0_found mccMaster_volume_0_found
#define scattered_flag mccMaster_scattered_flag
#define scattered_flag_VP mccMaster_scattered_flag_VP
#define master_transposed_rotation_matrix mccMaster_master_transposed_rotation_matrix
#define temp_rotation_matrix mccMaster_temp_rotation_matrix
#define non_rotated_position mccMaster_non_rotated_position
#define rotated_position mccMaster_rotated_position
#define enable_tagging mccMaster_enable_tagging
#define stop_tagging_ray mccMaster_stop_tagging_ray
#define stop_creating_nodes mccMaster_stop_creating_nodes
#define enable_tagging_check mccMaster_enable_tagging_check
#define master_tagging_node_list mccMaster_master_tagging_node_list
#define current_tagging_node mccMaster_current_tagging_node
#define tagging_leaf_counter mccMaster_tagging_leaf_counter
#define number_of_scattering_events mccMaster_number_of_scattering_events
#define real_transmission_probability mccMaster_real_transmission_probability
#define mc_transmission_probability mccMaster_mc_transmission_probability
#define number_of_masks mccMaster_number_of_masks
#define number_of_masked_volumes mccMaster_number_of_masked_volumes
#define need_to_run_within_which_volume mccMaster_need_to_run_within_which_volume
#define mask_index_main mccMaster_mask_index_main
#define mask_iterate mccMaster_mask_iterate
#define mask_status_list mccMaster_mask_status_list
#define current_mask_intersect_list_status mccMaster_current_mask_intersect_list_status
#define mask_volume_index_list mccMaster_mask_volume_index_list
#define geometry_component_index_list mccMaster_geometry_component_index_list
#define Volume_copies_allocated mccMaster_Volume_copies_allocated
#define loggers_with_data_array mccMaster_loggers_with_data_array
#define p_old mccMaster_p_old
#define this_logger mccMaster_this_logger
#define conditional_status mccMaster_conditional_status
#define tagging_conditional_list mccMaster_tagging_conditional_list
#define free_tagging_conditioanl_list mccMaster_free_tagging_conditioanl_list
#define logger_conditional_extend_array mccMaster_logger_conditional_extend_array
#define tagging_conditional_extend mccMaster_tagging_conditional_extend
#define max_conditional_extend_index mccMaster_max_conditional_extend_index
#define safty_distance mccMaster_safty_distance
#define safty_distance2 mccMaster_safty_distance2
#define number_of_processes_array mccMaster_number_of_processes_array
{   /* Declarations of Master=Union_master() SETTING parameters. */
MCNUM allow_inside_start = mccMaster_allow_inside_start;
MCNUM history_limit = mccMaster_history_limit;
#line 896 "/usr/share/mcstas/2.5/contrib/union/Union_master.comp"
{
  
  if (trace_verbal) printf("\n\n\n\n\n----------- NEW RAY -------------------------------------------------\n");
  if (trace_verbal) printf("Union_master component name: %s \n \n",NAME_CURRENT_COMP);
  
  // Initialize logic
  done = 0;
  error_msg = 0;
  clear_intersection_table(&intersection_time_table);
  
  
  time_propagated_without_scattering = 0;
  v_length = sqrt(vx*vx+vy*vy+vz*vz);
  
  // Initialize logger system / Statistics
  number_of_scattering_events = 0;
  for (volume_index = 1;volume_index<number_of_volumes;volume_index++) { // No reason to update volume 0, as scattering doesn't happen there
    scattered_flag[volume_index] = 0;
    for (process_index=0;process_index<number_of_processes_array[volume_index];process_index++)
      scattered_flag_VP[volume_index][process_index] = 0;
  }
  loggers_with_data_array.used_elements = 0;
  tagging_conditional_extend = 0;
  for (iterate=0;iterate<max_conditional_extend_index+1;iterate++) {
    logger_conditional_extend_array[iterate] = 0;
  }

  
  r_start[0] = x; r_start[1] = y; r_start[2] = z;
  r[0]=x;r[1]=y;r[2]=z;v[0]=vx;v[1]=vy;v[2]=vz;
  
  ray_position = coords_set(x,y,z);
  ray_velocity = coords_set(vx,vy,vz);
  
  // Mask update: need to check the mask status for the initial position
  for (iterate=0;iterate<number_of_masks;iterate++) {
    if(Volumes[mask_volume_index_list.elements[iterate]]->geometry.within_function(ray_position,&Volumes[mask_volume_index_list.elements[iterate]]->geometry) == 1) {
      mask_status_list.elements[iterate] = 1;
    } else {
      mask_status_list.elements[iterate] = 0;
    }
  }
      
  if (trace_verbal) print_1d_int_list(mask_status_list,"Initial mask status list");
  
  // Now the initial current_volume can be found, which requires the up to date mask_status_list
  current_volume = within_which_volume(ray_position,starting_lists.reduced_start_list,starting_lists.starting_destinations_list,Volumes,&mask_status_list,number_of_volumes,pre_allocated1,pre_allocated2,pre_allocated3);
  
  // Using the mask_status_list and the current volume, the current_mask_intersect_list_status can be made
  //  it contains the effective mask status of all volumes on the current volumes mask intersect list, which needs to be calculated,
  //  but only when the current volume or mask status changes, not under for example scattering inside the current volume
  update_current_mask_intersect_status(&current_mask_intersect_list_status, &mask_status_list, Volumes, &current_volume);
  
  if (trace_verbal) printf("Starting current_volume = %d\n",current_volume);
  
  // This error is commented out for testing external components in Union systems
  if (allow_inside_start == 0 && starting_lists.allowed_starting_volume_logic_list.elements[current_volume] == 0) {
    printf("ERROR, ray ''teleported'' into Union component %s, if intentional, set allow_inside_start=1\n",NAME_CURRENT_COMP);
  }
  if (starting_volume_warning == 0 && current_volume != 0) {
        printf("WARNING: Ray started in volume ''%s'' rather than the surrounding vacuum in component %s. This warning is only shown once.\n",Volumes[current_volume]->name,NAME_CURRENT_COMP);
        starting_volume_warning = 1;
  }
  
  // Starting the tagging tree
  if (enable_tagging) {
    current_tagging_node = master_tagging_node_list.elements[current_volume];
    if (tagging_leaf_counter > history_limit) stop_creating_nodes = 1;
      stop_tagging_ray = 0;
  }
  
  if (trace_verbal && enable_tagging) printf("current_tagging_node->intensity = %f\n",current_tagging_node->intensity);
  if (trace_verbal && enable_tagging) printf("current_tagging_node->number_of_rays = %d \n",current_tagging_node->number_of_rays);
  
  // Multiple scattering loop
  limit = 100000;
  while (done == 0) {
    limit--;
    
    /*
    if (limit < 3) trace_verbal = 1;
    else trace_verbal = 0;
    */
    
    if (trace_verbal) printf("----------- START OF WHILE LOOP --------------------------------------\n");
    if (trace_verbal) print_intersection_table(intersection_time_table);
    if (trace_verbal) printf("current_volume = %d \n",current_volume);
    
    
    /*
    The first part of the main loop is ripe for being rewritten, as many lines of code are repated, and redundant work is being performed.
    The current structure is as follows:
    
    Calculation of intersections for volumes on the intersection list
    Calculation of intersections for volumes on active mask list (including work on finding if ALL/ANY is respected for the mask status)
    Calculation of intersections for the current volume, if no time no time have been found for one of it's children
    Selection of the lowest time among intersection list
    Selection of the lowest time among active intersection mask lists
    Selection of the lowest time for the current volume
    
    The most obvious change is that the time selection can be done together with the calculations.
    It is kept for now as the program needs to work as earlier when doing testing of new features.
    */
    
    // Checking intersections for all volumes in the intersect list.
    for (start=check=Volumes[current_volume]->geometry.intersect_check_list.elements;check-start<Volumes[current_volume]->geometry.intersect_check_list.num_elements;check++) {
    // This will leave check as a pointer to the intergers in the intersect_check_list and iccrement nicely
        if (trace_verbal) printf("Intersect_list = %d being checked \n",*check);
    
        if (intersection_time_table.calculated[*check] == 0) {
             if (trace_verbal) printf("running intersection for intersect_list with *check = %d \n",*check);
            // if (trace_verbal) printf("r = (%f,%f,%f) v = (%f,%f,%f) \n",r[0],r[1],r[2],v[0],v[1],v[2]);
            geometry_output = Volumes[*check]->geometry.intersect_function(intersection_time_table.intersection_times[*check],number_of_solutions,r_start,v,&Volumes[*check]->geometry);
            intersection_time_table.calculated[*check] = 1;
            // if (trace_verbal) printf("succesfully calculated intersection times for volume *check = %d \n",*check);
        }
    }
    
    // Mask update: add additional loop for checking intersections with masked volumes depending on mask statuses
    //for (start=check=Volumes[current_volume]->geometry.mask_intersect_list.elements;check-start<Volumes[current_volume]->geometry.mask_intersect_list.num_elements;check++) {
    // Can not use this method, as it is not certain the elements are even allocated
    for (mask_iterate=0;mask_iterate<Volumes[current_volume]->geometry.mask_intersect_list.num_elements;mask_iterate++) {
      if (current_mask_intersect_list_status.elements[mask_iterate] == 1) {
        if (trace_verbal) printf("Mask Intersect_list = %d being checked \n",Volumes[current_volume]->geometry.mask_intersect_list.elements[mask_iterate]);
        if (intersection_time_table.calculated[Volumes[current_volume]->geometry.mask_intersect_list.elements[mask_iterate]] == 0) {
          if (trace_verbal) printf("running intersection for mask_intersect_list element = %d \n",Volumes[current_volume]->geometry.mask_intersect_list.elements[mask_iterate]);
          // if (trace_verbal) printf("r = (%f,%f,%f) v = (%f,%f,%f) \n",r[0],r[1],r[2],v[0],v[1],v[2]);
          geometry_output = Volumes[Volumes[current_volume]->geometry.mask_intersect_list.elements[mask_iterate]]->geometry.intersect_function(intersection_time_table.intersection_times[Volumes[current_volume]->geometry.mask_intersect_list.elements[mask_iterate]],number_of_solutions,r_start,v,&Volumes[Volumes[current_volume]->geometry.mask_intersect_list.elements[mask_iterate]]->geometry);
          intersection_time_table.calculated[Volumes[current_volume]->geometry.mask_intersect_list.elements[mask_iterate]] = 1;
          // if (trace_verbal) printf("succesfully calculated intersection times for volume *check = %d \n",*check);
        }
      }
    }
    
    // removed this small optimization for debug purpose, and it fixed all problems
    // Checking if there are intersections with children of current volume, which means there is an intersection before current_volume, and thus can be skipped. But only if they have not been overwritten. In case current_volume is 0, there is no need to do this regardless.
    if (current_volume != 0 && intersection_time_table.calculated[current_volume] == 0) {
        if (trace_verbal) printf("Checking if children of current_volume = %d have intersections. \n",current_volume);
        intersection_with_children = 0;
        for (start = check = Volumes[current_volume]->geometry.children.elements;check - start < Volumes[current_volume]->geometry.children.num_elements;check++) {
            if (trace_verbal) printf("Checking if child %d of current_volume = %d have intersections. \n",*check,current_volume);
            // Not all children will be calculated, but if they are not, there would never be an intersection to begin with
            // Could also make a list of direct children (removing childrens children), which would cut the search time down when the volumes are nested
            // Check if there is a viable above time_propagated_without_scattering solution
            // Only check the first of the two results in the intersection table, as they are ordered, and the second is of no interest
            if (intersection_time_table.calculated[*check] == 1 && intersection_time_table.intersection_times[*check][0] > time_propagated_without_scattering) {
                // If this child is masked, it's mask status need to be 1 in order to be taken into account
                if (Volumes[*check]->geometry.is_masked_volume == 0) {
                  if (trace_verbal) printf("Found an child of current_volume with an intersection. Skips calculating for current_volume \n");
                  intersection_with_children = 1;
                  break; // No need to check more, if there is just one it is not necessary to calculate intersection with current_volume yet
                } else {
                  if (trace_verbal) printf("Found an child of current_volume with an intersection, but it is masked. Check to see if it can skip calculating for current_volume \n");
                  
                  if (Volumes[*check]->geometry.mask_mode == 2) { // ANY mask mode
                    tree_next_volume = 0;
                    for (mask_start=mask_check=Volumes[*check]->geometry.masked_by_mask_index_list.elements;mask_check-mask_start<Volumes[*check]->geometry.masked_by_mask_index_list.num_elements;mask_check++) {
                       if (mask_status_list.elements[*mask_check] == 1) {
                         intersection_with_children = 1;
                         break;
                       }
                    }
                  } else { // ALL mask mode
                    intersection_with_children = 1;
                    for (mask_start=mask_check=Volumes[*check]->geometry.masked_by_mask_index_list.elements;mask_check-mask_start<Volumes[*check]->geometry.masked_by_mask_index_list.num_elements;mask_check++) {
                      if (mask_status_list.elements[*mask_check] == 0) {
                        intersection_with_children = 0;
                        break;
                      }
                    }
                  }
                  if (trace_verbal) printf("The mask status was 1, can actually skip intersection calculation for current volume \n");
                  if (intersection_with_children == 1) break;
                }
            }
        }
        if (trace_verbal) printf("intersection_with_children = %d \n",intersection_with_children);
        if (intersection_with_children == 0) {
            geometry_output = Volumes[current_volume]->geometry.intersect_function(intersection_time_table.intersection_times[current_volume],number_of_solutions,r_start,v,&Volumes[current_volume]->geometry);
            //geometry_output = Volumes[*check        ]->geometry.intersect_function(intersection_time_table.intersection_times[*check],number_of_solutions,r_start,v,&Volumes[*check        ]->geometry);
            intersection_time_table.calculated[current_volume] = 1;
        }
    }

    // At this point, intersection_time_table is updated with intersection times of all possible intersections.
    if (trace_verbal) print_intersection_table(intersection_time_table);
    // done = 1;
    
    // Find the lowest non negative intersection time on the intersection time lists (but only look in intersection_check_list and active mask_intersection_lists)
    // Declare:
    time_found = 0;
    for (start=check=Volumes[current_volume]->geometry.intersect_check_list.elements;check-start<Volumes[current_volume]->geometry.intersect_check_list.num_elements;check++) {
        for (solution = 0;solution<intersection_time_table.n_elements[*check];solution++) {
            if (time_found) {
                if ((intersection_time = intersection_time_table.intersection_times[*check][solution]) > time_propagated_without_scattering &&  intersection_time < min_intersection_time) {
                    min_intersection_time = intersection_time;min_solution = solution;min_volume = *check;
                }
            } else {
                if ((intersection_time = intersection_time_table.intersection_times[*check][solution]) > time_propagated_without_scattering) {
                    min_intersection_time = intersection_time;min_solution = solution;min_volume = *check;
                    time_found = 1;
                }
            }
        }
    }
    
    
    // Check the mask_intersect_list
    //for (start=check=Volumes[current_volume]->geometry.mask_intersect_list.elements;check-start<Volumes[current_volume]->geometry.mask_intersect_list.num_elements;check++) {
    // cant use the above line as elements may not even be allocated
    for (mask_iterate=0;mask_iterate<Volumes[current_volume]->geometry.mask_intersect_list.num_elements;mask_iterate++) {
      if (current_mask_intersect_list_status.elements[mask_iterate] == 1) {
        for (solution = 0;solution<intersection_time_table.n_elements[Volumes[current_volume]->geometry.mask_intersect_list.elements[mask_iterate]];solution++) {
            if (time_found) {
                if ((intersection_time = intersection_time_table.intersection_times[Volumes[current_volume]->geometry.mask_intersect_list.elements[mask_iterate]][solution]) > time_propagated_without_scattering &&  intersection_time < min_intersection_time) {
                    min_intersection_time = intersection_time;min_solution = solution;min_volume = Volumes[current_volume]->geometry.mask_intersect_list.elements[mask_iterate];
                }
            } else {
                if ((intersection_time = intersection_time_table.intersection_times[Volumes[current_volume]->geometry.mask_intersect_list.elements[mask_iterate]][solution]) > time_propagated_without_scattering) {
                    min_intersection_time = intersection_time;min_solution = solution;min_volume = Volumes[current_volume]->geometry.mask_intersect_list.elements[mask_iterate];
                    time_found = 1;
                }
            }
        }
      }
    }
    
    // And check the current_volume
    for (solution = 0;solution<intersection_time_table.n_elements[current_volume];solution++) {
        if (time_found) {
            if ((intersection_time = intersection_time_table.intersection_times[current_volume][solution]) > time_propagated_without_scattering && intersection_time < min_intersection_time) {
                min_intersection_time = intersection_time;min_solution = solution;min_volume = current_volume;
            }
        } else {
            if ((intersection_time = intersection_time_table.intersection_times[current_volume][solution]) > time_propagated_without_scattering) {
                min_intersection_time = intersection_time;min_solution = solution;min_volume = current_volume;
                time_found = 1;
            }
        }
    }
    
    if (trace_verbal) printf("min_intersection_time = %f \n",min_intersection_time);
    if (trace_verbal) printf("min_solution          = %d \n",min_solution);
    if (trace_verbal) printf("min_volume            = %d \n",min_volume);
    if (trace_verbal) printf("time_found            = %d \n",time_found);


    if (time_found) {
        // Check for errors (debugging phase)
        // Run physics function (include propagating to scattering point, and new velocity)
        time_to_boundery = min_intersection_time - time_propagated_without_scattering;
        scattering_event = 0;
        if (current_volume != 0) {
          if (Volumes[current_volume]->p_physics->number_of_processes == 0) {
            // If there are no processes, the volume could be vacuum or an absorber,
            if (Volumes[current_volume]->p_physics->is_vacuum == 0)
              p *= exp(-Volumes[current_volume]->p_physics->my_a*2200*time_to_boundery);
                    
            //printf("name of material: %s \n",Volumes[current_volume]->name);
            //printf("length to boundery  = %f\n",length_to_boundery);
            //printf("absorption cross section  = %f\n",Volumes[current_volume]->p_physics->my_a);
            //printf("chance to get through this length of absorber: %f \%\n",100*exp(-Volumes[current_volume]->p_physics->my_a*length_to_boundery));
                    
          } else {
            // Since there is a non-zero number of processes in this material, all the scattering cross section for these are calculated
            my_sum = 0; k[0] = V2K*vx; k[1] = V2K*vy; k[2] = V2K*vz; p_my_trace = my_trace; wavevector = coords_set(k[0],k[1],k[2]);
            for (process_start = process = Volumes[current_volume]->p_physics->p_scattering_array;process - process_start < Volumes[current_volume]->p_physics->number_of_processes;process++) {
            
              if (Volumes[current_volume]->p_physics->p_scattering_array[process - process_start].non_isotropic_rot_index != -1) {
                // If the process is not isotropic, the wavevector is transformed into the local coordinate system of the process
                wavevector_rotated = rot_apply(Volumes[current_volume]->geometry.process_rot_matrix_array[Volumes[current_volume]->p_physics->p_scattering_array[process - process_start].non_isotropic_rot_index],wavevector);
                
                coords_get(wavevector_rotated,&k_rotated[0],&k_rotated[1],&k_rotated[2]);

              } else {
                k_rotated[0] = k[0]; k_rotated[1] = k[1]; k_rotated[2] = k[2];
              }
              
              // Call the probability for scattering function assighed to this specific procress (the process pointer is updated in the for loop)
              process->probability_for_scattering_function(p_my_trace,k_rotated,process->data_transfer,&Volumes[current_volume]->geometry.focus_data);
              
              my_sum += *p_my_trace;
              if (trace_verbal) printf("my_trace = %f, my_sum = %f\n",*p_my_trace,my_sum);
              p_my_trace++;
            }
            
            if (trace_verbal) printf("time_propagated_without_scattering = %f.\n",time_propagated_without_scattering);
            if (trace_verbal) printf("v_length                           = %f.\n",v_length);
            
            length_to_boundery = time_to_boundery * v_length;
              
            if (trace_verbal) printf("exp(- length_to_boundery*my_sum) = %f. length_to_boundery = %f. my_sum = %f.\n",exp(-length_to_boundery*my_sum),length_to_boundery,my_sum);
            //printf("exp(- length_to_boundery*my_sum) = %f. length_to_boundery = %f. my_sum = %f.\n",exp(-length_to_boundery*my_sum),length_to_boundery,my_sum);
            
            // Selecting if a scattering takes place, and what scattering process.
            // This section have too many if statements, and unessecary calculations
            // Could make seperate functions for p_interact on/off and interact_fraction on/off,
            //   and set function pointers to these in initialize, thus avoiding many unessecary if statements calculations of x/x.
            
            
            my_sum_plus_abs = my_sum + Volumes[current_volume]->p_physics->my_a*(2200/v_length);
            
            // May need a safeguard for a situation with my_sum = 0
            // NEW SYSTEM
            if (my_sum < 1E-18) {
                // The scattering cross section is basicly zero, no scattering should occur.
                scattering_event = 0;
                p *= exp(-length_to_boundery*my_sum_plus_abs); // Correct for absorption
            } else if (length_to_boundery < safty_distance2) {
                // Experiment 28/11/2016
                // Too close to boundery to safly make another scattering, attenuate
                p *= exp(-length_to_boundery*my_sum_plus_abs); // Attentuate the beam for the small distance
                scattering_event = 0;
            } else {
                // The scattering cross section is above zero
                if (Volumes[current_volume]->geometry.geometry_p_interact != 0) {
                    // a fraction of the beam (geometry_p_interact) is forced to scatter
                    //
                    //real_transmission_probability = exp(-length_to_boundery*my_sum_plus_abs);
                    real_transmission_probability = exp(-length_to_boundery*my_sum_plus_abs);
                    mc_transmission_probability = (1 - Volumes[current_volume]->geometry.geometry_p_interact);
                    if ((scattering_event = (rand01() > mc_transmission_probability))) {
                        // Scattering event happens, this is the correction for the weight
                        p *= (1-real_transmission_probability)/(1-mc_transmission_probability);
                        // Find length to next scattering knowing the ray will scatter.
                        length_to_scattering = safty_distance -log(1 - rand0max((1 - exp(-my_sum_plus_abs*(length_to_boundery-safty_distance2))))) / my_sum_plus_abs;
                    } else
                        // Scattering event does not happen, this is the appropriate correction
                        p *= real_transmission_probability/mc_transmission_probability;
                    
                } else {
                    // probability to scatter is the natural value
                    if(my_sum*length_to_boundery < 1e-6) {
                      if (length_to_boundery > safty_distance2) {
                        if (rand01() < exp(-length_to_boundery*my_sum_plus_abs)) {
                          length_to_scattering = safty_distance + rand0max(length_to_boundery - safty_distance2);
                          p *= length_to_boundery*my_sum*exp(-length_to_scattering*my_sum_plus_abs);
                          scattering_event = 1;
                        } else scattering_event = 0;
                      } else {
                        // The distance is too short to reliably make a scattering event (in comparison to accuraccy of intersect functions)
                        p *= exp(-length_to_boundery*my_sum_plus_abs); // Attentuate the beam for the small distance
                        scattering_event = 0;
                      }
                    } else {
                        length_to_scattering = safty_distance -log(1 - rand01() ) / my_sum_plus_abs;
                        if (length_to_scattering < length_to_boundery - safty_distance) scattering_event = 1;
                        else scattering_event = 0;
                    }
                }
                
                if (scattering_event == 1) {
                  // Adjust weight for absorption
                  p *= my_sum/my_sum_plus_abs;
                  if (my_sum/my_sum_plus_abs > 1) printf("weight factor above 1! Should not happen! \n");
                  // Select process
                  if (Volumes[current_volume]->p_physics->number_of_processes == 1) {
                    // trivial case
                    selected_process = 0;
                  } else {
                    if (Volumes[current_volume]->p_physics->interact_control == 1) {
                      // Interact corrected
                      mc_prop = rand01();culmative_probability=0;total_process_interact=1;
                  
                      // If any of the processes have probability 0, they should be excluded from the selection
                      for (iterate = 0;iterate < Volumes[current_volume]->p_physics->number_of_processes;iterate++) {
                        //printf("my_trace = %E \n",my_trace[iterate]);
                        if (my_trace[iterate] < 1E-18) {
                          //printf("Identified my_trace = 0 \n");
                          // When this happens, the total force probability is corrected and the probability for this particular instance is set to 0
                          total_process_interact -= Volumes[current_volume]->p_physics->p_scattering_array[iterate].process_p_interact;
                          my_trace_fraction_control[iterate] = 0;
                          // In cases where my_trace is not zero, the forced fraction is still used.
                        } else my_trace_fraction_control[iterate] = Volumes[current_volume]->p_physics->p_scattering_array[iterate].process_p_interact;
                      }
                      for (iterate = 0;iterate < Volumes[current_volume]->p_physics->number_of_processes;iterate++) {
                        culmative_probability += my_trace_fraction_control[iterate]/total_process_interact;
                        if (culmative_probability > mc_prop) {
                          selected_process = iterate;
                          p *= (my_trace[iterate]/my_sum)*(total_process_interact/my_trace_fraction_control[iterate]);
                          break;
                        }
                      }
                    
                    } else {
                      // Natural probability distribution
                      mc_prop = rand01();culmative_probability=0;
                      for (iterate = 0;iterate < Volumes[current_volume]->p_physics->number_of_processes;iterate++) {
                        culmative_probability += my_trace[iterate]/my_sum;
                        if (culmative_probability > mc_prop) {
                          selected_process = iterate;
                          break;
                        }
                      }
                    }
                  }
                } // end of select process
            }
            
          }
        } // Done checking for scttering event and in case of scattering selecting process
        
        // Safty measure for testing
        // Test shows it was sufficient to stop problems!
        //if (fabs(length_to_scattering - length_to_boundery) < 1E-15) scattering_event = 0;
     
        if (scattering_event) {
            if (trace_verbal) printf("SCATTERING EVENT \n");
            if (trace_verbal) printf("current_volume            = %d \n",current_volume);
            if (trace_verbal) printf("r = (%f,%f,%f) v = (%f,%f,%f) \n",r[0],r[1],r[2],v[0],v[1],v[2]);
            // printf("did scatter: my_trace = %E = %f \n",my_trace[selected_process],my_trace[selected_process]);
            
            // Select time to propagate before scattering (linear)
            //time_to_scattering = (min_intersection_time - time_propagated_without_scattering) * rand01();
            
            // For debugging the true exponential distribution is tested:
            
            /*
            if(my_trace[selected_process]*length_to_boundery < 1e-6)
                time_to_scattering = rand0max(length_to_boundery);
            else
                //time_to_scattering = -log(1 - rand0max((1 - exp(-my_trace[selected_process]*length_to_boundery)))) / my_trace[selected_process];
                // Updated on the 16/6, as it was noticed the scattering position should depend on the total my, not the one for the selected process
                // It remains to be investigated if it should even contain the absorption my.
                time_to_scattering = -log(1 - rand0max((1 - exp(-my_sum*length_to_boundery)))) / my_sum;
            
            time_to_scattering /= v_length;
            */
            
            time_to_scattering = length_to_scattering/v_length;
            
            // time_to_scattering = (min_intersection_time - time_propagated_without_scattering) * 0.5;
            if (trace_verbal) printf("time to scattering        = %2.20f \n",time_to_scattering);
            //length_to_scattering = time_to_scattering*v_length;
            
            //printf("length to boundery = %f, length to scattering = %f \n",length_to_boundery,length_to_scattering);
            
            // Correcting for using a linear distribution instead of exponential
            //p *= 5*length_to_boundery*my_trace[selected_process]*exp(-length_to_scattering*my_trace[selected_process]);
            
            // Absorption
            //p *= exp(-Volumes[current_volume]->p_physics->my_a*(2200/v_length)*length_to_scattering);
            //p *= exp(-Volumes[current_volume]->p_physics->my_a*2200*time_to_scattering);
            //PROP_DT(time_to_scattering); // May be replace by version without gravity
            
            //if (Volumes[current_volume]->geometry.within_function(ray_position,&Volumes[current_volume]->geometry) == 0)
            //  printf("ERROR, out of volume before propagating to scattering point!\n");
            
            //x += time_to_scattering*vx; y += time_to_scattering*vy; z += time_to_scattering*vz; t += time_to_scattering;
            x += time_to_scattering*vx; y += time_to_scattering*vy; z += time_to_scattering*vz; t += time_to_scattering;
            r_start[0] = x; r_start[1] = y; r_start[2] = z;
            r[0] = x; r[1] = y; r[2] = z;
            ray_position = coords_set(x,y,z);
            ray_velocity = coords_set(vx,vy,vz);
            
            // Safe check that should be unecessary. Used to fine tune how close to the edge of a volume a scattering event is allowed to take place (1E-14 m away currently).
            if (Volumes[current_volume]->geometry.within_function(ray_position,&Volumes[current_volume]->geometry) == 0) {
              printf("\nERROR, propagated out of current volume instead of to a point within!\n");
              printf("length_to_scattering_specified = %2.20f\n             length propagated = %2.20f\n            length_to_boundery = %2.20f \n   current_position = (%lf,%lf,%lf) \n",length_to_scattering,sqrt(time_to_scattering*time_to_scattering*vx*vx+time_to_scattering*time_to_scattering*vy*vy+time_to_scattering*time_to_scattering*vz*vz),length_to_boundery,x,y,z);
              
              volume_index = within_which_volume(ray_position,starting_lists.reduced_start_list,starting_lists.starting_destinations_list,Volumes,&mask_status_list,number_of_volumes,pre_allocated1,pre_allocated2,pre_allocated3);
              
              printf("Debug info: Volumes[current_volume]->name = %s, but now inside volume number %d named %s.\n",Volumes[current_volume]->name,volume_index,Volumes[volume_index]->name);
              printf("Ray absorbed \n");
              ABSORB;
            }
            
            // Save information before scattering event needed in logging section
            p_old = p;
            
        
            // Rotation to local process coordinate system (for non isotropic processes)
            if (Volumes[current_volume]->p_physics->p_scattering_array[selected_process].non_isotropic_rot_index != -1) {
                ray_velocity_rotated = rot_apply(Volumes[current_volume]->geometry.process_rot_matrix_array[Volumes[current_volume]->p_physics->p_scattering_array[selected_process].non_isotropic_rot_index],ray_velocity);
            } else {
                ray_velocity_rotated = ray_velocity;
            }
                
            // test_physics_scattering(double *k_final, double *k_initial, union data_transfer_union data_transfer) {
            //k[0] = V2K*ray_velocity.x; k[1] = V2K*ray_velocity.y; k[2] = V2K*ray_velocity.z;
            coords_get(coords_scalar_mult(ray_velocity_rotated,V2K),&k[0],&k[1],&k[2]);
            
            // I may replace a intial and final k with one instance that serves as both input and output
            if (0 == Volumes[current_volume]->p_physics->p_scattering_array[selected_process].scattering_function(k_new,k,&p,Volumes[current_volume]->p_physics->p_scattering_array[selected_process].data_transfer,&Volumes[current_volume]->geometry.focus_data)) {
              /*
              // PowderN and Single_crystal requires the option of absorbing the neutron, which is weird. If there is a scattering probability, there should be a new direction.
              // It can arise from need to simplify sampling process and end up in cases where weight factor is 0, and the ray should be absorbed in these cases
              printf("ERROR: Union_master: %s.Absorbed ray because scattering function returned 0 (error/absorb)\n",NAME_CURRENT_COMP);
              component_error_msg++;
              if (component_error_msg > 100) {
                printf("To many errors encountered, exiting. \n");
                exit(1);
              }
              */
              ABSORB;
            }
            
            // Update velocity using k
            ray_velocity_rotated = coords_set(K2V*k_new[0],K2V*k_new[1],K2V*k_new[2]);
            
            // Transformation back to main coordinate system (maybe one should only do this when multiple scattering in that volume was over, especially if there is only one non isotropic frame)
            if (Volumes[current_volume]->p_physics->p_scattering_array[selected_process].non_isotropic_rot_index != -1) {
                ray_velocity_final = rot_apply(Volumes[current_volume]->geometry.transpose_process_rot_matrix_array[Volumes[current_volume]->p_physics->p_scattering_array[selected_process].non_isotropic_rot_index],ray_velocity_rotated);
            } else {
               ray_velocity_final = ray_velocity_rotated;
            }
            
            // Write velocity to global variable (temp, only really necessary at final)
            //vx = ray_velocity.x; vy = ray_velocity.y; vz = ray_velocity.z;
            coords_get(ray_velocity_final,&vx,&vy,&vz);
            
            // Write velocity in array format as it is still used by intersect functions (temp, they need to be updated to ray_position / ray_velocity)
            v[0] = vx; v[1] = vy; v[2] = vz;
            v_length = sqrt(vx*vx+vy*vy+vz*vz);
            if (verbal) if (v_length < 1) printf("velocity set to less than 1\n");
            
            
            if (trace_verbal) printf("Running logger system for specific volumes \n");
            // Logging for detector components assosiated with this volume
            for (log_index=0;log_index<Volumes[current_volume]->loggers.num_elements;log_index++) {
              //printf("logging time! Volume specific version. Volume name = %s \n",Volumes[current_volume]->name);
              //printf("  log_index = %d \n",log_index);
              if (Volumes[current_volume]->loggers.p_logger_volume[log_index].p_logger_process[selected_process] != NULL) {
                // Technically the scattering function could edit k, the wavevector before the scattering, even though there would be little point to doing that.
                // Could save a secure copy and pass that instead to be certain that no scattering process accidently tampers with the logging.
                // printf("  the logging function pointer was not NULL \n");
                // PV (number of time scattered in this volume/process combination is not recorded. Need to expand scattering_flag to contain 2D, volume and process
                
                // This function calls a logger function which in turn stores some data among the passed, and possibly performs some basic data analysis
                Volumes[current_volume]->loggers.p_logger_volume[log_index].p_logger_process[selected_process]->function_pointers.active_record_function(&ray_position,k_new,k,p,p_old,t,scattered_flag[current_volume],scattered_flag_VP[current_volume][selected_process],number_of_scattering_events,Volumes[current_volume]->loggers.p_logger_volume[log_index].p_logger_process[selected_process],&loggers_with_data_array);
                // If the logging component have a conditional attatched, the collected data will be written to a temporary place
                // At the end of the rays life, it will be checked if the condition is met
                //  if it is met, the temporary data is transfered to permanent, and temp is cleared.
                //  if it is not met, the temporary data is cleared.
              }
            }
            
            if (trace_verbal) printf("Running logger system for all volumes \n");
            for (log_index=0;log_index<global_all_volume_logger_list.num_elements;log_index++) {
              //printf("logging time! Global version. log_index = %d \n",log_index);
              // As above, but on a global scale, meaning scattering in all volumes are logged
              
              // Problems with VN, PV, as there is no assosiated volume or process. The functions however need to have the same input to make the logger components general.
              // Could be interesting to have a monitor that just globally measurres the second scattering event in any volume (must be two in the same). Weird but not meaningless.
              global_all_volume_logger_list.elements[log_index].logger->function_pointers.active_record_function(&ray_position,k_new,k,p,p_old,t,scattered_flag[current_volume],scattered_flag_VP[current_volume][selected_process],number_of_scattering_events,global_all_volume_logger_list.elements[log_index].logger,&loggers_with_data_array);
            }
            
            
            SCATTER;
            ++number_of_scattering_events;
            ++scattered_flag[current_volume];
            ++scattered_flag_VP[current_volume][selected_process];
            

            // Empty intersection time lists
            clear_intersection_table(&intersection_time_table);
            time_propagated_without_scattering = 0.0;
            if (trace_verbal) printf("SCATTERED SUCSSESFULLY \n");
            if (trace_verbal) printf("r = (%f,%f,%f) v = (%f,%f,%f) \n",x,y,z,vx,vy,vz);
            
            if (trace_verbal && enable_tagging && stop_tagging_ray == 0) printf("Before new process node: current_tagging_node->intensity = %f\n",current_tagging_node->intensity);
            if (trace_verbal && enable_tagging && stop_tagging_ray == 0) printf("Before new process node: current_tagging_node->number_of_rays = %d\n",current_tagging_node->number_of_rays);
            
            if (enable_tagging && stop_tagging_ray == 0)
                current_tagging_node = goto_process_node(current_tagging_node,selected_process,Volumes[current_volume],&stop_tagging_ray,stop_creating_nodes);
            
            if (trace_verbal && enable_tagging && stop_tagging_ray == 0) printf("After new process node: current_tagging_node->intensity = %f\n",current_tagging_node->intensity);
            if (trace_verbal && enable_tagging && stop_tagging_ray == 0) printf("After new process node: current_tagging_node->number_of_rays = %d\n",current_tagging_node->number_of_rays);
            
        } else {
        
            if (trace_verbal) printf("Propagate out of volume %d\n", current_volume);
            if (trace_verbal) printf("r = (%f,%f,%f) v = (%f,%f,%f) \n",x,y,z,vx,vy,vz);
            // Propagate neutron to found minimum time
            // PROP_DT(min_intersection_time - time_propagated_without_scattering);
            //time_to_boundery = min_intersection_time - time_propagated_without_scattering;
            x += time_to_boundery*vx;
            y += time_to_boundery*vy;
            z += time_to_boundery*vz;
            t += time_to_boundery;
            r[0] = x; r[1] = y; r[2] = z;
            ray_position = coords_set(x,y,z);
            ray_velocity = coords_set(vx,vy,vz);
            
            /*
            // Absorption moved to before testing if scattering occurs
            if (current_volume != 0) {
                if (Volumes[current_volume]->p_physics->is_vacuum == 0) {
                    // Absorption is done explicitly when propagating out of a volume, but between all scattering events is done implicitly
                   
                    // Old version
                    //length_to_boundery = (min_intersection_time - time_propagated_without_scattering) * v_length;
                    //p *= exp(-Volumes[current_volume]->p_physics->my_a*(2200/v_length)*length_to_boundery);
                    
                    if (Volumes[current_volume]->p_physics->number_of_processes == 0) {
                      // Optimized version
                      //p *= exp(-Volumes[current_volume]->p_physics->my_a*2200*time_to_boundery);
                    
                    //printf("name of material: %s \n",Volumes[current_volume]->name);
                    //printf("length to boundery  = %f\n",length_to_boundery);
                    //printf("absorption cross section  = %f\n",Volumes[current_volume]->p_physics->my_a);
                    //printf("chance to get through this length of absorber: %f \%\n",100*exp(-Volumes[current_volume]->p_physics->my_a*length_to_boundery));
                    }
                }
            }
            */
            
            time_propagated_without_scattering = min_intersection_time;
            SCATTER; // For debugging purposes
            if (trace_verbal) printf("r = (%f,%f,%f) v = (%f,%f,%f) \n",x,y,z,vx,vy,vz);
            // Remove this entry from the intersection_time_table
            intersection_time_table.intersection_times[min_volume][min_solution] = -1;
            
            // Use destination list for corresponding intersection entry n,i) to find next volume
            if (trace_verbal) printf("PROPAGATION FROM VOLUME %d \n",current_volume);
            if (min_volume == current_volume) {
            
                if (trace_verbal) printf("min_volume == current_volume \n");
                // List approach to finding the next volume.
                // When the ray intersects the current volume, the next volume must be on the destination list of the current volume
                // However, the reduced_destination_list can be investigated first, and depending on the results, the
                //  direct children of the volumes on the reduced destination list are investigated.
                // In the worst case, all direct children are investigated, which is eqvivalent to the entire destination list.
                // There is however a certain overhead in the logic needed to set up this tree, avoid duplicates of direct children, and so on.
                // This method is only faster than just checking the destination list when there are direct children (nested structures),
                //  but in general the tree method scales better with complexity, and is only slightly slower in simple cases.
                
                if (Volumes[current_volume]->geometry.destinations_list.num_elements == 1)
                    tree_next_volume = Volumes[current_volume]->geometry.destinations_list.elements[0];
                else {
                    ray_position = coords_set(x,y,z);
                    ray_velocity = coords_set(vx,vy,vz);
                    tree_next_volume = within_which_volume(ray_position,Volumes[current_volume]->geometry.reduced_destinations_list,Volumes[current_volume]->geometry.destinations_list,Volumes,&mask_status_list,number_of_volumes,pre_allocated1,pre_allocated2,pre_allocated3);
                }

                if (trace_verbal) printf("tree method moves from %d to %d\n",current_volume,tree_next_volume);
                
                if (trace_verbal && enable_tagging && stop_tagging_ray == 0) printf("Before new tree volume node: current_tagging_node->intensity = %f\n",current_tagging_node->intensity);
                if (trace_verbal && enable_tagging && stop_tagging_ray == 0) printf("Before new tree volume node: current_tagging_node->number_of_rays = %d\n",current_tagging_node->number_of_rays);
                
                if (enable_tagging && stop_tagging_ray == 0)
                    current_tagging_node = goto_volume_node(current_tagging_node, current_volume, tree_next_volume, Volumes,&stop_tagging_ray,stop_creating_nodes);
                
                if (trace_verbal && enable_tagging && stop_tagging_ray == 0) printf("After new tree volume node: current_tagging_node->intensity = %f\n",current_tagging_node->intensity);
                if (trace_verbal && enable_tagging && stop_tagging_ray == 0) printf("After new tree volume node: current_tagging_node->number_of_rays = %d\n",current_tagging_node->number_of_rays);
                
                // Set next volume to the solution found in the tree method
                current_volume = tree_next_volume;
                update_current_mask_intersect_status(&current_mask_intersect_list_status, &mask_status_list, Volumes, &current_volume);
                if (trace_verbal) print_1d_int_list(current_mask_intersect_list_status,"Updated current_mask_intersect_list_status");
                
                
                // Debugging phase
                /*
                if (tree_next_volume == 0) {
                    volume_0_found=0;
                    for (start = check = Volumes[current_volume]->geometry.destinations_list.elements;check - start < Volumes[current_volume]->geometry.destinations_list.num_elements;check++) {
                        if (*check == 0) {
                            volume_0_found = 1;
                        }
                    }
                    if (volume_0_found==0) printf("ERROR The within_which_volume function returned 0 for a volume where volume 0 is not on the destination list!\n");
                }
                */
                
            } else {
                if (trace_verbal && enable_tagging && stop_tagging_ray == 0) printf("Before new intersection volume node: current_tagging_node->intensity = %f\n",current_tagging_node->intensity);
                if (trace_verbal && enable_tagging && stop_tagging_ray == 0) printf("Before new intersection volume node: current_tagging_node->number_of_rays = %d\n",current_tagging_node->number_of_rays);
            
                //if (enable_tagging) current_tagging_node = goto_volume_node(current_tagging_node, current_volume, min_volume, Volumes);
                
                
                
                // Mask update: If the min_volume is not a mask, things are simple, current_volume = min_volume.
                //               however, if it is a mask, the mask status will switch.
                //               if the mask status becomes one, the masked volumes inside may be the next volume (unless they are children of the mask)
                //               if the mask status becomes zero (and the current volume is masked by min_volume), the destinations list of the mask is searched
                //               if the mask status becomes zero (and the current volume is NOT masked by min volume), the current volume doesn't change
                
                
                if (Volumes[min_volume]->geometry.is_mask_volume == 0) {
                  if (trace_verbal) printf("Min volume is not a mask, next volume = min volume\n");
                  if (enable_tagging && stop_tagging_ray == 0) current_tagging_node = goto_volume_node(current_tagging_node, current_volume, min_volume, Volumes,&stop_tagging_ray,stop_creating_nodes);
                  current_volume = min_volume;
                  //update_current_mask_intersect_status(&current_mask_intersect_list_status, &mask_status_list, Volumes, &current_volume);
                } else {
                  if (trace_verbal) printf("Current volume is not a mask, complex decision tree\n");
                  if (mask_status_list.elements[Volumes[min_volume]->geometry.mask_index] == 1) {
                    // We are leaving the mask, change the status
                    if (trace_verbal) printf("mask status changed from 1 to 0 as a mask is left\n");
                    mask_status_list.elements[Volumes[min_volume]->geometry.mask_index] = 0;
                    // If the current volume is masked by this mask, run within_which_volume using the masks destination list, otherwise keep the current volume
                    //if (on_int_list(Volumes[min_volume]->geometry.mask_list,current_volume))
                    if (on_int_list(Volumes[current_volume]->geometry.masked_by_list,min_volume) == 1) {
                      if (trace_verbal) printf("The current volume was masked by this mask, and my need updating\n");
                      // In case of ANY mode, need to see if another mask on the masked_by list of the current volume is active, and if so, nothing happens
                      need_to_run_within_which_volume = 1;
                      if (Volumes[current_volume]->geometry.mask_mode == 2) {
                        for (mask_start=mask_check=Volumes[current_volume]->geometry.masked_by_mask_index_list.elements;mask_check-mask_start<Volumes[current_volume]->geometry.masked_by_mask_index_list.num_elements;mask_check++) {
                          if (mask_status_list.elements[*mask_check] == 1) {
                            // Nothing needs to be done, the effective mask status of the current volume is still 1
                            need_to_run_within_which_volume = 0;
                            break;
                          }
                        }
                      }
                      if (need_to_run_within_which_volume == 1) {
                        if (trace_verbal) printf("The current volume was masked by this mask, and does need updating\n");
                        if (Volumes[min_volume]->geometry.destinations_list.num_elements == 1) {
                          if (trace_verbal) printf("Only one element in the destination tree of the mask\n");
                          // If there is only one element on the destinations list (quite common) there is no reason to run within_which_volume
                          // Instead the mask status is calculated here
                          if (Volumes[Volumes[min_volume]->geometry.destinations_list.elements[0]]->geometry.is_masked_volume == 1) {
                            if (trace_verbal) printf("The one element is however masked, so the mask status need to be calculated\n");
                            // figure out the effective mask status of this volume
                            if (Volumes[Volumes[min_volume]->geometry.destinations_list.elements[0]]->geometry.mask_mode == 2) { // ANY mask mode
                              tree_next_volume = 0;
                              for (mask_start=mask_check=Volumes[Volumes[min_volume]->geometry.destinations_list.elements[0]]->geometry.masked_by_mask_index_list.elements;mask_check-mask_start<Volumes[Volumes[min_volume]->geometry.destinations_list.elements[0]]->geometry.masked_by_mask_index_list.num_elements;mask_check++) {
                                if (mask_status_list.elements[*mask_check] == 1) {
                                  tree_next_volume = Volumes[min_volume]->geometry.destinations_list.elements[0];
                                  break;
                                }
                              }
                            } else { // ALL mask mode
                              tree_next_volume = Volumes[min_volume]->geometry.destinations_list.elements[0];
                              for (mask_start=mask_check=Volumes[Volumes[min_volume]->geometry.destinations_list.elements[0]]->geometry.masked_by_mask_index_list.elements;mask_check-mask_start<Volumes[Volumes[min_volume]->geometry.destinations_list.elements[0]]->geometry.masked_by_mask_index_list.num_elements;mask_check++) {
                                if (mask_status_list.elements[*mask_check] == 0) {
                                  tree_next_volume = 0;
                                  break;
                                }
                              }
                            }
                          } else tree_next_volume = Volumes[min_volume]->geometry.destinations_list.elements[0];
                          if (trace_verbal) printf("The method found the next tree volume to be %d\n",tree_next_volume);
                          if (enable_tagging && stop_tagging_ray == 0) current_tagging_node = goto_volume_node(current_tagging_node, current_volume, tree_next_volume, Volumes,&stop_tagging_ray,stop_creating_nodes);
                          current_volume = tree_next_volume;
                          //update_current_mask_intersect_status(&current_mask_intersect_list_status, &mask_status_list, Volumes, &current_volume);
                        } else {
                          if (trace_verbal) printf("Many elements in destinations list, use within_which_volume\n");
                          ray_position = coords_set(x,y,z);
                          ray_velocity = coords_set(vx,vy,vz);
                          tree_next_volume = within_which_volume(ray_position,Volumes[min_volume]->geometry.reduced_destinations_list,Volumes[min_volume]->geometry.destinations_list,Volumes,&mask_status_list,number_of_volumes,pre_allocated1,pre_allocated2,pre_allocated3);
                        // } Bug fixed on 27/11/2016
                          if (enable_tagging && stop_tagging_ray == 0) current_tagging_node = goto_volume_node(current_tagging_node, current_volume, tree_next_volume, Volumes,&stop_tagging_ray,stop_creating_nodes);
                          current_volume = tree_next_volume;
                          if (trace_verbal) printf("Set new new volume to %d\n",tree_next_volume);
                        } // Moved here on 27/11/2016, problem was two calls to current_tagging_node when only one element in destinations list
                        //update_current_mask_intersect_status(&current_mask_intersect_list_status, &mask_status_list, Volumes, &current_volume);
                      } else {
                        if (trace_verbal) printf("Did not need updating, as another mask was covering the volume\n");
                      }
                    }

                  } else {
                    // Here beccause the mask status of the mask that is intersected was 0, and it is thus switched to 1
                    mask_status_list.elements[Volumes[min_volume]->geometry.mask_index] = 1;
                    // When entering a mask, the new highest priority volume may be one of the masked volumes, if not we keep the current volume
                    ray_position = coords_set(x,y,z);
                    ray_velocity = coords_set(vx,vy,vz);
                    //tree_next_volume = within_which_volume(ray_position,Volumes[min_volume]->geometry.mask_list,Volumes[min_volume]->geometry.destinations_list,Volumes,&mask_status_list,number_of_volumes,pre_allocated1,pre_allocated2,pre_allocated3);
                    // Bug found on the 2/9 2016, the destinations_list of a mask does not contain the volumes inside it. Could make an additional list for this.
                    // The temporary fix will be to use the mask list for both reduced destinations list and destinations list.
                    tree_next_volume = within_which_volume(ray_position,Volumes[min_volume]->geometry.mask_list,Volumes[min_volume]->geometry.mask_list,Volumes,&mask_status_list,number_of_volumes,pre_allocated1,pre_allocated2,pre_allocated3);
                    // if within_which_volume returns 0, no result was found (volume 0 can not be masked, so it could not be on the mask list)
                    if (tree_next_volume != 0) {
                      if (Volumes[tree_next_volume]->geometry.priority_value > Volumes[current_volume]->geometry.priority_value) {
                        // In case the current volume has a higher priority, nothing happens, otherwise change current volume
                        if (enable_tagging && stop_tagging_ray == 0) current_tagging_node = goto_volume_node(current_tagging_node, current_volume, tree_next_volume, Volumes,&stop_tagging_ray,stop_creating_nodes);
                        current_volume = tree_next_volume;
                      }
                    }
                    //update_current_mask_intersect_status(&current_mask_intersect_list_status, &mask_status_list, Volumes, &current_volume);
                  }
                }
                
                // Regardless of the outcome of the above code, either the mask status or current volume have changed, and thus a effective mask update is needed.
                update_current_mask_intersect_status(&current_mask_intersect_list_status, &mask_status_list, Volumes, &current_volume);
                if (trace_verbal) print_1d_int_list(mask_status_list,"Updated mask status list");
                if (trace_verbal) print_1d_int_list(current_mask_intersect_list_status,"Updated current_mask_intersect_list_status");
                
                if (trace_verbal && enable_tagging) printf("After new intersection volume node: current_tagging_node->intensity = %f\n",current_tagging_node->intensity);
                if (trace_verbal && enable_tagging) printf("After new intersection volume node: current_tagging_node->number_of_rays = %d\n",current_tagging_node->number_of_rays);
                
            }
            if (Volumes[current_volume]->geometry.is_exit_volume==1) {
                    done = 1; // Exit volumes allow the ray to escape the component
                    ray_sucseeded = 1; // Allows the ray to
                    /*
                    // Moved to after while loop to collect code
                    if (enable_tagging && stop_tagging_ray == 0)
                        add_statistics_to_node(current_tagging_node,&ray_position, &ray_velocity, &p, &tagging_leaf_counter);
                
                    x += vx*1E-9; y += vy*1E-9; z += vz*1E-9; t += 1E-9;
                    */
            }
            if (trace_verbal) printf(" TO VOLUME %d \n",current_volume);
        }
    } else { // Here because a shortest time is not found
        if (current_volume == 0) {
            done = 1;
            ray_sucseeded = 1;
            /*
            // Moved to after while loop to collect code
            if (trace_verbal && enable_tagging && stop_tagging_ray == 0) printf("Before new intersection volume node: current_tagging_nodbe->intensity = %f\n",current_tagging_node->intensity);
            if (trace_verbal && enable_tagging && stop_tagging_ray == 0) printf("Before new intersection volume node: current_tagging_node->number_of_rays = %d\n",current_tagging_node->number_of_rays);
            
            if (enable_tagging && stop_tagging_ray == 0)
                add_statistics_to_node(current_tagging_node,&ray_position, &ray_velocity, &p, &tagging_leaf_counter);
            
            if (trace_verbal && enable_tagging && stop_tagging_ray == 0) printf("After new intersection volume node: current_tagging_node->intensity = %f\n",current_tagging_node->intensity);
            if (trace_verbal && enable_tagging && stop_tagging_ray == 0) printf("After new intersection volume node: current_tagging_node->number_of_rays = %d\n",current_tagging_node->number_of_rays);
            
            x += vx*1E-9; y += vy*1E-9; z += vz*1E-9; t += 1E-9;
            */
            
        } else { // Check for errors (debugging phase)
            if (error_msg == 0) {
              component_error_msg++;
              ray_sucseeded = 0;
              done = 1; // stop the loop
              printf("\n----------------------------------------------------------------------------------------------------\n");
              printf("Union_master %s: Somehow reached a situation with no intersection time found, but still inside volume %d instead of 0\n",NAME_CURRENT_COMP,current_volume);
              for (volume_index = 1; volume_index < number_of_volumes; volume_index++) {
                if (Volumes[volume_index]->geometry.within_function(ray_position,&Volumes[volume_index]->geometry) == 1)
                    printf("The ray is in volume       %d\n",volume_index);
              }
              
              print_1d_int_list(mask_status_list,"mask status list");
              for (iterate=0;iterate<number_of_volumes;iterate++)
                 printf("%d:%d - ",iterate,scattered_flag[iterate]);
              printf("\n");
              printf("r = (%f,%f,%f) v = (%f,%f,%f) \n",x,y,z,vx,vy,vz);
              
              printf("Trace error number (%d/100) \n",component_error_msg);
              //print_intersection_table(intersection_time_table);
              //printf("Cluster is loosing the thread, debug for this behavior\n");
              //printf("global_all_volume_logger_list.num_elements = %d\n",global_all_volume_logger_list.num_elements);
              //printf("global_specific_volumes_logger_list.num_elements = %d\n",global_specific_volumes_logger_list.num_elements);
              
              //exit(1); // Temporary
            }
            error_msg++;
            
            //exit(1); // temp for debug
            
            if (component_error_msg > 100) {
                printf("To many errors encountered, exiting. \n");
                exit(1);
            }
        }
    }
    /*
    */
    if (limit == 0) {done = 1; ray_sucseeded = 0; printf("Reached limit on number of interactions, and discarded the neutron, was in volume %d\n",current_volume);ABSORB;}
    if (trace_verbal) printf("----------- END OF WHILE LOOP --------------------------------------\n");
    //printf("This ray did %d iterations in the while loop\n",1000-limit);
    
  }
  // Could move all add_statistics and similar to this point, but need to filter for failed rays
  if (ray_sucseeded == 1) {
    
    /*
    // Instead of keeping global and specific loggers apart, lets do them in one loop using the loggers_with_data_array
    // Only needed for conditionals
    // Loop over global loggers, may or may not have a conditional, only necessary to do anything here when they do.
    for (log_index=0;log_index<global_all_volume_logger_list.num_elements;log_index++) {
      if (global_all_volume_logger_list.elements[log_index].logger->has_conditional == 1) {
        if (1 == conditional_return = global_all_volume_logger_list.elements[log_index].logger->conditional(scattered_flag,scattered_flag_VP,global_loggers.elements[log_index].conditional_data_union,k_final,current_volume,x,y,z)) {
          global_all_volume_logger_list.elements[log_index].logger->function_pointers->temp_to_perm(global_loggers.elements[log_index].logger_data_union);
          
        }
        if (global_all_volume_logger_list.elements[log_index].logger->conditional_extend_index != -1)
          // The user can set a condtional_extend_index, so that the evaluation of this specific conditional can be taken easily from extend
          conditional_extend_array.elements[global_all_volume_logger_list.elements[log_index].logger->conditional_extend_index] = conditional_return;
          // Do not need to reset these to 0, as they will be manually set to 0 if not fulfilled
      }
    }
    */

    if (trace_verbal) printf("----------- logger loop --------------------------------------\n");
    // Loggers attatched to specific volumes need to be handled with care to avoid looping over all loggers for every ray
    for (log_index=0;log_index<loggers_with_data_array.used_elements;log_index++) {
      // Check all conditionals attatched to the current logger
      this_logger = loggers_with_data_array.logger_pointers[log_index];
      conditional_status = 1;
      for (iterate=0;iterate<loggers_with_data_array.logger_pointers[log_index]->conditional_list.num_elements;iterate++) {
        // Call this particular conditional. If it fails, report the status and break
        //printf("checking conditional! \n");
        if (trace_verbal) printf("Checking conditional number %d for logger named %s \n",iterate,loggers_with_data_array.logger_pointers[log_index]->name);
        if (0 == this_logger->conditional_list.conditional_functions[iterate](
                         this_logger->conditional_list.p_data_unions[iterate],
                         &ray_position,&ray_velocity,&p,&t,&current_volume,
                         &number_of_scattering_events,scattered_flag,scattered_flag_VP)) {
          conditional_status = 0;
          break;
        }
      }
      if (conditional_status == 1) {
        // If a logger does not have a conditional, it will write directly to perm, and not even add it to the loggers_with_data_array, thus we know the temp_to_perm function needs to be called
        // The input for the temp_to_perm function is a pointer to the logger_data_union for the appropriate logger
        loggers_with_data_array.logger_pointers[log_index]->function_pointers.temp_to_perm(&loggers_with_data_array.logger_pointers[log_index]->data_union);
        
        // Luxury feature to be added later
        if (loggers_with_data_array.logger_pointers[log_index]->logger_extend_index != -1) {
          if (trace_verbal) printf("Updating logger_conditional_extend_array[%d] to 1 (max length = %d)\n",loggers_with_data_array.logger_pointers[log_index]->logger_extend_index,max_conditional_extend_index);
          // The user can set a condtional_extend_index, so that the evaluation of this specific conditional can be taken easily from extend
          logger_conditional_extend_array[loggers_with_data_array.logger_pointers[log_index]->logger_extend_index] = 1;
          // Are all reset to 0 for each new ray
          if (trace_verbal) printf("Updated extend index sucessfully\n");
          
          //printf("extend_array[%d] = 1 \n",loggers_with_data_array.logger_pointers[log_index]->logger_extend_index);
        }
        
      } else {
         // Conditional status was 0, clear temp data for this logger
         loggers_with_data_array.logger_pointers[log_index]->function_pointers.clear_temp(&loggers_with_data_array.logger_pointers[log_index]->data_union);
      
      }
    }
    
    if (enable_tagging && stop_tagging_ray == 0) {
      conditional_status = 1;
      for (iterate=0;iterate<tagging_conditional_list->num_elements;iterate++) {
        // Call this particular conditional. If it fails, report the status and break
        // Since a conditional can work for a logger and master_tagging at the same time, it may be evaluated twice
        //printf("checking conditional! \n");
        if (trace_verbal) printf("Checking tagging conditional number %d\n",iterate);
        if (0 == tagging_conditional_list->conditional_functions[iterate](
                         tagging_conditional_list->p_data_unions[iterate],
                         &ray_position,&ray_velocity,&p,&t,&current_volume,
                         &number_of_scattering_events,scattered_flag,scattered_flag_VP)) {
          conditional_status = 0;
          break;
        }
      }
      if (conditional_status == 1) {
        tagging_conditional_extend = 1;
        if (trace_verbal) printf("Before new intersection volume node: current_tagging_nodbe->intensity = %f\n",current_tagging_node->intensity);
        if (trace_verbal) printf("Before new intersection volume node: current_tagging_node->number_of_rays = %d\n",current_tagging_node->number_of_rays);
        
          add_statistics_to_node(current_tagging_node,&ray_position, &ray_velocity, &p, &tagging_leaf_counter);
            
        if (trace_verbal) printf("After new intersection volume node: current_tagging_node->intensity = %f\n",current_tagging_node->intensity);
        if (trace_verbal) printf("After new intersection volume node: current_tagging_node->number_of_rays = %d\n",current_tagging_node->number_of_rays);
      }
    }
    
    // Move the rays a nano meter away from the surface it left, in case activation counter > 1, this will prevent the ray from starting on a volume boundery
    x += vx*1E-9; y += vy*1E-9; z += vz*1E-9; t += 1E-9;
        
  } else {
    ABSORB; // Absorb rays that didn't exit correctly for whatever reason
    // Could error log here
  }
  
  
}
#line 18971 "Union_manual_example.c"
}   /* End of Master=Union_master() SETTING parameter declarations. */
#undef number_of_processes_array
#undef safty_distance2
#undef safty_distance
#undef max_conditional_extend_index
#undef tagging_conditional_extend
#undef logger_conditional_extend_array
#undef free_tagging_conditioanl_list
#undef tagging_conditional_list
#undef conditional_status
#undef this_logger
#undef p_old
#undef loggers_with_data_array
#undef Volume_copies_allocated
#undef geometry_component_index_list
#undef mask_volume_index_list
#undef current_mask_intersect_list_status
#undef mask_status_list
#undef mask_iterate
#undef mask_index_main
#undef need_to_run_within_which_volume
#undef number_of_masked_volumes
#undef number_of_masks
#undef mc_transmission_probability
#undef real_transmission_probability
#undef number_of_scattering_events
#undef tagging_leaf_counter
#undef current_tagging_node
#undef master_tagging_node_list
#undef enable_tagging_check
#undef stop_creating_nodes
#undef stop_tagging_ray
#undef enable_tagging
#undef rotated_position
#undef non_rotated_position
#undef temp_rotation_matrix
#undef master_transposed_rotation_matrix
#undef scattered_flag_VP
#undef scattered_flag
#undef volume_0_found
#undef ray_velocity_final
#undef ray_velocity
#undef ray_position
#undef pre_allocated3
#undef pre_allocated2
#undef pre_allocated1
#undef tree_next_volume
#undef geometry_output
#undef intersection_with_children
#undef start
#undef check
#undef number_of_solutions_static
#undef number_of_solutions
#undef current_volume
#undef done
#undef next_volume_priority
#undef next_volume
#undef a_next_volume_found
#undef time_propagated_without_scattering
#undef scattering_event
#undef selected_process
#undef time_to_boundery
#undef length_to_boundery
#undef length_to_scattering
#undef time_to_scattering
#undef mc_prop
#undef culmative_probability
#undef my_sum_plus_abs
#undef my_sum
#undef v_length
#undef k_new
#undef k
#undef my_trace_fraction_control
#undef p_my_trace
#undef my_trace
#undef process_start
#undef process
#undef min_intersection_time
#undef intersection_time
#undef time_found
#undef min_volume
#undef min_solution
#undef solution
#undef limit
#undef max_number_of_processes
#undef solutions
#undef process_index
#undef volume_index
#undef number_of_volumes
#undef string_output
#undef component_error_msg
#undef error_msg
#undef v
#undef r_start
#undef r
#undef starting_lists
#undef Geometries
#undef Volumes
#undef intersection_time_table
#undef geometry_list_index
#undef previous_master_index
#undef global_master_element
#undef starting_volume_warning
#undef finally_verbal
#undef trace_verbal
#undef list_verbal
#undef verbal
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompMaster:
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

  /* TRACE Component Banana_monitor [15] */
  mccoordschange(mcposrBanana_monitor, mcrotrBanana_monitor,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Banana_monitor (without coords transformations) */
  mcJumpTrace_Banana_monitor:
  SIG_MESSAGE("Banana_monitor (Trace)");
  mcDEBUG_COMP("Banana_monitor")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompBanana_monitor
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
#define mccompcurname  Banana_monitor
#define mccompcurtype  Monitor_nD
#define mccompcurindex 15
#define user1 mccBanana_monitor_user1
#define user2 mccBanana_monitor_user2
#define user3 mccBanana_monitor_user3
#define DEFS mccBanana_monitor_DEFS
#define Vars mccBanana_monitor_Vars
#define detector mccBanana_monitor_detector
#define offdata mccBanana_monitor_offdata
{   /* Declarations of Banana_monitor=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccBanana_monitor_xwidth;
MCNUM yheight = mccBanana_monitor_yheight;
MCNUM zdepth = mccBanana_monitor_zdepth;
MCNUM xmin = mccBanana_monitor_xmin;
MCNUM xmax = mccBanana_monitor_xmax;
MCNUM ymin = mccBanana_monitor_ymin;
MCNUM ymax = mccBanana_monitor_ymax;
MCNUM zmin = mccBanana_monitor_zmin;
MCNUM zmax = mccBanana_monitor_zmax;
MCNUM bins = mccBanana_monitor_bins;
MCNUM min = mccBanana_monitor_min;
MCNUM max = mccBanana_monitor_max;
MCNUM restore_neutron = mccBanana_monitor_restore_neutron;
MCNUM radius = mccBanana_monitor_radius;
char* options = mccBanana_monitor_options;
char* filename = mccBanana_monitor_filename;
char* geometry = mccBanana_monitor_geometry;
char* username1 = mccBanana_monitor_username1;
char* username2 = mccBanana_monitor_username2;
char* username3 = mccBanana_monitor_username3;
int nowritefile = mccBanana_monitor_nowritefile;
#line 310 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
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
#line 19380 "Union_manual_example.c"
}   /* End of Banana_monitor=Monitor_nD() SETTING parameter declarations. */
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
  mcabsorbCompBanana_monitor:
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

  /* User SAVE code for component 'Progress'. */
  SIG_MESSAGE("Progress (Save)");
#define mccompcurname  Progress
#define mccompcurtype  Progress_bar
#define mccompcurindex 7
#define IntermediateCnts mccProgress_IntermediateCnts
#define StartTime mccProgress_StartTime
#define EndTime mccProgress_EndTime
#define CurrentTime mccProgress_CurrentTime
{   /* Declarations of Progress=Progress_bar() SETTING parameters. */
char* profile = mccProgress_profile;
MCNUM percent = mccProgress_percent;
MCNUM flag_save = mccProgress_flag_save;
MCNUM minutes = mccProgress_minutes;
#line 115 "/usr/share/mcstas/2.5/misc/Progress_bar.comp"
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
#line 19496 "Union_manual_example.c"
}   /* End of Progress=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Banana_monitor'. */
  SIG_MESSAGE("Banana_monitor (Save)");
#define mccompcurname  Banana_monitor
#define mccompcurtype  Monitor_nD
#define mccompcurindex 15
#define user1 mccBanana_monitor_user1
#define user2 mccBanana_monitor_user2
#define user3 mccBanana_monitor_user3
#define DEFS mccBanana_monitor_DEFS
#define Vars mccBanana_monitor_Vars
#define detector mccBanana_monitor_detector
#define offdata mccBanana_monitor_offdata
{   /* Declarations of Banana_monitor=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccBanana_monitor_xwidth;
MCNUM yheight = mccBanana_monitor_yheight;
MCNUM zdepth = mccBanana_monitor_zdepth;
MCNUM xmin = mccBanana_monitor_xmin;
MCNUM xmax = mccBanana_monitor_xmax;
MCNUM ymin = mccBanana_monitor_ymin;
MCNUM ymax = mccBanana_monitor_ymax;
MCNUM zmin = mccBanana_monitor_zmin;
MCNUM zmax = mccBanana_monitor_zmax;
MCNUM bins = mccBanana_monitor_bins;
MCNUM min = mccBanana_monitor_min;
MCNUM max = mccBanana_monitor_max;
MCNUM restore_neutron = mccBanana_monitor_restore_neutron;
MCNUM radius = mccBanana_monitor_radius;
char* options = mccBanana_monitor_options;
char* filename = mccBanana_monitor_filename;
char* geometry = mccBanana_monitor_geometry;
char* username1 = mccBanana_monitor_username1;
char* username2 = mccBanana_monitor_username2;
char* username3 = mccBanana_monitor_username3;
int nowritefile = mccBanana_monitor_nowritefile;
#line 480 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
{
  /* save results, but do not free pointers */
  detector = Monitor_nD_Save(&DEFS, &Vars);
}
#line 19545 "Union_manual_example.c"
}   /* End of Banana_monitor=Monitor_nD() SETTING parameter declarations. */
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

  /* User FINALLY code for component 'Al_inc'. */
  SIG_MESSAGE("Al_inc (Finally)");
#define mccompcurname  Al_inc
#define mccompcurtype  Incoherent_process
#define mccompcurindex 1
#define This_process mccAl_inc_This_process
#define Incoherent_storage mccAl_inc_Incoherent_storage
#define effective_my_scattering mccAl_inc_effective_my_scattering
{   /* Declarations of Al_inc=Incoherent_process() SETTING parameters. */
MCNUM sigma = mccAl_inc_sigma;
MCNUM packing_factor = mccAl_inc_packing_factor;
MCNUM unit_cell_volume = mccAl_inc_unit_cell_volume;
MCNUM interact_fraction = mccAl_inc_interact_fraction;
#line 148 "/usr/share/mcstas/2.5/contrib/union/Incoherent_process.comp"
{
// Since the process and it's storage is a static allocation, there is nothing to deallocate

}
#line 19583 "Union_manual_example.c"
}   /* End of Al_inc=Incoherent_process() SETTING parameter declarations. */
#undef effective_my_scattering
#undef Incoherent_storage
#undef This_process
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[1]) fprintf(stderr, "Warning: No neutron could reach Component[1] Al_inc\n");
    if (mcAbsorbProp[1]) fprintf(stderr, "Warning: %g events were removed in Component[1] Al_inc=Incoherent_process()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[1]);
  /* User FINALLY code for component 'Al_powder'. */
  SIG_MESSAGE("Al_powder (Finally)");
#define mccompcurname  Al_powder
#define mccompcurtype  Powder_process
#define mccompcurindex 2
#define format mccAl_powder_format
#define This_process mccAl_powder_This_process
#define Powder_storage mccAl_powder_Powder_storage
#define effective_my_scattering mccAl_powder_effective_my_scattering
#define d_phi mccAl_powder_d_phi
#define line_info mccAl_powder_line_info
#define columns mccAl_powder_columns
{   /* Declarations of Al_powder=Powder_process() SETTING parameters. */
char* reflections = mccAl_powder_reflections;
MCNUM packing_factor = mccAl_powder_packing_factor;
MCNUM Vc = mccAl_powder_Vc;
MCNUM delta_d_d = mccAl_powder_delta_d_d;
MCNUM DW = mccAl_powder_DW;
MCNUM nb_atoms = mccAl_powder_nb_atoms;
MCNUM d_phi = mccAl_powder_d_phi;
MCNUM density = mccAl_powder_density;
MCNUM weight = mccAl_powder_weight;
MCNUM barns = mccAl_powder_barns;
MCNUM Strain = mccAl_powder_Strain;
MCNUM interact_fraction = mccAl_powder_interact_fraction;
#line 765 "/usr/share/mcstas/2.5/contrib/union/Powder_process.comp"
{
  free(line_info.list);
  free(line_info.q_v);
  free(line_info.w_v);
  free(line_info.my_s_v2);
}
#line 19626 "Union_manual_example.c"
}   /* End of Al_powder=Powder_process() SETTING parameter declarations. */
#undef columns
#undef line_info
#undef d_phi
#undef effective_my_scattering
#undef Powder_storage
#undef This_process
#undef format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[2]) fprintf(stderr, "Warning: No neutron could reach Component[2] Al_powder\n");
    if (mcAbsorbProp[2]) fprintf(stderr, "Warning: %g events were removed in Component[2] Al_powder=Powder_process()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[2]);
  /* User FINALLY code for component 'Al'. */
  SIG_MESSAGE("Al (Finally)");
#define mccompcurname  Al
#define mccompcurtype  Union_make_material
#define mccompcurindex 3
#define loop_index mccAl_loop_index
#define this_material mccAl_this_material
#define accepted_processes mccAl_accepted_processes
#define global_material_element mccAl_global_material_element
{   /* Declarations of Al=Union_make_material() SETTING parameters. */
char* process_string = mccAl_process_string;
MCNUM my_absorption = mccAl_my_absorption;
MCNUM absorber = mccAl_absorber;
#line 259 "/usr/share/mcstas/2.5/contrib/union/Union_make_material.comp"
{
// The elements of the scattering array used static allocation and is thus deallocated automatically
if (this_material.number_of_processes > 0) free(this_material.p_scattering_array);
if (accepted_processes.num_elements > 0) free(accepted_processes.elements);

// Checking if any Union volumes are defined after the master component
#ifdef MASTER_DETECTOR
  #ifdef ANY_GEOMETRY_DETECTOR_DECLARE
    #ifndef MASTER_DETECTOR_WARNING
      for (loop_index=0;loop_index<global_geometry_list.num_elements;loop_index++) {
        if (global_geometry_list.elements[loop_index].component_index > global_master_list.elements[global_master_list.num_elements-1].component_index) {
          printf("WARNING: No Union_master component defined after Union volume named %s, this components did not affect the simulation in any way.\n",global_geometry_list.elements[loop_index].name);
        }
      }
      // Decided to have this as a warning without exiting the simulation
      // In order to only show this warning once, the MASTER_DETECTOR_WARNING is defined
      #define MASTER_DETECTOR_WARNING dummy
    #endif
  #endif
#endif

// Checking if the user remembered to put in a Union_master
#ifndef MASTER_DETECTOR
  #ifdef ANY_GEOMETRY_DETECTOR_DECLARE
    #ifndef MASTER_DETECTOR_WARNING
      printf("\nWARNING: No Union_master component used, these components did not affect the simulation in any way:\n");
      for (loop_index=0;loop_index<global_geometry_list.num_elements;loop_index++)
        printf("  %s\n",global_geometry_list.elements[loop_index].name);
      printf("\n");
      // Decided to have this as a warning without exiting the simulation
      // In order to only show this warning once, the MASTER_DETECTOR_WARNING is defined
      #define MASTER_DETECTOR_WARNING dummy
    #endif
  #endif
#endif


}
#line 19693 "Union_manual_example.c"
}   /* End of Al=Union_make_material() SETTING parameter declarations. */
#undef global_material_element
#undef accepted_processes
#undef this_material
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[3]) fprintf(stderr, "Warning: No neutron could reach Component[3] Al\n");
    if (mcAbsorbProp[3]) fprintf(stderr, "Warning: %g events were removed in Component[3] Al=Union_make_material()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[3]);
  /* User FINALLY code for component 'Fe_inc'. */
  SIG_MESSAGE("Fe_inc (Finally)");
#define mccompcurname  Fe_inc
#define mccompcurtype  Incoherent_process
#define mccompcurindex 4
#define This_process mccFe_inc_This_process
#define Incoherent_storage mccFe_inc_Incoherent_storage
#define effective_my_scattering mccFe_inc_effective_my_scattering
{   /* Declarations of Fe_inc=Incoherent_process() SETTING parameters. */
MCNUM sigma = mccFe_inc_sigma;
MCNUM packing_factor = mccFe_inc_packing_factor;
MCNUM unit_cell_volume = mccFe_inc_unit_cell_volume;
MCNUM interact_fraction = mccFe_inc_interact_fraction;
#line 148 "/usr/share/mcstas/2.5/contrib/union/Incoherent_process.comp"
{
// Since the process and it's storage is a static allocation, there is nothing to deallocate

}
#line 19723 "Union_manual_example.c"
}   /* End of Fe_inc=Incoherent_process() SETTING parameter declarations. */
#undef effective_my_scattering
#undef Incoherent_storage
#undef This_process
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[4]) fprintf(stderr, "Warning: No neutron could reach Component[4] Fe_inc\n");
    if (mcAbsorbProp[4]) fprintf(stderr, "Warning: %g events were removed in Component[4] Fe_inc=Incoherent_process()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[4]);
  /* User FINALLY code for component 'Fe_powder'. */
  SIG_MESSAGE("Fe_powder (Finally)");
#define mccompcurname  Fe_powder
#define mccompcurtype  Powder_process
#define mccompcurindex 5
#define format mccFe_powder_format
#define This_process mccFe_powder_This_process
#define Powder_storage mccFe_powder_Powder_storage
#define effective_my_scattering mccFe_powder_effective_my_scattering
#define d_phi mccFe_powder_d_phi
#define line_info mccFe_powder_line_info
#define columns mccFe_powder_columns
{   /* Declarations of Fe_powder=Powder_process() SETTING parameters. */
char* reflections = mccFe_powder_reflections;
MCNUM packing_factor = mccFe_powder_packing_factor;
MCNUM Vc = mccFe_powder_Vc;
MCNUM delta_d_d = mccFe_powder_delta_d_d;
MCNUM DW = mccFe_powder_DW;
MCNUM nb_atoms = mccFe_powder_nb_atoms;
MCNUM d_phi = mccFe_powder_d_phi;
MCNUM density = mccFe_powder_density;
MCNUM weight = mccFe_powder_weight;
MCNUM barns = mccFe_powder_barns;
MCNUM Strain = mccFe_powder_Strain;
MCNUM interact_fraction = mccFe_powder_interact_fraction;
#line 765 "/usr/share/mcstas/2.5/contrib/union/Powder_process.comp"
{
  free(line_info.list);
  free(line_info.q_v);
  free(line_info.w_v);
  free(line_info.my_s_v2);
}
#line 19766 "Union_manual_example.c"
}   /* End of Fe_powder=Powder_process() SETTING parameter declarations. */
#undef columns
#undef line_info
#undef d_phi
#undef effective_my_scattering
#undef Powder_storage
#undef This_process
#undef format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[5]) fprintf(stderr, "Warning: No neutron could reach Component[5] Fe_powder\n");
    if (mcAbsorbProp[5]) fprintf(stderr, "Warning: %g events were removed in Component[5] Fe_powder=Powder_process()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[5]);
  /* User FINALLY code for component 'Fe'. */
  SIG_MESSAGE("Fe (Finally)");
#define mccompcurname  Fe
#define mccompcurtype  Union_make_material
#define mccompcurindex 6
#define loop_index mccFe_loop_index
#define this_material mccFe_this_material
#define accepted_processes mccFe_accepted_processes
#define global_material_element mccFe_global_material_element
{   /* Declarations of Fe=Union_make_material() SETTING parameters. */
char* process_string = mccFe_process_string;
MCNUM my_absorption = mccFe_my_absorption;
MCNUM absorber = mccFe_absorber;
#line 259 "/usr/share/mcstas/2.5/contrib/union/Union_make_material.comp"
{
// The elements of the scattering array used static allocation and is thus deallocated automatically
if (this_material.number_of_processes > 0) free(this_material.p_scattering_array);
if (accepted_processes.num_elements > 0) free(accepted_processes.elements);

// Checking if any Union volumes are defined after the master component
#ifdef MASTER_DETECTOR
  #ifdef ANY_GEOMETRY_DETECTOR_DECLARE
    #ifndef MASTER_DETECTOR_WARNING
      for (loop_index=0;loop_index<global_geometry_list.num_elements;loop_index++) {
        if (global_geometry_list.elements[loop_index].component_index > global_master_list.elements[global_master_list.num_elements-1].component_index) {
          printf("WARNING: No Union_master component defined after Union volume named %s, this components did not affect the simulation in any way.\n",global_geometry_list.elements[loop_index].name);
        }
      }
      // Decided to have this as a warning without exiting the simulation
      // In order to only show this warning once, the MASTER_DETECTOR_WARNING is defined
      #define MASTER_DETECTOR_WARNING dummy
    #endif
  #endif
#endif

// Checking if the user remembered to put in a Union_master
#ifndef MASTER_DETECTOR
  #ifdef ANY_GEOMETRY_DETECTOR_DECLARE
    #ifndef MASTER_DETECTOR_WARNING
      printf("\nWARNING: No Union_master component used, these components did not affect the simulation in any way:\n");
      for (loop_index=0;loop_index<global_geometry_list.num_elements;loop_index++)
        printf("  %s\n",global_geometry_list.elements[loop_index].name);
      printf("\n");
      // Decided to have this as a warning without exiting the simulation
      // In order to only show this warning once, the MASTER_DETECTOR_WARNING is defined
      #define MASTER_DETECTOR_WARNING dummy
    #endif
  #endif
#endif


}
#line 19833 "Union_manual_example.c"
}   /* End of Fe=Union_make_material() SETTING parameter declarations. */
#undef global_material_element
#undef accepted_processes
#undef this_material
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[6]) fprintf(stderr, "Warning: No neutron could reach Component[6] Fe\n");
    if (mcAbsorbProp[6]) fprintf(stderr, "Warning: %g events were removed in Component[6] Fe=Union_make_material()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[6]);
  /* User FINALLY code for component 'Progress'. */
  SIG_MESSAGE("Progress (Finally)");
#define mccompcurname  Progress
#define mccompcurtype  Progress_bar
#define mccompcurindex 7
#define IntermediateCnts mccProgress_IntermediateCnts
#define StartTime mccProgress_StartTime
#define EndTime mccProgress_EndTime
#define CurrentTime mccProgress_CurrentTime
{   /* Declarations of Progress=Progress_bar() SETTING parameters. */
char* profile = mccProgress_profile;
MCNUM percent = mccProgress_percent;
MCNUM flag_save = mccProgress_flag_save;
MCNUM minutes = mccProgress_minutes;
#line 133 "/usr/share/mcstas/2.5/misc/Progress_bar.comp"
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
#line 19872 "Union_manual_example.c"
}   /* End of Progress=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[7]) fprintf(stderr, "Warning: No neutron could reach Component[7] Progress\n");
    if (mcAbsorbProp[7]) fprintf(stderr, "Warning: %g events were removed in Component[7] Progress=Progress_bar()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[7]);
    if (!mcNCounter[8]) fprintf(stderr, "Warning: No neutron could reach Component[8] Source\n");
    if (mcAbsorbProp[8]) fprintf(stderr, "Warning: %g events were removed in Component[8] Source=Source_simple()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[8]);
  /* User FINALLY code for component 'Guide'. */
  SIG_MESSAGE("Guide (Finally)");
#define mccompcurname  Guide
#define mccompcurtype  Guide_gravity
#define mccompcurindex 9
#define GVars mccGuide_GVars
#define pTable mccGuide_pTable
{   /* Declarations of Guide=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccGuide_w1;
MCNUM h1 = mccGuide_h1;
MCNUM w2 = mccGuide_w2;
MCNUM h2 = mccGuide_h2;
MCNUM l = mccGuide_l;
MCNUM R0 = mccGuide_R0;
MCNUM Qc = mccGuide_Qc;
MCNUM alpha = mccGuide_alpha;
MCNUM m = mccGuide_m;
MCNUM W = mccGuide_W;
MCNUM nslit = mccGuide_nslit;
MCNUM d = mccGuide_d;
MCNUM mleft = mccGuide_mleft;
MCNUM mright = mccGuide_mright;
MCNUM mtop = mccGuide_mtop;
MCNUM mbottom = mccGuide_mbottom;
MCNUM nhslit = mccGuide_nhslit;
MCNUM G = mccGuide_G;
MCNUM aleft = mccGuide_aleft;
MCNUM aright = mccGuide_aright;
MCNUM atop = mccGuide_atop;
MCNUM abottom = mccGuide_abottom;
MCNUM wavy = mccGuide_wavy;
MCNUM wavy_z = mccGuide_wavy_z;
MCNUM wavy_tb = mccGuide_wavy_tb;
MCNUM wavy_lr = mccGuide_wavy_lr;
MCNUM chamfers = mccGuide_chamfers;
MCNUM chamfers_z = mccGuide_chamfers_z;
MCNUM chamfers_lr = mccGuide_chamfers_lr;
MCNUM chamfers_tb = mccGuide_chamfers_tb;
MCNUM nelements = mccGuide_nelements;
MCNUM nu = mccGuide_nu;
MCNUM phase = mccGuide_phase;
char* reflect = mccGuide_reflect;
#line 562 "/usr/share/mcstas/2.5/optics/Guide_gravity.comp"
{
if (GVars.warnings > 100) {
  fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
  fprintf(stderr,"%s: warning: This message has been repeated %g times\n", GVars.compcurname, GVars.warnings);
}
}
#line 19935 "Union_manual_example.c"
}   /* End of Guide=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[9]) fprintf(stderr, "Warning: No neutron could reach Component[9] Guide\n");
    if (mcAbsorbProp[9]) fprintf(stderr, "Warning: %g events were removed in Component[9] Guide=Guide_gravity()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[9]);
    if (!mcNCounter[10]) fprintf(stderr, "Warning: No neutron could reach Component[10] sample_position\n");
    if (mcAbsorbProp[10]) fprintf(stderr, "Warning: %g events were removed in Component[10] sample_position=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[10]);
    if (!mcNCounter[11]) fprintf(stderr, "Warning: No neutron could reach Component[11] Hexagonal_container_1\n");
    if (mcAbsorbProp[11]) fprintf(stderr, "Warning: %g events were removed in Component[11] Hexagonal_container_1=Union_box()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[11]);
    if (!mcNCounter[12]) fprintf(stderr, "Warning: No neutron could reach Component[12] Hexagonal_container_2\n");
    if (mcAbsorbProp[12]) fprintf(stderr, "Warning: %g events were removed in Component[12] Hexagonal_container_2=Union_box()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[12]);
    if (!mcNCounter[13]) fprintf(stderr, "Warning: No neutron could reach Component[13] Sample\n");
    if (mcAbsorbProp[13]) fprintf(stderr, "Warning: %g events were removed in Component[13] Sample=Union_cylinder()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[13]);
  /* User FINALLY code for component 'Master'. */
  SIG_MESSAGE("Master (Finally)");
#define mccompcurname  Master
#define mccompcurtype  Union_master
#define mccompcurindex 14
#define verbal mccMaster_verbal
#define list_verbal mccMaster_list_verbal
#define trace_verbal mccMaster_trace_verbal
#define finally_verbal mccMaster_finally_verbal
#define starting_volume_warning mccMaster_starting_volume_warning
#define global_master_element mccMaster_global_master_element
#define previous_master_index mccMaster_previous_master_index
#define geometry_list_index mccMaster_geometry_list_index
#define intersection_time_table mccMaster_intersection_time_table
#define Volumes mccMaster_Volumes
#define Geometries mccMaster_Geometries
#define starting_lists mccMaster_starting_lists
#define r mccMaster_r
#define r_start mccMaster_r_start
#define v mccMaster_v
#define error_msg mccMaster_error_msg
#define component_error_msg mccMaster_component_error_msg
#define string_output mccMaster_string_output
#define number_of_volumes mccMaster_number_of_volumes
#define volume_index mccMaster_volume_index
#define process_index mccMaster_process_index
#define solutions mccMaster_solutions
#define max_number_of_processes mccMaster_max_number_of_processes
#define limit mccMaster_limit
#define solution mccMaster_solution
#define min_solution mccMaster_min_solution
#define min_volume mccMaster_min_volume
#define time_found mccMaster_time_found
#define intersection_time mccMaster_intersection_time
#define min_intersection_time mccMaster_min_intersection_time
#define process mccMaster_process
#define process_start mccMaster_process_start
#define my_trace mccMaster_my_trace
#define p_my_trace mccMaster_p_my_trace
#define my_trace_fraction_control mccMaster_my_trace_fraction_control
#define k mccMaster_k
#define k_new mccMaster_k_new
#define v_length mccMaster_v_length
#define my_sum mccMaster_my_sum
#define my_sum_plus_abs mccMaster_my_sum_plus_abs
#define culmative_probability mccMaster_culmative_probability
#define mc_prop mccMaster_mc_prop
#define time_to_scattering mccMaster_time_to_scattering
#define length_to_scattering mccMaster_length_to_scattering
#define length_to_boundery mccMaster_length_to_boundery
#define time_to_boundery mccMaster_time_to_boundery
#define selected_process mccMaster_selected_process
#define scattering_event mccMaster_scattering_event
#define time_propagated_without_scattering mccMaster_time_propagated_without_scattering
#define a_next_volume_found mccMaster_a_next_volume_found
#define next_volume mccMaster_next_volume
#define next_volume_priority mccMaster_next_volume_priority
#define done mccMaster_done
#define current_volume mccMaster_current_volume
#define number_of_solutions mccMaster_number_of_solutions
#define number_of_solutions_static mccMaster_number_of_solutions_static
#define check mccMaster_check
#define start mccMaster_start
#define intersection_with_children mccMaster_intersection_with_children
#define geometry_output mccMaster_geometry_output
#define tree_next_volume mccMaster_tree_next_volume
#define pre_allocated1 mccMaster_pre_allocated1
#define pre_allocated2 mccMaster_pre_allocated2
#define pre_allocated3 mccMaster_pre_allocated3
#define ray_position mccMaster_ray_position
#define ray_velocity mccMaster_ray_velocity
#define ray_velocity_final mccMaster_ray_velocity_final
#define volume_0_found mccMaster_volume_0_found
#define scattered_flag mccMaster_scattered_flag
#define scattered_flag_VP mccMaster_scattered_flag_VP
#define master_transposed_rotation_matrix mccMaster_master_transposed_rotation_matrix
#define temp_rotation_matrix mccMaster_temp_rotation_matrix
#define non_rotated_position mccMaster_non_rotated_position
#define rotated_position mccMaster_rotated_position
#define enable_tagging mccMaster_enable_tagging
#define stop_tagging_ray mccMaster_stop_tagging_ray
#define stop_creating_nodes mccMaster_stop_creating_nodes
#define enable_tagging_check mccMaster_enable_tagging_check
#define master_tagging_node_list mccMaster_master_tagging_node_list
#define current_tagging_node mccMaster_current_tagging_node
#define tagging_leaf_counter mccMaster_tagging_leaf_counter
#define number_of_scattering_events mccMaster_number_of_scattering_events
#define real_transmission_probability mccMaster_real_transmission_probability
#define mc_transmission_probability mccMaster_mc_transmission_probability
#define number_of_masks mccMaster_number_of_masks
#define number_of_masked_volumes mccMaster_number_of_masked_volumes
#define need_to_run_within_which_volume mccMaster_need_to_run_within_which_volume
#define mask_index_main mccMaster_mask_index_main
#define mask_iterate mccMaster_mask_iterate
#define mask_status_list mccMaster_mask_status_list
#define current_mask_intersect_list_status mccMaster_current_mask_intersect_list_status
#define mask_volume_index_list mccMaster_mask_volume_index_list
#define geometry_component_index_list mccMaster_geometry_component_index_list
#define Volume_copies_allocated mccMaster_Volume_copies_allocated
#define loggers_with_data_array mccMaster_loggers_with_data_array
#define p_old mccMaster_p_old
#define this_logger mccMaster_this_logger
#define conditional_status mccMaster_conditional_status
#define tagging_conditional_list mccMaster_tagging_conditional_list
#define free_tagging_conditioanl_list mccMaster_free_tagging_conditioanl_list
#define logger_conditional_extend_array mccMaster_logger_conditional_extend_array
#define tagging_conditional_extend mccMaster_tagging_conditional_extend
#define max_conditional_extend_index mccMaster_max_conditional_extend_index
#define safty_distance mccMaster_safty_distance
#define safty_distance2 mccMaster_safty_distance2
#define number_of_processes_array mccMaster_number_of_processes_array
{   /* Declarations of Master=Union_master() SETTING parameters. */
MCNUM allow_inside_start = mccMaster_allow_inside_start;
MCNUM history_limit = mccMaster_history_limit;
#line 1861 "/usr/share/mcstas/2.5/contrib/union/Union_master.comp"
{
// write out histories from tagging system if enabled
if (enable_tagging) {
    if (finally_verbal) printf("Writing tagging tree to disk \n");
    if (finally_verbal) printf("Number of leafs = %d \n",tagging_leaf_counter);
    // While writing the tagging tree to disk, all the leafs are deallocated
    write_tagging_tree(&master_tagging_node_list, Volumes, tagging_leaf_counter, number_of_volumes);
}
if (master_tagging_node_list.num_elements > 0) free(master_tagging_node_list.elements);




if (finally_verbal) printf("Freeing variables which are always allocated \n");
// free allocated arrays specific to this master union component
free(scattered_flag);
free(my_trace);
free(pre_allocated1);
free(pre_allocated2);
free(pre_allocated3);
free(number_of_processes_array);

if (finally_verbal) printf("Freeing intersection_time_table \n");
for (iterate = 1;iterate < intersection_time_table.num_volumes;iterate++){
    free(intersection_time_table.intersection_times[iterate]);
}

free(intersection_time_table.n_elements);
free(intersection_time_table.calculated);
free(intersection_time_table.intersection_times);

if (free_tagging_conditioanl_list == 1) free(tagging_conditional_list);

/*
if (tagging_conditional_list->num_elements > 0) {
  free(tagging_conditional_list.conditional_functions);
  free(tagging_conditional_list.p_data_unions);
}
*/

if (finally_verbal) printf("Freeing lists for individual volumes \n");
for (volume_index=0;volume_index<number_of_volumes;volume_index++) {
  
  if (Volumes[volume_index]->geometry.intersect_check_list.num_elements > 0) free(Volumes[volume_index]->geometry.intersect_check_list.elements);
  if (Volumes[volume_index]->geometry.destinations_list.num_elements > 0) free(Volumes[volume_index]->geometry.destinations_list.elements);
  if (Volumes[volume_index]->geometry.reduced_destinations_list.num_elements > 0) free(Volumes[volume_index]->geometry.reduced_destinations_list.elements);
  if (Volumes[volume_index]->geometry.children.num_elements > 0) free(Volumes[volume_index]->geometry.children.elements);
  if (Volumes[volume_index]->geometry.direct_children.num_elements > 0) free(Volumes[volume_index]->geometry.direct_children.elements);
  if (Volumes[volume_index]->geometry.masked_by_list.num_elements > 0) free(Volumes[volume_index]->geometry.masked_by_list.elements);
  if (Volumes[volume_index]->geometry.masked_by_mask_index_list.num_elements > 0) free(Volumes[volume_index]->geometry.masked_by_mask_index_list.elements);
  if (Volumes[volume_index]->geometry.mask_list.num_elements > 0) free(Volumes[volume_index]->geometry.mask_list.elements);
  if (Volumes[volume_index]->geometry.mask_intersect_list.num_elements > 0) free(Volumes[volume_index]->geometry.mask_intersect_list.elements);
  if (enable_tagging)
    if (Volumes[volume_index]->geometry.next_volume_list.num_elements > 0) free(Volumes[volume_index]->geometry.next_volume_list.elements);
  // Add dealocation of logging
  

  if (volume_index > 0) { // Volume 0 does not have physical properties allocated
    free(scattered_flag_VP[volume_index]);
    if (Volumes[volume_index]->geometry.process_rot_allocated == 1) {
          free(Volumes[volume_index]->geometry.process_rot_matrix_array);
          free(Volumes[volume_index]->geometry.transpose_process_rot_matrix_array);
          Volumes[volume_index]->geometry.process_rot_allocated = 0;
    }
    if (on_int_list(Volume_copies_allocated,volume_index))
      // This is a local copy of a volume, deallocate that local copy (all the allocated memory attachted to it was just deallocated, so this should not leave any leaks)
      free(Volumes[volume_index]);
    else
      // Only free p_physics for vacuum volumes for the original at the end (there is a p_physics allocated for each vacuum volume)
      if (Volumes[volume_index]->p_physics->is_vacuum == 1 ) free(Volumes[volume_index]->p_physics);
  }
  
  if (Volumes[volume_index]->loggers.num_elements >0) {
    for (iterate=0;iterate<Volumes[volume_index]->loggers.num_elements;iterate++) {
      free(Volumes[volume_index]->loggers.p_logger_volume[iterate].p_logger_process);
    }
    free(Volumes[volume_index]->loggers.p_logger_volume);
  }
  
}

free(scattered_flag_VP);

if (finally_verbal) printf("Freeing starting lists \n");
if (starting_lists.allowed_starting_volume_logic_list.num_elements > 0) free(starting_lists.allowed_starting_volume_logic_list.elements);
if (starting_lists.reduced_start_list.num_elements > 0) free(starting_lists.reduced_start_list.elements);
if (starting_lists.start_logic_list.num_elements > 0) free(starting_lists.start_logic_list.elements);

if (finally_verbal) printf("Freeing mask lists \n");
if (mask_status_list.num_elements>0) free(mask_status_list.elements);
if (current_mask_intersect_list_status.num_elements>0) free(current_mask_intersect_list_status.elements);
if (mask_volume_index_list.num_elements>0) free(mask_volume_index_list.elements);

if (finally_verbal) printf("Freeing component index list \n");
if (geometry_component_index_list.num_elements>0) free(geometry_component_index_list.elements);


if (finally_verbal) printf("Freeing Volumes \n");
free(Volumes);

// Free global allocated arrays if this is the last master union component in the instrument file

if (global_master_list.elements[global_master_list.num_elements-1].component_index == INDEX_CURRENT_COMP) {
    if (finally_verbal) printf("Freeing global arrays because this is the last Union master component\n");
    
    // Freeing lists allocated in Union_initialization
    //#ifdef PROCESS_DETECTOR
        if (finally_verbal) printf("Freeing global process list \n");
        if (global_process_list.num_elements > 0) free(global_process_list.elements);
    //#endif

    //#ifdef MATERIAL_DETECTOR
        if (finally_verbal) printf("Freeing global material list \n");
        if (global_material_list.num_elements > 0) free(global_material_list.elements);
    //#endif

    //#ifdef ANY_GEOMETRY_DETECTOR_DECLARE
        if (finally_verbal) printf("Freeing global geometry list \n");
        if (global_geometry_list.num_elements > 0) free(global_geometry_list.elements);
    //#endif
    
    //#ifdef MASTER_DETECTOR
        if (finally_verbal) printf("Freeing global master list \n");
        if (global_master_list.num_elements > 0) free(global_master_list.elements);
    //#endif
    
    
    //#ifdef UNION_LOGGER_DECLARE
    if (finally_verbal) printf("Freeing global logger lists \n");
    for (iterate=0;iterate<global_all_volume_logger_list.num_elements;iterate++) {
      if (global_all_volume_logger_list.elements[iterate].logger->conditional_list.num_elements > 0) {
        free(global_all_volume_logger_list.elements[iterate].logger->conditional_list.conditional_functions);
        free(global_all_volume_logger_list.elements[iterate].logger->conditional_list.p_data_unions);
      }
    }
    if (global_all_volume_logger_list.num_elements > 0) free(global_all_volume_logger_list.elements);
    
    
    for (iterate=0;iterate<global_specific_volumes_logger_list.num_elements;iterate++) {
      if (global_specific_volumes_logger_list.elements[iterate].logger->conditional_list.num_elements > 0) {
        free(global_specific_volumes_logger_list.elements[iterate].logger->conditional_list.conditional_functions);
        free(global_specific_volumes_logger_list.elements[iterate].logger->conditional_list.p_data_unions);
      }
    }
    if (global_specific_volumes_logger_list.num_elements > 0) free(global_specific_volumes_logger_list.elements);
    //#endif
    
    for (iterate=0;iterate<global_tagging_conditional_list.num_elements;iterate++) {
      if (global_tagging_conditional_list.elements[iterate].conditional_list.num_elements > 0) {
        free(global_tagging_conditional_list.elements[iterate].conditional_list.conditional_functions);
        free(global_tagging_conditional_list.elements[iterate].conditional_list.p_data_unions);
      }
    }
    if (global_tagging_conditional_list.num_elements>0) free(global_tagging_conditional_list.elements);
    
    /*
    if (finally_verbal) printf("Freeing global tagging conditional list \n");
    if (global_tagging_conditional_list.num_elements > 0) free(global_tagging_conditional_list.elements);
    */
}

}
#line 20230 "Union_manual_example.c"
}   /* End of Master=Union_master() SETTING parameter declarations. */
#undef number_of_processes_array
#undef safty_distance2
#undef safty_distance
#undef max_conditional_extend_index
#undef tagging_conditional_extend
#undef logger_conditional_extend_array
#undef free_tagging_conditioanl_list
#undef tagging_conditional_list
#undef conditional_status
#undef this_logger
#undef p_old
#undef loggers_with_data_array
#undef Volume_copies_allocated
#undef geometry_component_index_list
#undef mask_volume_index_list
#undef current_mask_intersect_list_status
#undef mask_status_list
#undef mask_iterate
#undef mask_index_main
#undef need_to_run_within_which_volume
#undef number_of_masked_volumes
#undef number_of_masks
#undef mc_transmission_probability
#undef real_transmission_probability
#undef number_of_scattering_events
#undef tagging_leaf_counter
#undef current_tagging_node
#undef master_tagging_node_list
#undef enable_tagging_check
#undef stop_creating_nodes
#undef stop_tagging_ray
#undef enable_tagging
#undef rotated_position
#undef non_rotated_position
#undef temp_rotation_matrix
#undef master_transposed_rotation_matrix
#undef scattered_flag_VP
#undef scattered_flag
#undef volume_0_found
#undef ray_velocity_final
#undef ray_velocity
#undef ray_position
#undef pre_allocated3
#undef pre_allocated2
#undef pre_allocated1
#undef tree_next_volume
#undef geometry_output
#undef intersection_with_children
#undef start
#undef check
#undef number_of_solutions_static
#undef number_of_solutions
#undef current_volume
#undef done
#undef next_volume_priority
#undef next_volume
#undef a_next_volume_found
#undef time_propagated_without_scattering
#undef scattering_event
#undef selected_process
#undef time_to_boundery
#undef length_to_boundery
#undef length_to_scattering
#undef time_to_scattering
#undef mc_prop
#undef culmative_probability
#undef my_sum_plus_abs
#undef my_sum
#undef v_length
#undef k_new
#undef k
#undef my_trace_fraction_control
#undef p_my_trace
#undef my_trace
#undef process_start
#undef process
#undef min_intersection_time
#undef intersection_time
#undef time_found
#undef min_volume
#undef min_solution
#undef solution
#undef limit
#undef max_number_of_processes
#undef solutions
#undef process_index
#undef volume_index
#undef number_of_volumes
#undef string_output
#undef component_error_msg
#undef error_msg
#undef v
#undef r_start
#undef r
#undef starting_lists
#undef Geometries
#undef Volumes
#undef intersection_time_table
#undef geometry_list_index
#undef previous_master_index
#undef global_master_element
#undef starting_volume_warning
#undef finally_verbal
#undef trace_verbal
#undef list_verbal
#undef verbal
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[14]) fprintf(stderr, "Warning: No neutron could reach Component[14] Master\n");
    if (mcAbsorbProp[14]) fprintf(stderr, "Warning: %g events were removed in Component[14] Master=Union_master()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[14]);
  /* User FINALLY code for component 'Banana_monitor'. */
  SIG_MESSAGE("Banana_monitor (Finally)");
#define mccompcurname  Banana_monitor
#define mccompcurtype  Monitor_nD
#define mccompcurindex 15
#define user1 mccBanana_monitor_user1
#define user2 mccBanana_monitor_user2
#define user3 mccBanana_monitor_user3
#define DEFS mccBanana_monitor_DEFS
#define Vars mccBanana_monitor_Vars
#define detector mccBanana_monitor_detector
#define offdata mccBanana_monitor_offdata
{   /* Declarations of Banana_monitor=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccBanana_monitor_xwidth;
MCNUM yheight = mccBanana_monitor_yheight;
MCNUM zdepth = mccBanana_monitor_zdepth;
MCNUM xmin = mccBanana_monitor_xmin;
MCNUM xmax = mccBanana_monitor_xmax;
MCNUM ymin = mccBanana_monitor_ymin;
MCNUM ymax = mccBanana_monitor_ymax;
MCNUM zmin = mccBanana_monitor_zmin;
MCNUM zmax = mccBanana_monitor_zmax;
MCNUM bins = mccBanana_monitor_bins;
MCNUM min = mccBanana_monitor_min;
MCNUM max = mccBanana_monitor_max;
MCNUM restore_neutron = mccBanana_monitor_restore_neutron;
MCNUM radius = mccBanana_monitor_radius;
char* options = mccBanana_monitor_options;
char* filename = mccBanana_monitor_filename;
char* geometry = mccBanana_monitor_geometry;
char* username1 = mccBanana_monitor_username1;
char* username2 = mccBanana_monitor_username2;
char* username3 = mccBanana_monitor_username3;
int nowritefile = mccBanana_monitor_nowritefile;
#line 486 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
{
  /* free pointers */
  if (!nowritefile) {
    Monitor_nD_Finally(&DEFS, &Vars);
  }
}
#line 20385 "Union_manual_example.c"
}   /* End of Banana_monitor=Monitor_nD() SETTING parameter declarations. */
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

    if (!mcNCounter[15]) fprintf(stderr, "Warning: No neutron could reach Component[15] Banana_monitor\n");
    if (mcAbsorbProp[15]) fprintf(stderr, "Warning: %g events were removed in Component[15] Banana_monitor=Monitor_nD()\n"
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

  /* MCDISPLAY code for component 'Progress'. */
  SIG_MESSAGE("Progress (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Progress");
#define mccompcurname  Progress
#define mccompcurtype  Progress_bar
#define mccompcurindex 7
#define IntermediateCnts mccProgress_IntermediateCnts
#define StartTime mccProgress_StartTime
#define EndTime mccProgress_EndTime
#define CurrentTime mccProgress_CurrentTime
{   /* Declarations of Progress=Progress_bar() SETTING parameters. */
char* profile = mccProgress_profile;
MCNUM percent = mccProgress_percent;
MCNUM flag_save = mccProgress_flag_save;
MCNUM minutes = mccProgress_minutes;
#line 147 "/usr/share/mcstas/2.5/misc/Progress_bar.comp"
{
  
}
#line 20434 "Union_manual_example.c"
}   /* End of Progress=Progress_bar() SETTING parameter declarations. */
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
#define mccompcurindex 8
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
#line 171 "/usr/share/mcstas/2.5/sources/Source_simple.comp"
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
#line 20483 "Union_manual_example.c"
}   /* End of Source=Source_simple() SETTING parameter declarations. */
#undef srcArea
#undef square
#undef pmul
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Guide'. */
  SIG_MESSAGE("Guide (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Guide");
#define mccompcurname  Guide
#define mccompcurtype  Guide_gravity
#define mccompcurindex 9
#define GVars mccGuide_GVars
#define pTable mccGuide_pTable
{   /* Declarations of Guide=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccGuide_w1;
MCNUM h1 = mccGuide_h1;
MCNUM w2 = mccGuide_w2;
MCNUM h2 = mccGuide_h2;
MCNUM l = mccGuide_l;
MCNUM R0 = mccGuide_R0;
MCNUM Qc = mccGuide_Qc;
MCNUM alpha = mccGuide_alpha;
MCNUM m = mccGuide_m;
MCNUM W = mccGuide_W;
MCNUM nslit = mccGuide_nslit;
MCNUM d = mccGuide_d;
MCNUM mleft = mccGuide_mleft;
MCNUM mright = mccGuide_mright;
MCNUM mtop = mccGuide_mtop;
MCNUM mbottom = mccGuide_mbottom;
MCNUM nhslit = mccGuide_nhslit;
MCNUM G = mccGuide_G;
MCNUM aleft = mccGuide_aleft;
MCNUM aright = mccGuide_aright;
MCNUM atop = mccGuide_atop;
MCNUM abottom = mccGuide_abottom;
MCNUM wavy = mccGuide_wavy;
MCNUM wavy_z = mccGuide_wavy_z;
MCNUM wavy_tb = mccGuide_wavy_tb;
MCNUM wavy_lr = mccGuide_wavy_lr;
MCNUM chamfers = mccGuide_chamfers;
MCNUM chamfers_z = mccGuide_chamfers_z;
MCNUM chamfers_lr = mccGuide_chamfers_lr;
MCNUM chamfers_tb = mccGuide_chamfers_tb;
MCNUM nelements = mccGuide_nelements;
MCNUM nu = mccGuide_nu;
MCNUM phase = mccGuide_phase;
char* reflect = mccGuide_reflect;
#line 571 "/usr/share/mcstas/2.5/optics/Guide_gravity.comp"
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
#line 20598 "Union_manual_example.c"
}   /* End of Guide=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'sample_position'. */
  SIG_MESSAGE("sample_position (McDisplay)");
  printf("MCDISPLAY: component %s\n", "sample_position");
#define mccompcurname  sample_position
#define mccompcurtype  Arm
#define mccompcurindex 10
#line 40 "/usr/share/mcstas/2.5/optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 20620 "Union_manual_example.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Master'. */
  SIG_MESSAGE("Master (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Master");
#define mccompcurname  Master
#define mccompcurtype  Union_master
#define mccompcurindex 14
#define verbal mccMaster_verbal
#define list_verbal mccMaster_list_verbal
#define trace_verbal mccMaster_trace_verbal
#define finally_verbal mccMaster_finally_verbal
#define starting_volume_warning mccMaster_starting_volume_warning
#define global_master_element mccMaster_global_master_element
#define previous_master_index mccMaster_previous_master_index
#define geometry_list_index mccMaster_geometry_list_index
#define intersection_time_table mccMaster_intersection_time_table
#define Volumes mccMaster_Volumes
#define Geometries mccMaster_Geometries
#define starting_lists mccMaster_starting_lists
#define r mccMaster_r
#define r_start mccMaster_r_start
#define v mccMaster_v
#define error_msg mccMaster_error_msg
#define component_error_msg mccMaster_component_error_msg
#define string_output mccMaster_string_output
#define number_of_volumes mccMaster_number_of_volumes
#define volume_index mccMaster_volume_index
#define process_index mccMaster_process_index
#define solutions mccMaster_solutions
#define max_number_of_processes mccMaster_max_number_of_processes
#define limit mccMaster_limit
#define solution mccMaster_solution
#define min_solution mccMaster_min_solution
#define min_volume mccMaster_min_volume
#define time_found mccMaster_time_found
#define intersection_time mccMaster_intersection_time
#define min_intersection_time mccMaster_min_intersection_time
#define process mccMaster_process
#define process_start mccMaster_process_start
#define my_trace mccMaster_my_trace
#define p_my_trace mccMaster_p_my_trace
#define my_trace_fraction_control mccMaster_my_trace_fraction_control
#define k mccMaster_k
#define k_new mccMaster_k_new
#define v_length mccMaster_v_length
#define my_sum mccMaster_my_sum
#define my_sum_plus_abs mccMaster_my_sum_plus_abs
#define culmative_probability mccMaster_culmative_probability
#define mc_prop mccMaster_mc_prop
#define time_to_scattering mccMaster_time_to_scattering
#define length_to_scattering mccMaster_length_to_scattering
#define length_to_boundery mccMaster_length_to_boundery
#define time_to_boundery mccMaster_time_to_boundery
#define selected_process mccMaster_selected_process
#define scattering_event mccMaster_scattering_event
#define time_propagated_without_scattering mccMaster_time_propagated_without_scattering
#define a_next_volume_found mccMaster_a_next_volume_found
#define next_volume mccMaster_next_volume
#define next_volume_priority mccMaster_next_volume_priority
#define done mccMaster_done
#define current_volume mccMaster_current_volume
#define number_of_solutions mccMaster_number_of_solutions
#define number_of_solutions_static mccMaster_number_of_solutions_static
#define check mccMaster_check
#define start mccMaster_start
#define intersection_with_children mccMaster_intersection_with_children
#define geometry_output mccMaster_geometry_output
#define tree_next_volume mccMaster_tree_next_volume
#define pre_allocated1 mccMaster_pre_allocated1
#define pre_allocated2 mccMaster_pre_allocated2
#define pre_allocated3 mccMaster_pre_allocated3
#define ray_position mccMaster_ray_position
#define ray_velocity mccMaster_ray_velocity
#define ray_velocity_final mccMaster_ray_velocity_final
#define volume_0_found mccMaster_volume_0_found
#define scattered_flag mccMaster_scattered_flag
#define scattered_flag_VP mccMaster_scattered_flag_VP
#define master_transposed_rotation_matrix mccMaster_master_transposed_rotation_matrix
#define temp_rotation_matrix mccMaster_temp_rotation_matrix
#define non_rotated_position mccMaster_non_rotated_position
#define rotated_position mccMaster_rotated_position
#define enable_tagging mccMaster_enable_tagging
#define stop_tagging_ray mccMaster_stop_tagging_ray
#define stop_creating_nodes mccMaster_stop_creating_nodes
#define enable_tagging_check mccMaster_enable_tagging_check
#define master_tagging_node_list mccMaster_master_tagging_node_list
#define current_tagging_node mccMaster_current_tagging_node
#define tagging_leaf_counter mccMaster_tagging_leaf_counter
#define number_of_scattering_events mccMaster_number_of_scattering_events
#define real_transmission_probability mccMaster_real_transmission_probability
#define mc_transmission_probability mccMaster_mc_transmission_probability
#define number_of_masks mccMaster_number_of_masks
#define number_of_masked_volumes mccMaster_number_of_masked_volumes
#define need_to_run_within_which_volume mccMaster_need_to_run_within_which_volume
#define mask_index_main mccMaster_mask_index_main
#define mask_iterate mccMaster_mask_iterate
#define mask_status_list mccMaster_mask_status_list
#define current_mask_intersect_list_status mccMaster_current_mask_intersect_list_status
#define mask_volume_index_list mccMaster_mask_volume_index_list
#define geometry_component_index_list mccMaster_geometry_component_index_list
#define Volume_copies_allocated mccMaster_Volume_copies_allocated
#define loggers_with_data_array mccMaster_loggers_with_data_array
#define p_old mccMaster_p_old
#define this_logger mccMaster_this_logger
#define conditional_status mccMaster_conditional_status
#define tagging_conditional_list mccMaster_tagging_conditional_list
#define free_tagging_conditioanl_list mccMaster_free_tagging_conditioanl_list
#define logger_conditional_extend_array mccMaster_logger_conditional_extend_array
#define tagging_conditional_extend mccMaster_tagging_conditional_extend
#define max_conditional_extend_index mccMaster_max_conditional_extend_index
#define safty_distance mccMaster_safty_distance
#define safty_distance2 mccMaster_safty_distance2
#define number_of_processes_array mccMaster_number_of_processes_array
{   /* Declarations of Master=Union_master() SETTING parameters. */
MCNUM allow_inside_start = mccMaster_allow_inside_start;
MCNUM history_limit = mccMaster_history_limit;
#line 2026 "/usr/share/mcstas/2.5/contrib/union/Union_master.comp"
{
  // mcdisplay is handled in the component files for each geometry and called here. The line function is only available in this section, and not through functions,
  //   so all the lines to be drawn for each volume are collected in a structure that is then drawn here.
  
  struct lines_to_draw lines_to_draw_master;
  for (volume_index = 1; volume_index < number_of_volumes; volume_index++) {
        if (Volumes[volume_index]->geometry.visualization_on == 1) {
            lines_to_draw_master.number_of_lines = 0;
            
            Volumes[volume_index]->geometry.mcdisplay_function(&lines_to_draw_master,volume_index,Geometries,number_of_volumes);
            
            for (iterate = 0;iterate<lines_to_draw_master.number_of_lines;iterate++) {
               if (lines_to_draw_master.lines[iterate].number_of_dashes == 1) {
                 line(lines_to_draw_master.lines[iterate].point1.x,lines_to_draw_master.lines[iterate].point1.y,lines_to_draw_master.lines[iterate].point1.z,lines_to_draw_master.lines[iterate].point2.x,lines_to_draw_master.lines[iterate].point2.y,lines_to_draw_master.lines[iterate].point2.z);
               }
               else {
                 dashed_line(lines_to_draw_master.lines[iterate].point1.x,lines_to_draw_master.lines[iterate].point1.y,lines_to_draw_master.lines[iterate].point1.z,lines_to_draw_master.lines[iterate].point2.x,lines_to_draw_master.lines[iterate].point2.y,lines_to_draw_master.lines[iterate].point2.z,lines_to_draw_master.lines[iterate].number_of_dashes);
               }
            }
            
            if (lines_to_draw_master.number_of_lines>0) free(lines_to_draw_master.lines);
        }
   }

}
#line 20766 "Union_manual_example.c"
}   /* End of Master=Union_master() SETTING parameter declarations. */
#undef number_of_processes_array
#undef safty_distance2
#undef safty_distance
#undef max_conditional_extend_index
#undef tagging_conditional_extend
#undef logger_conditional_extend_array
#undef free_tagging_conditioanl_list
#undef tagging_conditional_list
#undef conditional_status
#undef this_logger
#undef p_old
#undef loggers_with_data_array
#undef Volume_copies_allocated
#undef geometry_component_index_list
#undef mask_volume_index_list
#undef current_mask_intersect_list_status
#undef mask_status_list
#undef mask_iterate
#undef mask_index_main
#undef need_to_run_within_which_volume
#undef number_of_masked_volumes
#undef number_of_masks
#undef mc_transmission_probability
#undef real_transmission_probability
#undef number_of_scattering_events
#undef tagging_leaf_counter
#undef current_tagging_node
#undef master_tagging_node_list
#undef enable_tagging_check
#undef stop_creating_nodes
#undef stop_tagging_ray
#undef enable_tagging
#undef rotated_position
#undef non_rotated_position
#undef temp_rotation_matrix
#undef master_transposed_rotation_matrix
#undef scattered_flag_VP
#undef scattered_flag
#undef volume_0_found
#undef ray_velocity_final
#undef ray_velocity
#undef ray_position
#undef pre_allocated3
#undef pre_allocated2
#undef pre_allocated1
#undef tree_next_volume
#undef geometry_output
#undef intersection_with_children
#undef start
#undef check
#undef number_of_solutions_static
#undef number_of_solutions
#undef current_volume
#undef done
#undef next_volume_priority
#undef next_volume
#undef a_next_volume_found
#undef time_propagated_without_scattering
#undef scattering_event
#undef selected_process
#undef time_to_boundery
#undef length_to_boundery
#undef length_to_scattering
#undef time_to_scattering
#undef mc_prop
#undef culmative_probability
#undef my_sum_plus_abs
#undef my_sum
#undef v_length
#undef k_new
#undef k
#undef my_trace_fraction_control
#undef p_my_trace
#undef my_trace
#undef process_start
#undef process
#undef min_intersection_time
#undef intersection_time
#undef time_found
#undef min_volume
#undef min_solution
#undef solution
#undef limit
#undef max_number_of_processes
#undef solutions
#undef process_index
#undef volume_index
#undef number_of_volumes
#undef string_output
#undef component_error_msg
#undef error_msg
#undef v
#undef r_start
#undef r
#undef starting_lists
#undef Geometries
#undef Volumes
#undef intersection_time_table
#undef geometry_list_index
#undef previous_master_index
#undef global_master_element
#undef starting_volume_warning
#undef finally_verbal
#undef trace_verbal
#undef list_verbal
#undef verbal
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Banana_monitor'. */
  SIG_MESSAGE("Banana_monitor (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Banana_monitor");
#define mccompcurname  Banana_monitor
#define mccompcurtype  Monitor_nD
#define mccompcurindex 15
#define user1 mccBanana_monitor_user1
#define user2 mccBanana_monitor_user2
#define user3 mccBanana_monitor_user3
#define DEFS mccBanana_monitor_DEFS
#define Vars mccBanana_monitor_Vars
#define detector mccBanana_monitor_detector
#define offdata mccBanana_monitor_offdata
{   /* Declarations of Banana_monitor=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccBanana_monitor_xwidth;
MCNUM yheight = mccBanana_monitor_yheight;
MCNUM zdepth = mccBanana_monitor_zdepth;
MCNUM xmin = mccBanana_monitor_xmin;
MCNUM xmax = mccBanana_monitor_xmax;
MCNUM ymin = mccBanana_monitor_ymin;
MCNUM ymax = mccBanana_monitor_ymax;
MCNUM zmin = mccBanana_monitor_zmin;
MCNUM zmax = mccBanana_monitor_zmax;
MCNUM bins = mccBanana_monitor_bins;
MCNUM min = mccBanana_monitor_min;
MCNUM max = mccBanana_monitor_max;
MCNUM restore_neutron = mccBanana_monitor_restore_neutron;
MCNUM radius = mccBanana_monitor_radius;
char* options = mccBanana_monitor_options;
char* filename = mccBanana_monitor_filename;
char* geometry = mccBanana_monitor_geometry;
char* username1 = mccBanana_monitor_username1;
char* username2 = mccBanana_monitor_username2;
char* username3 = mccBanana_monitor_username3;
int nowritefile = mccBanana_monitor_nowritefile;
#line 494 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
{
  if (geometry && strlen(geometry) && strcmp(geometry,"0") && strcmp(geometry, "NULL"))
  {
    off_display(offdata);
  } else {
    Monitor_nD_McDisplay(&DEFS, &Vars);
  }
}
#line 20922 "Union_manual_example.c"
}   /* End of Banana_monitor=Monitor_nD() SETTING parameter declarations. */
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
/* end of generated C code Union_manual_example.c */
