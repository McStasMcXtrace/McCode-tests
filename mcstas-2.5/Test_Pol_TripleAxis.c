/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: Test_Pol_TripleAxis.instr (Test_Pol_TripleAxis)
 * Date:       Fri Oct 11 20:25:20 2019
 * File:       Test_Pol_TripleAxis.c
 * Compile:    cc -o Test_Pol_TripleAxis.out Test_Pol_TripleAxis.c 
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

#line 695 "Test_Pol_TripleAxis.c"

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

#line 928 "Test_Pol_TripleAxis.c"

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

#line 4947 "Test_Pol_TripleAxis.c"

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

#line 5307 "Test_Pol_TripleAxis.c"
#ifdef MC_TRACE_ENABLED
int mctraceenabled = 1;
#else
int mctraceenabled = 0;
#endif
#define MCSTAS "/usr/share/mcstas/2.5/"
int mcdefaultmain = 1;
char mcinstrument_name[] = "Test_Pol_TripleAxis";
char mcinstrument_source[] = "Test_Pol_TripleAxis.instr";
char *mcinstrument_exe=NULL; /* will be set to argv[0] in main */
int main(int argc, char *argv[]){return mccode_main(argc, argv);}
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);

/* Shared user declarations for all components 'Monochromator_flat'. */
#line 71 "/usr/share/mcstas/2.5/optics/Monochromator_flat.comp"
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
#line 5347 "Test_Pol_TripleAxis.c"

/* Shared user declarations for all components 'Monochromator_pol'. */
#line 95 "/usr/share/mcstas/2.5/optics/Monochromator_pol.comp"
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
    MAGNET_OFF; \
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
    MAGNET_ON;
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
    MAGNET_ON;
  }else{
    MAGNET_OFF;
  }
  return (void*) t;
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

#line 8163 "Test_Pol_TripleAxis.c"

/* Shared user declarations for all components 'V_sample'. */
#line 100 "/usr/share/mcstas/2.5/obsolete/V_sample.comp"
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
#line 8180 "Test_Pol_TripleAxis.c"

/* Shared user declarations for all components 'Pol_constBfield'. */
#line 49 "/usr/share/mcstas/2.5/optics/Pol_constBfield.comp"
double IntersectWall(double pos, double vel, double wallpos) {
    /* Function to calculate where the neutron hit the wall */

    if(vel==0)
      return -1;

    if(vel>0)
      return (wallpos-pos)/vel;
    else
      return (-wallpos-pos)/vel;
  }
#line 8195 "Test_Pol_TripleAxis.c"

/* Instrument parameters. */
int mcipOPTION;
MCNUM mcipLAMBDA;
MCNUM mcipMOZ;

#define mcNUMIPAR 3
int mcnumipar = 3;
struct mcinputtable_struct mcinputtable[mcNUMIPAR+1] = {
  "OPTION", &mcipOPTION, instr_type_int, "", 
  "LAMBDA", &mcipLAMBDA, instr_type_double, "2.0", 
  "MOZ", &mcipMOZ, instr_type_double, "40", 
  NULL, NULL, instr_type_double, ""
};

/* User declarations from instrument definition. */
#define mccompcurname  Test_Pol_TripleAxis
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaTest_Pol_TripleAxis coords_set(0,0,0)
#define OPTION mcipOPTION
#define LAMBDA mcipLAMBDA
#define MOZ mcipMOZ
#line 42 "Test_Pol_TripleAxis.instr"
  double DM  = 3.3539; /* Monochromator d-spacing in Angs */
                       /* PG002 Orders : 1st 3.355 2e 1.6775, 3e 1.1183 */

  /* to compute */
  double A1,A2;
  double A3,A4;
  double A5,A6;
  double mono_q, Ei;

  /* This variable helps us switch on and of the spin flipper*/
  double filterProb;
#line 8231 "Test_Pol_TripleAxis.c"
#undef MOZ
#undef LAMBDA
#undef OPTION
#undef mcposaTest_Pol_TripleAxis
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

/* Setting parameters for component 'Mono' [5]. */
MCNUM mccMono_zmin;
MCNUM mccMono_zmax;
MCNUM mccMono_ymin;
MCNUM mccMono_ymax;
MCNUM mccMono_zwidth;
MCNUM mccMono_yheight;
MCNUM mccMono_mosaich;
MCNUM mccMono_mosaicv;
MCNUM mccMono_r0;
MCNUM mccMono_Q;
MCNUM mccMono_DM;

/* Setting parameters for component 'Mono_pol' [6]. */
MCNUM mccMono_pol_zwidth;
MCNUM mccMono_pol_yheight;
MCNUM mccMono_pol_mosaic;
MCNUM mccMono_pol_dspread;
MCNUM mccMono_pol_Q;
MCNUM mccMono_pol_DM;
MCNUM mccMono_pol_pThreshold;
MCNUM mccMono_pol_Rup;
MCNUM mccMono_pol_Rdown;
int mccMono_pol_debug;

/* Definition parameters for component 'pollambdaMonitorMono' [7]. */
#define mccpollambdaMonitorMono_xwidth 0.1
#define mccpollambdaMonitorMono_yheight 0.1
#define mccpollambdaMonitorMono_nL 40
#define mccpollambdaMonitorMono_npol 11
#define mccpollambdaMonitorMono_restore_neutron 0
#define mccpollambdaMonitorMono_nowritefile 0
/* Setting parameters for component 'pollambdaMonitorMono' [7]. */
char mccpollambdaMonitorMono_filename[16384];
MCNUM mccpollambdaMonitorMono_mx;
MCNUM mccpollambdaMonitorMono_my;
MCNUM mccpollambdaMonitorMono_mz;
MCNUM mccpollambdaMonitorMono_Lmin;
MCNUM mccpollambdaMonitorMono_Lmax;

/* Setting parameters for component 'Sample' [10]. */
MCNUM mccSample_radius;
MCNUM mccSample_thickness;
MCNUM mccSample_zdepth;
MCNUM mccSample_Vc;
MCNUM mccSample_sigma_abs;
MCNUM mccSample_sigma_inc;
MCNUM mccSample_radius_i;
MCNUM mccSample_radius_o;
MCNUM mccSample_h;
MCNUM mccSample_focus_r;
MCNUM mccSample_pack;
MCNUM mccSample_frac;
MCNUM mccSample_f_QE;
MCNUM mccSample_gamma;
MCNUM mccSample_target_x;
MCNUM mccSample_target_y;
MCNUM mccSample_target_z;
MCNUM mccSample_focus_xw;
MCNUM mccSample_focus_yh;
MCNUM mccSample_focus_aw;
MCNUM mccSample_focus_ah;
MCNUM mccSample_xwidth;
MCNUM mccSample_yheight;
MCNUM mccSample_zthick;
MCNUM mccSample_rad_sphere;
MCNUM mccSample_sig_a;
MCNUM mccSample_sig_i;
MCNUM mccSample_V0;
int mccSample_target_index;
MCNUM mccSample_multiples;

/* Definition parameters for component 'pollambdaMonitorSample' [11]. */
#define mccpollambdaMonitorSample_xwidth 0.10
#define mccpollambdaMonitorSample_yheight 0.10
#define mccpollambdaMonitorSample_nL 40
#define mccpollambdaMonitorSample_npol 11
#define mccpollambdaMonitorSample_restore_neutron 0
#define mccpollambdaMonitorSample_nowritefile 0
/* Setting parameters for component 'pollambdaMonitorSample' [11]. */
char mccpollambdaMonitorSample_filename[16384];
MCNUM mccpollambdaMonitorSample_mx;
MCNUM mccpollambdaMonitorSample_my;
MCNUM mccpollambdaMonitorSample_mz;
MCNUM mccpollambdaMonitorSample_Lmin;
MCNUM mccpollambdaMonitorSample_Lmax;

/* Definition parameters for component 'msf' [12]. */
#define mccmsf_xwidth 0.10
#define mccmsf_yheight 0.10
#define mccmsf_length 0.2
/* Setting parameters for component 'msf' [12]. */
MCNUM mccmsf_B;
MCNUM mccmsf_fliplambda;
MCNUM mccmsf_flipangle;

/* Definition parameters for component 'pollambdaMonitorMSF' [13]. */
#define mccpollambdaMonitorMSF_xwidth 0.1
#define mccpollambdaMonitorMSF_yheight 0.1
#define mccpollambdaMonitorMSF_nL 40
#define mccpollambdaMonitorMSF_npol 11
#define mccpollambdaMonitorMSF_restore_neutron 0
#define mccpollambdaMonitorMSF_nowritefile 0
/* Setting parameters for component 'pollambdaMonitorMSF' [13]. */
char mccpollambdaMonitorMSF_filename[16384];
MCNUM mccpollambdaMonitorMSF_mx;
MCNUM mccpollambdaMonitorMSF_my;
MCNUM mccpollambdaMonitorMSF_mz;
MCNUM mccpollambdaMonitorMSF_Lmin;
MCNUM mccpollambdaMonitorMSF_Lmax;

/* Setting parameters for component 'Ana' [16]. */
MCNUM mccAna_zmin;
MCNUM mccAna_zmax;
MCNUM mccAna_ymin;
MCNUM mccAna_ymax;
MCNUM mccAna_zwidth;
MCNUM mccAna_yheight;
MCNUM mccAna_mosaich;
MCNUM mccAna_mosaicv;
MCNUM mccAna_r0;
MCNUM mccAna_Q;
MCNUM mccAna_DM;

/* Setting parameters for component 'Ana_pol' [17]. */
MCNUM mccAna_pol_zwidth;
MCNUM mccAna_pol_yheight;
MCNUM mccAna_pol_mosaic;
MCNUM mccAna_pol_dspread;
MCNUM mccAna_pol_Q;
MCNUM mccAna_pol_DM;
MCNUM mccAna_pol_pThreshold;
MCNUM mccAna_pol_Rup;
MCNUM mccAna_pol_Rdown;
int mccAna_pol_debug;

/* Definition parameters for component 'pollambdaMonitorDet' [18]. */
#define mccpollambdaMonitorDet_xwidth 0.1
#define mccpollambdaMonitorDet_yheight 0.1
#define mccpollambdaMonitorDet_nL 40
#define mccpollambdaMonitorDet_npol 11
#define mccpollambdaMonitorDet_restore_neutron 0
#define mccpollambdaMonitorDet_nowritefile 0
/* Setting parameters for component 'pollambdaMonitorDet' [18]. */
char mccpollambdaMonitorDet_filename[16384];
MCNUM mccpollambdaMonitorDet_mx;
MCNUM mccpollambdaMonitorDet_my;
MCNUM mccpollambdaMonitorDet_mz;
MCNUM mccpollambdaMonitorDet_Lmin;
MCNUM mccpollambdaMonitorDet_Lmax;

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
#line 8456 "Test_Pol_TripleAxis.c"
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
#line 60 "/usr/share/mcstas/2.5/sources/Source_simple.comp"
double pmul, srcArea;
int square;
double tx,ty,tz;
#line 8493 "Test_Pol_TripleAxis.c"
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

/* User declarations for component 'Mono_Arm' [3]. */
#define mccompcurname  Mono_Arm
#define mccompcurtype  Arm
#define mccompcurindex 3
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Mono_Out' [4]. */
#define mccompcurname  Mono_Out
#define mccompcurtype  Arm
#define mccompcurindex 4
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Mono' [5]. */
#define mccompcurname  Mono
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 5
#define mos_rms_y mccMono_mos_rms_y
#define mos_rms_z mccMono_mos_rms_z
#define mos_rms_max mccMono_mos_rms_max
#define mono_Q mccMono_mono_Q
#define zmin mccMono_zmin
#define zmax mccMono_zmax
#define ymin mccMono_ymin
#define ymax mccMono_ymax
#define zwidth mccMono_zwidth
#define yheight mccMono_yheight
#define mosaich mccMono_mosaich
#define mosaicv mccMono_mosaicv
#define r0 mccMono_r0
#define Q mccMono_Q
#define DM mccMono_DM
#line 95 "/usr/share/mcstas/2.5/optics/Monochromator_flat.comp"
  double mos_rms_y; /* root-mean-square of mosaic, in radians */
  double mos_rms_z;
  double mos_rms_max;
  double mono_Q;
#line 8554 "Test_Pol_TripleAxis.c"
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

/* User declarations for component 'Mono_pol' [6]. */
#define mccompcurname  Mono_pol
#define mccompcurtype  Monochromator_pol
#define mccompcurindex 6
#define mos_rms mccMono_pol_mos_rms
#define d_rms mccMono_pol_d_rms
#define mono_Q mccMono_pol_mono_Q
#define FN mccMono_pol_FN
#define FM mccMono_pol_FM
#define zwidth mccMono_pol_zwidth
#define yheight mccMono_pol_yheight
#define mosaic mccMono_pol_mosaic
#define dspread mccMono_pol_dspread
#define Q mccMono_pol_Q
#define DM mccMono_pol_DM
#define pThreshold mccMono_pol_pThreshold
#define Rup mccMono_pol_Rup
#define Rdown mccMono_pol_Rdown
#define debug mccMono_pol_debug
#line 100 "/usr/share/mcstas/2.5/optics/Monochromator_pol.comp"
double mos_rms; /* root-mean-square of mosaic, in radians */
double d_rms;   /* root-mean-square of d-spread, in AA */
double mono_Q;

double FN; /* Unit cell nuclear structure factor */
double FM; /* Unit cell magnetic structure factor */
#line 8600 "Test_Pol_TripleAxis.c"
#undef debug
#undef Rdown
#undef Rup
#undef pThreshold
#undef DM
#undef Q
#undef dspread
#undef mosaic
#undef yheight
#undef zwidth
#undef FM
#undef FN
#undef mono_Q
#undef d_rms
#undef mos_rms
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'pollambdaMonitorMono' [7]. */
#define mccompcurname  pollambdaMonitorMono
#define mccompcurtype  PolLambda_monitor
#define mccompcurindex 7
#define xwidth mccpollambdaMonitorMono_xwidth
#define yheight mccpollambdaMonitorMono_yheight
#define nL mccpollambdaMonitorMono_nL
#define npol mccpollambdaMonitorMono_npol
#define restore_neutron mccpollambdaMonitorMono_restore_neutron
#define nowritefile mccpollambdaMonitorMono_nowritefile
#define PolL_N mccpollambdaMonitorMono_PolL_N
#define PolL_p mccpollambdaMonitorMono_PolL_p
#define PolL_p2 mccpollambdaMonitorMono_PolL_p2
#define filename mccpollambdaMonitorMono_filename
#define mx mccpollambdaMonitorMono_mx
#define my mccpollambdaMonitorMono_my
#define mz mccpollambdaMonitorMono_mz
#define Lmin mccpollambdaMonitorMono_Lmin
#define Lmax mccpollambdaMonitorMono_Lmax
#line 56 "/usr/share/mcstas/2.5/monitors/PolLambda_monitor.comp"
double PolL_N[nL][npol];
double PolL_p[nL][npol];
double PolL_p2[nL][npol];
#line 8643 "Test_Pol_TripleAxis.c"
#undef Lmax
#undef Lmin
#undef mz
#undef my
#undef mx
#undef filename
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef npol
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Sample_Arm' [8]. */
#define mccompcurname  Sample_Arm
#define mccompcurtype  Arm
#define mccompcurindex 8
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Sample_Out' [9]. */
#define mccompcurname  Sample_Out
#define mccompcurtype  Arm
#define mccompcurindex 9
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Sample' [10]. */
#define mccompcurname  Sample
#define mccompcurtype  V_sample
#define mccompcurindex 10
#define VarsV mccSample_VarsV
#define radius mccSample_radius
#define thickness mccSample_thickness
#define zdepth mccSample_zdepth
#define Vc mccSample_Vc
#define sigma_abs mccSample_sigma_abs
#define sigma_inc mccSample_sigma_inc
#define radius_i mccSample_radius_i
#define radius_o mccSample_radius_o
#define h mccSample_h
#define focus_r mccSample_focus_r
#define pack mccSample_pack
#define frac mccSample_frac
#define f_QE mccSample_f_QE
#define gamma mccSample_gamma
#define target_x mccSample_target_x
#define target_y mccSample_target_y
#define target_z mccSample_target_z
#define focus_xw mccSample_focus_xw
#define focus_yh mccSample_focus_yh
#define focus_aw mccSample_focus_aw
#define focus_ah mccSample_focus_ah
#define xwidth mccSample_xwidth
#define yheight mccSample_yheight
#define zthick mccSample_zthick
#define rad_sphere mccSample_rad_sphere
#define sig_a mccSample_sig_a
#define sig_i mccSample_sig_i
#define V0 mccSample_V0
#define target_index mccSample_target_index
#define multiples mccSample_multiples
#line 117 "/usr/share/mcstas/2.5/obsolete/V_sample.comp"
  struct StructVarsV VarsV;
#line 8716 "Test_Pol_TripleAxis.c"
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

/* User declarations for component 'pollambdaMonitorSample' [11]. */
#define mccompcurname  pollambdaMonitorSample
#define mccompcurtype  PolLambda_monitor
#define mccompcurindex 11
#define xwidth mccpollambdaMonitorSample_xwidth
#define yheight mccpollambdaMonitorSample_yheight
#define nL mccpollambdaMonitorSample_nL
#define npol mccpollambdaMonitorSample_npol
#define restore_neutron mccpollambdaMonitorSample_restore_neutron
#define nowritefile mccpollambdaMonitorSample_nowritefile
#define PolL_N mccpollambdaMonitorSample_PolL_N
#define PolL_p mccpollambdaMonitorSample_PolL_p
#define PolL_p2 mccpollambdaMonitorSample_PolL_p2
#define filename mccpollambdaMonitorSample_filename
#define mx mccpollambdaMonitorSample_mx
#define my mccpollambdaMonitorSample_my
#define mz mccpollambdaMonitorSample_mz
#define Lmin mccpollambdaMonitorSample_Lmin
#define Lmax mccpollambdaMonitorSample_Lmax
#line 56 "/usr/share/mcstas/2.5/monitors/PolLambda_monitor.comp"
double PolL_N[nL][npol];
double PolL_p[nL][npol];
double PolL_p2[nL][npol];
#line 8775 "Test_Pol_TripleAxis.c"
#undef Lmax
#undef Lmin
#undef mz
#undef my
#undef mx
#undef filename
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef npol
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'msf' [12]. */
#define mccompcurname  msf
#define mccompcurtype  Pol_constBfield
#define mccompcurindex 12
#define xwidth mccmsf_xwidth
#define yheight mccmsf_yheight
#define length mccmsf_length
#define omegaL mccmsf_omegaL
#define B mccmsf_B
#define fliplambda mccmsf_fliplambda
#define flipangle mccmsf_flipangle
#line 64 "/usr/share/mcstas/2.5/optics/Pol_constBfield.comp"
  // Larmor frequency
  double omegaL;
#line 8809 "Test_Pol_TripleAxis.c"
#undef flipangle
#undef fliplambda
#undef B
#undef omegaL
#undef length
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'pollambdaMonitorMSF' [13]. */
#define mccompcurname  pollambdaMonitorMSF
#define mccompcurtype  PolLambda_monitor
#define mccompcurindex 13
#define xwidth mccpollambdaMonitorMSF_xwidth
#define yheight mccpollambdaMonitorMSF_yheight
#define nL mccpollambdaMonitorMSF_nL
#define npol mccpollambdaMonitorMSF_npol
#define restore_neutron mccpollambdaMonitorMSF_restore_neutron
#define nowritefile mccpollambdaMonitorMSF_nowritefile
#define PolL_N mccpollambdaMonitorMSF_PolL_N
#define PolL_p mccpollambdaMonitorMSF_PolL_p
#define PolL_p2 mccpollambdaMonitorMSF_PolL_p2
#define filename mccpollambdaMonitorMSF_filename
#define mx mccpollambdaMonitorMSF_mx
#define my mccpollambdaMonitorMSF_my
#define mz mccpollambdaMonitorMSF_mz
#define Lmin mccpollambdaMonitorMSF_Lmin
#define Lmax mccpollambdaMonitorMSF_Lmax
#line 56 "/usr/share/mcstas/2.5/monitors/PolLambda_monitor.comp"
double PolL_N[nL][npol];
double PolL_p[nL][npol];
double PolL_p2[nL][npol];
#line 8844 "Test_Pol_TripleAxis.c"
#undef Lmax
#undef Lmin
#undef mz
#undef my
#undef mx
#undef filename
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef npol
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Ana_Arm' [14]. */
#define mccompcurname  Ana_Arm
#define mccompcurtype  Arm
#define mccompcurindex 14
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Ana_Out' [15]. */
#define mccompcurname  Ana_Out
#define mccompcurtype  Arm
#define mccompcurindex 15
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Ana' [16]. */
#define mccompcurname  Ana
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 16
#define mos_rms_y mccAna_mos_rms_y
#define mos_rms_z mccAna_mos_rms_z
#define mos_rms_max mccAna_mos_rms_max
#define mono_Q mccAna_mono_Q
#define zmin mccAna_zmin
#define zmax mccAna_zmax
#define ymin mccAna_ymin
#define ymax mccAna_ymax
#define zwidth mccAna_zwidth
#define yheight mccAna_yheight
#define mosaich mccAna_mosaich
#define mosaicv mccAna_mosaicv
#define r0 mccAna_r0
#define Q mccAna_Q
#define DM mccAna_DM
#line 95 "/usr/share/mcstas/2.5/optics/Monochromator_flat.comp"
  double mos_rms_y; /* root-mean-square of mosaic, in radians */
  double mos_rms_z;
  double mos_rms_max;
  double mono_Q;
#line 8904 "Test_Pol_TripleAxis.c"
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

/* User declarations for component 'Ana_pol' [17]. */
#define mccompcurname  Ana_pol
#define mccompcurtype  Monochromator_pol
#define mccompcurindex 17
#define mos_rms mccAna_pol_mos_rms
#define d_rms mccAna_pol_d_rms
#define mono_Q mccAna_pol_mono_Q
#define FN mccAna_pol_FN
#define FM mccAna_pol_FM
#define zwidth mccAna_pol_zwidth
#define yheight mccAna_pol_yheight
#define mosaic mccAna_pol_mosaic
#define dspread mccAna_pol_dspread
#define Q mccAna_pol_Q
#define DM mccAna_pol_DM
#define pThreshold mccAna_pol_pThreshold
#define Rup mccAna_pol_Rup
#define Rdown mccAna_pol_Rdown
#define debug mccAna_pol_debug
#line 100 "/usr/share/mcstas/2.5/optics/Monochromator_pol.comp"
double mos_rms; /* root-mean-square of mosaic, in radians */
double d_rms;   /* root-mean-square of d-spread, in AA */
double mono_Q;

double FN; /* Unit cell nuclear structure factor */
double FM; /* Unit cell magnetic structure factor */
#line 8950 "Test_Pol_TripleAxis.c"
#undef debug
#undef Rdown
#undef Rup
#undef pThreshold
#undef DM
#undef Q
#undef dspread
#undef mosaic
#undef yheight
#undef zwidth
#undef FM
#undef FN
#undef mono_Q
#undef d_rms
#undef mos_rms
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'pollambdaMonitorDet' [18]. */
#define mccompcurname  pollambdaMonitorDet
#define mccompcurtype  PolLambda_monitor
#define mccompcurindex 18
#define xwidth mccpollambdaMonitorDet_xwidth
#define yheight mccpollambdaMonitorDet_yheight
#define nL mccpollambdaMonitorDet_nL
#define npol mccpollambdaMonitorDet_npol
#define restore_neutron mccpollambdaMonitorDet_restore_neutron
#define nowritefile mccpollambdaMonitorDet_nowritefile
#define PolL_N mccpollambdaMonitorDet_PolL_N
#define PolL_p mccpollambdaMonitorDet_PolL_p
#define PolL_p2 mccpollambdaMonitorDet_PolL_p2
#define filename mccpollambdaMonitorDet_filename
#define mx mccpollambdaMonitorDet_mx
#define my mccpollambdaMonitorDet_my
#define mz mccpollambdaMonitorDet_mz
#define Lmin mccpollambdaMonitorDet_Lmin
#define Lmax mccpollambdaMonitorDet_Lmax
#line 56 "/usr/share/mcstas/2.5/monitors/PolLambda_monitor.comp"
double PolL_N[nL][npol];
double PolL_p[nL][npol];
double PolL_p2[nL][npol];
#line 8993 "Test_Pol_TripleAxis.c"
#undef Lmax
#undef Lmin
#undef mz
#undef my
#undef mx
#undef filename
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef npol
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

Coords mcposaOrigin, mcposrOrigin;
Rotation mcrotaOrigin, mcrotrOrigin;
Coords mcposaSource, mcposrSource;
Rotation mcrotaSource, mcrotrSource;
Coords mcposaMono_Arm, mcposrMono_Arm;
Rotation mcrotaMono_Arm, mcrotrMono_Arm;
Coords mcposaMono_Out, mcposrMono_Out;
Rotation mcrotaMono_Out, mcrotrMono_Out;
Coords mcposaMono, mcposrMono;
Rotation mcrotaMono, mcrotrMono;
Coords mcposaMono_pol, mcposrMono_pol;
Rotation mcrotaMono_pol, mcrotrMono_pol;
Coords mcposapollambdaMonitorMono, mcposrpollambdaMonitorMono;
Rotation mcrotapollambdaMonitorMono, mcrotrpollambdaMonitorMono;
Coords mcposaSample_Arm, mcposrSample_Arm;
Rotation mcrotaSample_Arm, mcrotrSample_Arm;
Coords mcposaSample_Out, mcposrSample_Out;
Rotation mcrotaSample_Out, mcrotrSample_Out;
Coords mcposaSample, mcposrSample;
Rotation mcrotaSample, mcrotrSample;
Coords mcposapollambdaMonitorSample, mcposrpollambdaMonitorSample;
Rotation mcrotapollambdaMonitorSample, mcrotrpollambdaMonitorSample;
Coords mcposamsf, mcposrmsf;
Rotation mcrotamsf, mcrotrmsf;
Coords mcposapollambdaMonitorMSF, mcposrpollambdaMonitorMSF;
Rotation mcrotapollambdaMonitorMSF, mcrotrpollambdaMonitorMSF;
Coords mcposaAna_Arm, mcposrAna_Arm;
Rotation mcrotaAna_Arm, mcrotrAna_Arm;
Coords mcposaAna_Out, mcposrAna_Out;
Rotation mcrotaAna_Out, mcrotrAna_Out;
Coords mcposaAna, mcposrAna;
Rotation mcrotaAna, mcrotrAna;
Coords mcposaAna_pol, mcposrAna_pol;
Rotation mcrotaAna_pol, mcrotrAna_pol;
Coords mcposapollambdaMonitorDet, mcposrpollambdaMonitorDet;
Rotation mcrotapollambdaMonitorDet, mcrotrpollambdaMonitorDet;

MCNUM mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz, mcnt, mcnsx, mcnsy, mcnsz, mcnp;

/* end declare */

void mcinit(void) {
#define mccompcurname  Test_Pol_TripleAxis
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaTest_Pol_TripleAxis coords_set(0,0,0)
#define OPTION mcipOPTION
#define LAMBDA mcipLAMBDA
#define MOZ mcipMOZ
#line 57 "Test_Pol_TripleAxis.instr"
{
  int    ORDER = 1;
  double vi, Ki, Kf;
  int    SM, SS, SA;

  if(OPTION<-1 || OPTION>1) {

    printf("OPTION %d not permitted.", OPTION);
    exit(1);
  }

  /* SM : scattering at mono to the right (-1)/left(+1) */
  /* SS : scattering at sample to the right (-1)/left(+1) */
  /* SA : scattering at analyser to the right (-1)/left(+1) */
  SM = 1; SS = -1; SA = 1;

  mono_q = 2*PI*ORDER/DM;  /* Q mono in Angs-1 */

  Ki = 2*PI/LAMBDA;
  vi = K2V*fabs(Ki);
  Ei = VS2E*vi*vi;

  Kf=Ki;

  A2 = asin(mono_q/2/Ki)*RAD2DEG*2;
  A6 = asin(mono_q/2/Kf)*RAD2DEG*2;
  A4 = A2;

  A2 *= SM;   /* A1 : mono theta (crystal) */
  A1 = A2/2;  /* A2 : mono 2 theta (arm to sample) */
  A4 *= SS;   /* A3 : sample theta */
  A3 = A4/2;  /* A4 : sample 2 theta (arm to analyser) */
  A6 *= SA;   /* A5 : analyser theta (crystal) */
  A5 = A6/2;  /* A6 : analyser 2 theta (arm to Dector) */
}
#line 9098 "Test_Pol_TripleAxis.c"
#undef MOZ
#undef LAMBDA
#undef OPTION
#undef mcposaTest_Pol_TripleAxis
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
#line 39 "Test_Pol_TripleAxis.instr"
  if("NULL") strncpy(mccOrigin_profile, "NULL" ? "NULL" : "", 16384); else mccOrigin_profile[0]='\0';
#line 39 "Test_Pol_TripleAxis.instr"
  mccOrigin_percent = 10;
#line 39 "Test_Pol_TripleAxis.instr"
  mccOrigin_flag_save = 0;
#line 39 "Test_Pol_TripleAxis.instr"
  mccOrigin_minutes = 0;
#line 9127 "Test_Pol_TripleAxis.c"

  SIG_MESSAGE("Origin (Init:Place/Rotate)");
  rot_set_rotation(mcrotaOrigin,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9134 "Test_Pol_TripleAxis.c"
  rot_copy(mcrotrOrigin, mcrotaOrigin);
  mcposaOrigin = coords_set(
#line 98 "Test_Pol_TripleAxis.instr"
    0,
#line 98 "Test_Pol_TripleAxis.instr"
    0,
#line 98 "Test_Pol_TripleAxis.instr"
    0);
#line 9143 "Test_Pol_TripleAxis.c"
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
#line 102 "Test_Pol_TripleAxis.instr"
  mccSource_radius = 0.01;
#line 52 "Test_Pol_TripleAxis.instr"
  mccSource_yheight = 0;
#line 52 "Test_Pol_TripleAxis.instr"
  mccSource_xwidth = 0;
#line 102 "Test_Pol_TripleAxis.instr"
  mccSource_dist = 5.0;
#line 103 "Test_Pol_TripleAxis.instr"
  mccSource_focus_xw = 0.02;
#line 103 "Test_Pol_TripleAxis.instr"
  mccSource_focus_yh = 0.02;
#line 104 "Test_Pol_TripleAxis.instr"
  mccSource_E0 = Ei;
#line 104 "Test_Pol_TripleAxis.instr"
  mccSource_dE = 0.5;
#line 54 "Test_Pol_TripleAxis.instr"
  mccSource_lambda0 = 0;
#line 54 "Test_Pol_TripleAxis.instr"
  mccSource_dlambda = 0;
#line 55 "Test_Pol_TripleAxis.instr"
  mccSource_flux = 1;
#line 55 "Test_Pol_TripleAxis.instr"
  mccSource_gauss = 0;
#line 55 "Test_Pol_TripleAxis.instr"
  mccSource_target_index = + 1;
#line 9180 "Test_Pol_TripleAxis.c"

  SIG_MESSAGE("Source (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9187 "Test_Pol_TripleAxis.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaSource);
  rot_transpose(mcrotaOrigin, mctr1);
  rot_mul(mcrotaSource, mctr1, mcrotrSource);
  mctc1 = coords_set(
#line 105 "Test_Pol_TripleAxis.instr"
    0,
#line 105 "Test_Pol_TripleAxis.instr"
    0,
#line 105 "Test_Pol_TripleAxis.instr"
    0);
#line 9198 "Test_Pol_TripleAxis.c"
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
    /* Component Mono_Arm. */
  /* Setting parameters for component Mono_Arm. */
  SIG_MESSAGE("Mono_Arm (Init:SetPar)");

  SIG_MESSAGE("Mono_Arm (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 108 "Test_Pol_TripleAxis.instr"
    (0)*DEG2RAD,
#line 108 "Test_Pol_TripleAxis.instr"
    (A1)*DEG2RAD,
#line 108 "Test_Pol_TripleAxis.instr"
    (0)*DEG2RAD);
#line 9221 "Test_Pol_TripleAxis.c"
  rot_mul(mctr1, mcrotaSource, mcrotaMono_Arm);
  rot_transpose(mcrotaSource, mctr1);
  rot_mul(mcrotaMono_Arm, mctr1, mcrotrMono_Arm);
  mctc1 = coords_set(
#line 108 "Test_Pol_TripleAxis.instr"
    0,
#line 108 "Test_Pol_TripleAxis.instr"
    0,
#line 108 "Test_Pol_TripleAxis.instr"
    5.0);
#line 9232 "Test_Pol_TripleAxis.c"
  rot_transpose(mcrotaSource, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMono_Arm = coords_add(mcposaSource, mctc2);
  mctc1 = coords_sub(mcposaSource, mcposaMono_Arm);
  mcposrMono_Arm = rot_apply(mcrotaMono_Arm, mctc1);
  mcDEBUG_COMPONENT("Mono_Arm", mcposaMono_Arm, mcrotaMono_Arm)
  mccomp_posa[3] = mcposaMono_Arm;
  mccomp_posr[3] = mcposrMono_Arm;
  mcNCounter[3]  = mcPCounter[3] = mcP2Counter[3] = 0;
  mcAbsorbProp[3]= 0;
    /* Component Mono_Out. */
  /* Setting parameters for component Mono_Out. */
  SIG_MESSAGE("Mono_Out (Init:SetPar)");

  SIG_MESSAGE("Mono_Out (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 111 "Test_Pol_TripleAxis.instr"
    (0)*DEG2RAD,
#line 111 "Test_Pol_TripleAxis.instr"
    (A2)*DEG2RAD,
#line 111 "Test_Pol_TripleAxis.instr"
    (0)*DEG2RAD);
#line 9255 "Test_Pol_TripleAxis.c"
  rot_mul(mctr1, mcrotaSource, mcrotaMono_Out);
  rot_transpose(mcrotaMono_Arm, mctr1);
  rot_mul(mcrotaMono_Out, mctr1, mcrotrMono_Out);
  mctc1 = coords_set(
#line 111 "Test_Pol_TripleAxis.instr"
    0,
#line 111 "Test_Pol_TripleAxis.instr"
    0,
#line 111 "Test_Pol_TripleAxis.instr"
    0);
#line 9266 "Test_Pol_TripleAxis.c"
  rot_transpose(mcrotaMono_Arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMono_Out = coords_add(mcposaMono_Arm, mctc2);
  mctc1 = coords_sub(mcposaMono_Arm, mcposaMono_Out);
  mcposrMono_Out = rot_apply(mcrotaMono_Out, mctc1);
  mcDEBUG_COMPONENT("Mono_Out", mcposaMono_Out, mcrotaMono_Out)
  mccomp_posa[4] = mcposaMono_Out;
  mccomp_posr[4] = mcposrMono_Out;
  mcNCounter[4]  = mcPCounter[4] = mcP2Counter[4] = 0;
  mcAbsorbProp[4]= 0;
    /* Component Mono. */
  /* Setting parameters for component Mono. */
  SIG_MESSAGE("Mono (Init:SetPar)");
#line 64 "Test_Pol_TripleAxis.instr"
  mccMono_zmin = -0.05;
#line 64 "Test_Pol_TripleAxis.instr"
  mccMono_zmax = 0.05;
#line 64 "Test_Pol_TripleAxis.instr"
  mccMono_ymin = -0.05;
#line 64 "Test_Pol_TripleAxis.instr"
  mccMono_ymax = 0.05;
#line 114 "Test_Pol_TripleAxis.instr"
  mccMono_zwidth = 0.10;
#line 114 "Test_Pol_TripleAxis.instr"
  mccMono_yheight = 0.08;
#line 115 "Test_Pol_TripleAxis.instr"
  mccMono_mosaich = mcipMOZ;
#line 115 "Test_Pol_TripleAxis.instr"
  mccMono_mosaicv = mcipMOZ;
#line 116 "Test_Pol_TripleAxis.instr"
  mccMono_r0 = 1.0;
#line 116 "Test_Pol_TripleAxis.instr"
  mccMono_Q = mono_q;
#line 66 "Test_Pol_TripleAxis.instr"
  mccMono_DM = 0;
#line 9302 "Test_Pol_TripleAxis.c"

  SIG_MESSAGE("Mono (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9309 "Test_Pol_TripleAxis.c"
  rot_mul(mctr1, mcrotaMono_Arm, mcrotaMono);
  rot_transpose(mcrotaMono_Out, mctr1);
  rot_mul(mcrotaMono, mctr1, mcrotrMono);
  mctc1 = coords_set(
#line 117 "Test_Pol_TripleAxis.instr"
    0,
#line 117 "Test_Pol_TripleAxis.instr"
    0,
#line 117 "Test_Pol_TripleAxis.instr"
    0);
#line 9320 "Test_Pol_TripleAxis.c"
  rot_transpose(mcrotaMono_Arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMono = coords_add(mcposaMono_Arm, mctc2);
  mctc1 = coords_sub(mcposaMono_Out, mcposaMono);
  mcposrMono = rot_apply(mcrotaMono, mctc1);
  mcDEBUG_COMPONENT("Mono", mcposaMono, mcrotaMono)
  mccomp_posa[5] = mcposaMono;
  mccomp_posr[5] = mcposrMono;
  mcNCounter[5]  = mcPCounter[5] = mcP2Counter[5] = 0;
  mcAbsorbProp[5]= 0;
    /* Component Mono_pol. */
  /* Setting parameters for component Mono_pol. */
  SIG_MESSAGE("Mono_pol (Init:SetPar)");
#line 120 "Test_Pol_TripleAxis.instr"
  mccMono_pol_zwidth = 0.10;
#line 120 "Test_Pol_TripleAxis.instr"
  mccMono_pol_yheight = 0.08;
#line 121 "Test_Pol_TripleAxis.instr"
  mccMono_pol_mosaic = mcipMOZ;
#line 121 "Test_Pol_TripleAxis.instr"
  mccMono_pol_dspread = 0.0;
#line 122 "Test_Pol_TripleAxis.instr"
  mccMono_pol_Q = mono_q;
#line 89 "Test_Pol_TripleAxis.instr"
  mccMono_pol_DM = 0;
#line 89 "Test_Pol_TripleAxis.instr"
  mccMono_pol_pThreshold = 0;
#line 122 "Test_Pol_TripleAxis.instr"
  mccMono_pol_Rup = 1.0;
#line 122 "Test_Pol_TripleAxis.instr"
  mccMono_pol_Rdown = 0.0;
#line 122 "Test_Pol_TripleAxis.instr"
  mccMono_pol_debug = 0;
#line 9354 "Test_Pol_TripleAxis.c"

  SIG_MESSAGE("Mono_pol (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9361 "Test_Pol_TripleAxis.c"
  rot_mul(mctr1, mcrotaMono_Arm, mcrotaMono_pol);
  rot_transpose(mcrotaMono, mctr1);
  rot_mul(mcrotaMono_pol, mctr1, mcrotrMono_pol);
  mctc1 = coords_set(
#line 123 "Test_Pol_TripleAxis.instr"
    0,
#line 123 "Test_Pol_TripleAxis.instr"
    0,
#line 123 "Test_Pol_TripleAxis.instr"
    0);
#line 9372 "Test_Pol_TripleAxis.c"
  rot_transpose(mcrotaMono_Arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMono_pol = coords_add(mcposaMono_Arm, mctc2);
  mctc1 = coords_sub(mcposaMono, mcposaMono_pol);
  mcposrMono_pol = rot_apply(mcrotaMono_pol, mctc1);
  mcDEBUG_COMPONENT("Mono_pol", mcposaMono_pol, mcrotaMono_pol)
  mccomp_posa[6] = mcposaMono_pol;
  mccomp_posr[6] = mcposrMono_pol;
  mcNCounter[6]  = mcPCounter[6] = mcP2Counter[6] = 0;
  mcAbsorbProp[6]= 0;
    /* Component pollambdaMonitorMono. */
  /* Setting parameters for component pollambdaMonitorMono. */
  SIG_MESSAGE("pollambdaMonitorMono (Init:SetPar)");
#line 129 "Test_Pol_TripleAxis.instr"
  if("pollambdaMonMono.data") strncpy(mccpollambdaMonitorMono_filename, "pollambdaMonMono.data" ? "pollambdaMonMono.data" : "", 16384); else mccpollambdaMonitorMono_filename[0]='\0';
#line 50 "Test_Pol_TripleAxis.instr"
  mccpollambdaMonitorMono_mx = 0;
#line 129 "Test_Pol_TripleAxis.instr"
  mccpollambdaMonitorMono_my = 1;
#line 50 "Test_Pol_TripleAxis.instr"
  mccpollambdaMonitorMono_mz = 0;
#line 128 "Test_Pol_TripleAxis.instr"
  mccpollambdaMonitorMono_Lmin = 0.97 * mcipLAMBDA;
#line 128 "Test_Pol_TripleAxis.instr"
  mccpollambdaMonitorMono_Lmax = 1.03 * mcipLAMBDA;
#line 9398 "Test_Pol_TripleAxis.c"

  SIG_MESSAGE("pollambdaMonitorMono (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9405 "Test_Pol_TripleAxis.c"
  rot_mul(mctr1, mcrotaMono_Out, mcrotapollambdaMonitorMono);
  rot_transpose(mcrotaMono_pol, mctr1);
  rot_mul(mcrotapollambdaMonitorMono, mctr1, mcrotrpollambdaMonitorMono);
  mctc1 = coords_set(
#line 130 "Test_Pol_TripleAxis.instr"
    0,
#line 130 "Test_Pol_TripleAxis.instr"
    0,
#line 130 "Test_Pol_TripleAxis.instr"
    0.5);
#line 9416 "Test_Pol_TripleAxis.c"
  rot_transpose(mcrotaMono_Out, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapollambdaMonitorMono = coords_add(mcposaMono_Out, mctc2);
  mctc1 = coords_sub(mcposaMono_pol, mcposapollambdaMonitorMono);
  mcposrpollambdaMonitorMono = rot_apply(mcrotapollambdaMonitorMono, mctc1);
  mcDEBUG_COMPONENT("pollambdaMonitorMono", mcposapollambdaMonitorMono, mcrotapollambdaMonitorMono)
  mccomp_posa[7] = mcposapollambdaMonitorMono;
  mccomp_posr[7] = mcposrpollambdaMonitorMono;
  mcNCounter[7]  = mcPCounter[7] = mcP2Counter[7] = 0;
  mcAbsorbProp[7]= 0;
    /* Component Sample_Arm. */
  /* Setting parameters for component Sample_Arm. */
  SIG_MESSAGE("Sample_Arm (Init:SetPar)");

  SIG_MESSAGE("Sample_Arm (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 133 "Test_Pol_TripleAxis.instr"
    (0)*DEG2RAD,
#line 133 "Test_Pol_TripleAxis.instr"
    (A3)*DEG2RAD,
#line 133 "Test_Pol_TripleAxis.instr"
    (0)*DEG2RAD);
#line 9439 "Test_Pol_TripleAxis.c"
  rot_mul(mctr1, mcrotaMono_Out, mcrotaSample_Arm);
  rot_transpose(mcrotapollambdaMonitorMono, mctr1);
  rot_mul(mcrotaSample_Arm, mctr1, mcrotrSample_Arm);
  mctc1 = coords_set(
#line 133 "Test_Pol_TripleAxis.instr"
    0,
#line 133 "Test_Pol_TripleAxis.instr"
    0,
#line 133 "Test_Pol_TripleAxis.instr"
    1.0);
#line 9450 "Test_Pol_TripleAxis.c"
  rot_transpose(mcrotaMono_Out, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaSample_Arm = coords_add(mcposaMono_Out, mctc2);
  mctc1 = coords_sub(mcposapollambdaMonitorMono, mcposaSample_Arm);
  mcposrSample_Arm = rot_apply(mcrotaSample_Arm, mctc1);
  mcDEBUG_COMPONENT("Sample_Arm", mcposaSample_Arm, mcrotaSample_Arm)
  mccomp_posa[8] = mcposaSample_Arm;
  mccomp_posr[8] = mcposrSample_Arm;
  mcNCounter[8]  = mcPCounter[8] = mcP2Counter[8] = 0;
  mcAbsorbProp[8]= 0;
    /* Component Sample_Out. */
  /* Setting parameters for component Sample_Out. */
  SIG_MESSAGE("Sample_Out (Init:SetPar)");

  SIG_MESSAGE("Sample_Out (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 136 "Test_Pol_TripleAxis.instr"
    (0)*DEG2RAD,
#line 136 "Test_Pol_TripleAxis.instr"
    (A4)*DEG2RAD,
#line 136 "Test_Pol_TripleAxis.instr"
    (0)*DEG2RAD);
#line 9473 "Test_Pol_TripleAxis.c"
  rot_mul(mctr1, mcrotaMono_Out, mcrotaSample_Out);
  rot_transpose(mcrotaSample_Arm, mctr1);
  rot_mul(mcrotaSample_Out, mctr1, mcrotrSample_Out);
  mctc1 = coords_set(
#line 136 "Test_Pol_TripleAxis.instr"
    0,
#line 136 "Test_Pol_TripleAxis.instr"
    0,
#line 136 "Test_Pol_TripleAxis.instr"
    0);
#line 9484 "Test_Pol_TripleAxis.c"
  rot_transpose(mcrotaSample_Arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaSample_Out = coords_add(mcposaSample_Arm, mctc2);
  mctc1 = coords_sub(mcposaSample_Arm, mcposaSample_Out);
  mcposrSample_Out = rot_apply(mcrotaSample_Out, mctc1);
  mcDEBUG_COMPONENT("Sample_Out", mcposaSample_Out, mcrotaSample_Out)
  mccomp_posa[9] = mcposaSample_Out;
  mccomp_posr[9] = mcposrSample_Out;
  mcNCounter[9]  = mcPCounter[9] = mcP2Counter[9] = 0;
  mcAbsorbProp[9]= 0;
    /* Component Sample. */
  /* Setting parameters for component Sample. */
  SIG_MESSAGE("Sample (Init:SetPar)");
#line 139 "Test_Pol_TripleAxis.instr"
  mccSample_radius = 0.024;
#line 91 "Test_Pol_TripleAxis.instr"
  mccSample_thickness = 0;
#line 91 "Test_Pol_TripleAxis.instr"
  mccSample_zdepth = 0;
#line 91 "Test_Pol_TripleAxis.instr"
  mccSample_Vc = 13.827;
#line 91 "Test_Pol_TripleAxis.instr"
  mccSample_sigma_abs = 5.08;
#line 91 "Test_Pol_TripleAxis.instr"
  mccSample_sigma_inc = 5.08;
#line 92 "Test_Pol_TripleAxis.instr"
  mccSample_radius_i = 0;
#line 92 "Test_Pol_TripleAxis.instr"
  mccSample_radius_o = 0;
#line 92 "Test_Pol_TripleAxis.instr"
  mccSample_h = 0;
#line 92 "Test_Pol_TripleAxis.instr"
  mccSample_focus_r = 0;
#line 141 "Test_Pol_TripleAxis.instr"
  mccSample_pack = 1;
#line 92 "Test_Pol_TripleAxis.instr"
  mccSample_frac = 1;
#line 92 "Test_Pol_TripleAxis.instr"
  mccSample_f_QE = 0;
#line 92 "Test_Pol_TripleAxis.instr"
  mccSample_gamma = 0;
#line 140 "Test_Pol_TripleAxis.instr"
  mccSample_target_x = 0;
#line 140 "Test_Pol_TripleAxis.instr"
  mccSample_target_y = 0;
#line 140 "Test_Pol_TripleAxis.instr"
  mccSample_target_z = 1.0;
#line 141 "Test_Pol_TripleAxis.instr"
  mccSample_focus_xw = 0.01;
#line 141 "Test_Pol_TripleAxis.instr"
  mccSample_focus_yh = 0.01;
#line 94 "Test_Pol_TripleAxis.instr"
  mccSample_focus_aw = 0;
#line 94 "Test_Pol_TripleAxis.instr"
  mccSample_focus_ah = 0;
#line 94 "Test_Pol_TripleAxis.instr"
  mccSample_xwidth = 0;
#line 139 "Test_Pol_TripleAxis.instr"
  mccSample_yheight = 0.0254;
#line 94 "Test_Pol_TripleAxis.instr"
  mccSample_zthick = 0;
#line 94 "Test_Pol_TripleAxis.instr"
  mccSample_rad_sphere = 0;
#line 94 "Test_Pol_TripleAxis.instr"
  mccSample_sig_a = 0;
#line 94 "Test_Pol_TripleAxis.instr"
  mccSample_sig_i = 0;
#line 94 "Test_Pol_TripleAxis.instr"
  mccSample_V0 = 0;
#line 94 "Test_Pol_TripleAxis.instr"
  mccSample_target_index = 0;
#line 94 "Test_Pol_TripleAxis.instr"
  mccSample_multiples = 1;
#line 9558 "Test_Pol_TripleAxis.c"

  SIG_MESSAGE("Sample (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9565 "Test_Pol_TripleAxis.c"
  rot_mul(mctr1, mcrotaSample_Out, mcrotaSample);
  rot_transpose(mcrotaSample_Out, mctr1);
  rot_mul(mcrotaSample, mctr1, mcrotrSample);
  mctc1 = coords_set(
#line 142 "Test_Pol_TripleAxis.instr"
    0,
#line 142 "Test_Pol_TripleAxis.instr"
    0,
#line 142 "Test_Pol_TripleAxis.instr"
    0);
#line 9576 "Test_Pol_TripleAxis.c"
  rot_transpose(mcrotaSample_Out, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaSample = coords_add(mcposaSample_Out, mctc2);
  mctc1 = coords_sub(mcposaSample_Out, mcposaSample);
  mcposrSample = rot_apply(mcrotaSample, mctc1);
  mcDEBUG_COMPONENT("Sample", mcposaSample, mcrotaSample)
  mccomp_posa[10] = mcposaSample;
  mccomp_posr[10] = mcposrSample;
  mcNCounter[10]  = mcPCounter[10] = mcP2Counter[10] = 0;
  mcAbsorbProp[10]= 0;
    /* Component pollambdaMonitorSample. */
  /* Setting parameters for component pollambdaMonitorSample. */
  SIG_MESSAGE("pollambdaMonitorSample (Init:SetPar)");
#line 147 "Test_Pol_TripleAxis.instr"
  if("pollambdaMonSample.data") strncpy(mccpollambdaMonitorSample_filename, "pollambdaMonSample.data" ? "pollambdaMonSample.data" : "", 16384); else mccpollambdaMonitorSample_filename[0]='\0';
#line 50 "Test_Pol_TripleAxis.instr"
  mccpollambdaMonitorSample_mx = 0;
#line 147 "Test_Pol_TripleAxis.instr"
  mccpollambdaMonitorSample_my = 1;
#line 50 "Test_Pol_TripleAxis.instr"
  mccpollambdaMonitorSample_mz = 0;
#line 146 "Test_Pol_TripleAxis.instr"
  mccpollambdaMonitorSample_Lmin = 0.97 * mcipLAMBDA;
#line 146 "Test_Pol_TripleAxis.instr"
  mccpollambdaMonitorSample_Lmax = 1.03 * mcipLAMBDA;
#line 9602 "Test_Pol_TripleAxis.c"

  SIG_MESSAGE("pollambdaMonitorSample (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9609 "Test_Pol_TripleAxis.c"
  rot_mul(mctr1, mcrotaSample_Out, mcrotapollambdaMonitorSample);
  rot_transpose(mcrotaSample, mctr1);
  rot_mul(mcrotapollambdaMonitorSample, mctr1, mcrotrpollambdaMonitorSample);
  mctc1 = coords_set(
#line 148 "Test_Pol_TripleAxis.instr"
    0,
#line 148 "Test_Pol_TripleAxis.instr"
    0,
#line 148 "Test_Pol_TripleAxis.instr"
    0.25);
#line 9620 "Test_Pol_TripleAxis.c"
  rot_transpose(mcrotaSample_Out, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapollambdaMonitorSample = coords_add(mcposaSample_Out, mctc2);
  mctc1 = coords_sub(mcposaSample, mcposapollambdaMonitorSample);
  mcposrpollambdaMonitorSample = rot_apply(mcrotapollambdaMonitorSample, mctc1);
  mcDEBUG_COMPONENT("pollambdaMonitorSample", mcposapollambdaMonitorSample, mcrotapollambdaMonitorSample)
  mccomp_posa[11] = mcposapollambdaMonitorSample;
  mccomp_posr[11] = mcposrpollambdaMonitorSample;
  mcNCounter[11]  = mcPCounter[11] = mcP2Counter[11] = 0;
  mcAbsorbProp[11]= 0;
    /* Component msf. */
  /* Setting parameters for component msf. */
  SIG_MESSAGE("msf (Init:SetPar)");
#line 43 "Test_Pol_TripleAxis.instr"
  mccmsf_B = 0;
#line 151 "Test_Pol_TripleAxis.instr"
  mccmsf_fliplambda = mcipLAMBDA;
#line 43 "Test_Pol_TripleAxis.instr"
  mccmsf_flipangle = 180;
#line 9640 "Test_Pol_TripleAxis.c"

  SIG_MESSAGE("msf (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 153 "Test_Pol_TripleAxis.instr"
    (0)*DEG2RAD,
#line 153 "Test_Pol_TripleAxis.instr"
    (0)*DEG2RAD,
#line 153 "Test_Pol_TripleAxis.instr"
    (90)*DEG2RAD);
#line 9650 "Test_Pol_TripleAxis.c"
  rot_mul(mctr1, mcrotaSample_Out, mcrotamsf);
  rot_transpose(mcrotapollambdaMonitorSample, mctr1);
  rot_mul(mcrotamsf, mctr1, mcrotrmsf);
  mctc1 = coords_set(
#line 152 "Test_Pol_TripleAxis.instr"
    0,
#line 152 "Test_Pol_TripleAxis.instr"
    0,
#line 152 "Test_Pol_TripleAxis.instr"
    0.4);
#line 9661 "Test_Pol_TripleAxis.c"
  rot_transpose(mcrotaSample_Out, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposamsf = coords_add(mcposaSample_Out, mctc2);
  mctc1 = coords_sub(mcposapollambdaMonitorSample, mcposamsf);
  mcposrmsf = rot_apply(mcrotamsf, mctc1);
  mcDEBUG_COMPONENT("msf", mcposamsf, mcrotamsf)
  mccomp_posa[12] = mcposamsf;
  mccomp_posr[12] = mcposrmsf;
  mcNCounter[12]  = mcPCounter[12] = mcP2Counter[12] = 0;
  mcAbsorbProp[12]= 0;
    /* Component pollambdaMonitorMSF. */
  /* Setting parameters for component pollambdaMonitorMSF. */
  SIG_MESSAGE("pollambdaMonitorMSF (Init:SetPar)");
#line 158 "Test_Pol_TripleAxis.instr"
  if("pollambdaMonMsf.data") strncpy(mccpollambdaMonitorMSF_filename, "pollambdaMonMsf.data" ? "pollambdaMonMsf.data" : "", 16384); else mccpollambdaMonitorMSF_filename[0]='\0';
#line 50 "Test_Pol_TripleAxis.instr"
  mccpollambdaMonitorMSF_mx = 0;
#line 158 "Test_Pol_TripleAxis.instr"
  mccpollambdaMonitorMSF_my = 1;
#line 50 "Test_Pol_TripleAxis.instr"
  mccpollambdaMonitorMSF_mz = 0;
#line 157 "Test_Pol_TripleAxis.instr"
  mccpollambdaMonitorMSF_Lmin = 0.97 * mcipLAMBDA;
#line 157 "Test_Pol_TripleAxis.instr"
  mccpollambdaMonitorMSF_Lmax = 1.03 * mcipLAMBDA;
#line 9687 "Test_Pol_TripleAxis.c"

  SIG_MESSAGE("pollambdaMonitorMSF (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9694 "Test_Pol_TripleAxis.c"
  rot_mul(mctr1, mcrotaSample_Out, mcrotapollambdaMonitorMSF);
  rot_transpose(mcrotamsf, mctr1);
  rot_mul(mcrotapollambdaMonitorMSF, mctr1, mcrotrpollambdaMonitorMSF);
  mctc1 = coords_set(
#line 159 "Test_Pol_TripleAxis.instr"
    0,
#line 159 "Test_Pol_TripleAxis.instr"
    0,
#line 159 "Test_Pol_TripleAxis.instr"
    0.75);
#line 9705 "Test_Pol_TripleAxis.c"
  rot_transpose(mcrotaSample_Out, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapollambdaMonitorMSF = coords_add(mcposaSample_Out, mctc2);
  mctc1 = coords_sub(mcposamsf, mcposapollambdaMonitorMSF);
  mcposrpollambdaMonitorMSF = rot_apply(mcrotapollambdaMonitorMSF, mctc1);
  mcDEBUG_COMPONENT("pollambdaMonitorMSF", mcposapollambdaMonitorMSF, mcrotapollambdaMonitorMSF)
  mccomp_posa[13] = mcposapollambdaMonitorMSF;
  mccomp_posr[13] = mcposrpollambdaMonitorMSF;
  mcNCounter[13]  = mcPCounter[13] = mcP2Counter[13] = 0;
  mcAbsorbProp[13]= 0;
    /* Component Ana_Arm. */
  /* Setting parameters for component Ana_Arm. */
  SIG_MESSAGE("Ana_Arm (Init:SetPar)");

  SIG_MESSAGE("Ana_Arm (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 162 "Test_Pol_TripleAxis.instr"
    (0)*DEG2RAD,
#line 162 "Test_Pol_TripleAxis.instr"
    (A5)*DEG2RAD,
#line 162 "Test_Pol_TripleAxis.instr"
    (0)*DEG2RAD);
#line 9728 "Test_Pol_TripleAxis.c"
  rot_mul(mctr1, mcrotaSample_Out, mcrotaAna_Arm);
  rot_transpose(mcrotapollambdaMonitorMSF, mctr1);
  rot_mul(mcrotaAna_Arm, mctr1, mcrotrAna_Arm);
  mctc1 = coords_set(
#line 162 "Test_Pol_TripleAxis.instr"
    0,
#line 162 "Test_Pol_TripleAxis.instr"
    0,
#line 162 "Test_Pol_TripleAxis.instr"
    1.0);
#line 9739 "Test_Pol_TripleAxis.c"
  rot_transpose(mcrotaSample_Out, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaAna_Arm = coords_add(mcposaSample_Out, mctc2);
  mctc1 = coords_sub(mcposapollambdaMonitorMSF, mcposaAna_Arm);
  mcposrAna_Arm = rot_apply(mcrotaAna_Arm, mctc1);
  mcDEBUG_COMPONENT("Ana_Arm", mcposaAna_Arm, mcrotaAna_Arm)
  mccomp_posa[14] = mcposaAna_Arm;
  mccomp_posr[14] = mcposrAna_Arm;
  mcNCounter[14]  = mcPCounter[14] = mcP2Counter[14] = 0;
  mcAbsorbProp[14]= 0;
    /* Component Ana_Out. */
  /* Setting parameters for component Ana_Out. */
  SIG_MESSAGE("Ana_Out (Init:SetPar)");

  SIG_MESSAGE("Ana_Out (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 165 "Test_Pol_TripleAxis.instr"
    (0)*DEG2RAD,
#line 165 "Test_Pol_TripleAxis.instr"
    (A6)*DEG2RAD,
#line 165 "Test_Pol_TripleAxis.instr"
    (0)*DEG2RAD);
#line 9762 "Test_Pol_TripleAxis.c"
  rot_mul(mctr1, mcrotaSample_Out, mcrotaAna_Out);
  rot_transpose(mcrotaAna_Arm, mctr1);
  rot_mul(mcrotaAna_Out, mctr1, mcrotrAna_Out);
  mctc1 = coords_set(
#line 165 "Test_Pol_TripleAxis.instr"
    0,
#line 165 "Test_Pol_TripleAxis.instr"
    0,
#line 165 "Test_Pol_TripleAxis.instr"
    0);
#line 9773 "Test_Pol_TripleAxis.c"
  rot_transpose(mcrotaAna_Arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaAna_Out = coords_add(mcposaAna_Arm, mctc2);
  mctc1 = coords_sub(mcposaAna_Arm, mcposaAna_Out);
  mcposrAna_Out = rot_apply(mcrotaAna_Out, mctc1);
  mcDEBUG_COMPONENT("Ana_Out", mcposaAna_Out, mcrotaAna_Out)
  mccomp_posa[15] = mcposaAna_Out;
  mccomp_posr[15] = mcposrAna_Out;
  mcNCounter[15]  = mcPCounter[15] = mcP2Counter[15] = 0;
  mcAbsorbProp[15]= 0;
    /* Component Ana. */
  /* Setting parameters for component Ana. */
  SIG_MESSAGE("Ana (Init:SetPar)");
#line 64 "Test_Pol_TripleAxis.instr"
  mccAna_zmin = -0.05;
#line 64 "Test_Pol_TripleAxis.instr"
  mccAna_zmax = 0.05;
#line 64 "Test_Pol_TripleAxis.instr"
  mccAna_ymin = -0.05;
#line 64 "Test_Pol_TripleAxis.instr"
  mccAna_ymax = 0.05;
#line 168 "Test_Pol_TripleAxis.instr"
  mccAna_zwidth = 0.10;
#line 168 "Test_Pol_TripleAxis.instr"
  mccAna_yheight = 0.08;
#line 169 "Test_Pol_TripleAxis.instr"
  mccAna_mosaich = mcipMOZ;
#line 169 "Test_Pol_TripleAxis.instr"
  mccAna_mosaicv = mcipMOZ;
#line 170 "Test_Pol_TripleAxis.instr"
  mccAna_r0 = 1.0;
#line 170 "Test_Pol_TripleAxis.instr"
  mccAna_Q = mono_q;
#line 66 "Test_Pol_TripleAxis.instr"
  mccAna_DM = 0;
#line 9809 "Test_Pol_TripleAxis.c"

  SIG_MESSAGE("Ana (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9816 "Test_Pol_TripleAxis.c"
  rot_mul(mctr1, mcrotaAna_Arm, mcrotaAna);
  rot_transpose(mcrotaAna_Out, mctr1);
  rot_mul(mcrotaAna, mctr1, mcrotrAna);
  mctc1 = coords_set(
#line 171 "Test_Pol_TripleAxis.instr"
    0,
#line 171 "Test_Pol_TripleAxis.instr"
    0,
#line 171 "Test_Pol_TripleAxis.instr"
    0);
#line 9827 "Test_Pol_TripleAxis.c"
  rot_transpose(mcrotaAna_Arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaAna = coords_add(mcposaAna_Arm, mctc2);
  mctc1 = coords_sub(mcposaAna_Out, mcposaAna);
  mcposrAna = rot_apply(mcrotaAna, mctc1);
  mcDEBUG_COMPONENT("Ana", mcposaAna, mcrotaAna)
  mccomp_posa[16] = mcposaAna;
  mccomp_posr[16] = mcposrAna;
  mcNCounter[16]  = mcPCounter[16] = mcP2Counter[16] = 0;
  mcAbsorbProp[16]= 0;
    /* Component Ana_pol. */
  /* Setting parameters for component Ana_pol. */
  SIG_MESSAGE("Ana_pol (Init:SetPar)");
#line 174 "Test_Pol_TripleAxis.instr"
  mccAna_pol_zwidth = 0.10;
#line 174 "Test_Pol_TripleAxis.instr"
  mccAna_pol_yheight = 0.08;
#line 175 "Test_Pol_TripleAxis.instr"
  mccAna_pol_mosaic = mcipMOZ;
#line 175 "Test_Pol_TripleAxis.instr"
  mccAna_pol_dspread = 0.0;
#line 176 "Test_Pol_TripleAxis.instr"
  mccAna_pol_Q = mono_q;
#line 89 "Test_Pol_TripleAxis.instr"
  mccAna_pol_DM = 0;
#line 89 "Test_Pol_TripleAxis.instr"
  mccAna_pol_pThreshold = 0;
#line 176 "Test_Pol_TripleAxis.instr"
  mccAna_pol_Rup = 1.0;
#line 176 "Test_Pol_TripleAxis.instr"
  mccAna_pol_Rdown = 0.0;
#line 176 "Test_Pol_TripleAxis.instr"
  mccAna_pol_debug = 0;
#line 9861 "Test_Pol_TripleAxis.c"

  SIG_MESSAGE("Ana_pol (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9868 "Test_Pol_TripleAxis.c"
  rot_mul(mctr1, mcrotaAna_Arm, mcrotaAna_pol);
  rot_transpose(mcrotaAna, mctr1);
  rot_mul(mcrotaAna_pol, mctr1, mcrotrAna_pol);
  mctc1 = coords_set(
#line 177 "Test_Pol_TripleAxis.instr"
    0,
#line 177 "Test_Pol_TripleAxis.instr"
    0,
#line 177 "Test_Pol_TripleAxis.instr"
    0);
#line 9879 "Test_Pol_TripleAxis.c"
  rot_transpose(mcrotaAna_Arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaAna_pol = coords_add(mcposaAna_Arm, mctc2);
  mctc1 = coords_sub(mcposaAna, mcposaAna_pol);
  mcposrAna_pol = rot_apply(mcrotaAna_pol, mctc1);
  mcDEBUG_COMPONENT("Ana_pol", mcposaAna_pol, mcrotaAna_pol)
  mccomp_posa[17] = mcposaAna_pol;
  mccomp_posr[17] = mcposrAna_pol;
  mcNCounter[17]  = mcPCounter[17] = mcP2Counter[17] = 0;
  mcAbsorbProp[17]= 0;
    /* Component pollambdaMonitorDet. */
  /* Setting parameters for component pollambdaMonitorDet. */
  SIG_MESSAGE("pollambdaMonitorDet (Init:SetPar)");
#line 183 "Test_Pol_TripleAxis.instr"
  if("pollambdaMonDet.data") strncpy(mccpollambdaMonitorDet_filename, "pollambdaMonDet.data" ? "pollambdaMonDet.data" : "", 16384); else mccpollambdaMonitorDet_filename[0]='\0';
#line 50 "Test_Pol_TripleAxis.instr"
  mccpollambdaMonitorDet_mx = 0;
#line 183 "Test_Pol_TripleAxis.instr"
  mccpollambdaMonitorDet_my = 1;
#line 50 "Test_Pol_TripleAxis.instr"
  mccpollambdaMonitorDet_mz = 0;
#line 182 "Test_Pol_TripleAxis.instr"
  mccpollambdaMonitorDet_Lmin = 0.97 * mcipLAMBDA;
#line 182 "Test_Pol_TripleAxis.instr"
  mccpollambdaMonitorDet_Lmax = 1.03 * mcipLAMBDA;
#line 9905 "Test_Pol_TripleAxis.c"

  SIG_MESSAGE("pollambdaMonitorDet (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9912 "Test_Pol_TripleAxis.c"
  rot_mul(mctr1, mcrotaAna_Out, mcrotapollambdaMonitorDet);
  rot_transpose(mcrotaAna_pol, mctr1);
  rot_mul(mcrotapollambdaMonitorDet, mctr1, mcrotrpollambdaMonitorDet);
  mctc1 = coords_set(
#line 184 "Test_Pol_TripleAxis.instr"
    0,
#line 184 "Test_Pol_TripleAxis.instr"
    0,
#line 184 "Test_Pol_TripleAxis.instr"
    0.15);
#line 9923 "Test_Pol_TripleAxis.c"
  rot_transpose(mcrotaAna_Out, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapollambdaMonitorDet = coords_add(mcposaAna_Out, mctc2);
  mctc1 = coords_sub(mcposaAna_pol, mcposapollambdaMonitorDet);
  mcposrpollambdaMonitorDet = rot_apply(mcrotapollambdaMonitorDet, mctc1);
  mcDEBUG_COMPONENT("pollambdaMonitorDet", mcposapollambdaMonitorDet, mcrotapollambdaMonitorDet)
  mccomp_posa[18] = mcposapollambdaMonitorDet;
  mccomp_posr[18] = mcposrpollambdaMonitorDet;
  mcNCounter[18]  = mcPCounter[18] = mcP2Counter[18] = 0;
  mcAbsorbProp[18]= 0;
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
#line 9960 "Test_Pol_TripleAxis.c"
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
#line 10054 "Test_Pol_TripleAxis.c"
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

  /* Initializations for component Mono_Arm. */
  SIG_MESSAGE("Mono_Arm (Init)");

  /* Initializations for component Mono_Out. */
  SIG_MESSAGE("Mono_Out (Init)");

  /* Initializations for component Mono. */
  SIG_MESSAGE("Mono (Init)");
#define mccompcurname  Mono
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 5
#define mos_rms_y mccMono_mos_rms_y
#define mos_rms_z mccMono_mos_rms_z
#define mos_rms_max mccMono_mos_rms_max
#define mono_Q mccMono_mono_Q
#define zmin mccMono_zmin
#define zmax mccMono_zmax
#define ymin mccMono_ymin
#define ymax mccMono_ymax
#define zwidth mccMono_zwidth
#define yheight mccMono_yheight
#define mosaich mccMono_mosaich
#define mosaicv mccMono_mosaicv
#define r0 mccMono_r0
#define Q mccMono_Q
#define DM mccMono_DM
#line 102 "/usr/share/mcstas/2.5/optics/Monochromator_flat.comp"
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
#line 10116 "Test_Pol_TripleAxis.c"
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

  /* Initializations for component Mono_pol. */
  SIG_MESSAGE("Mono_pol (Init)");
#define mccompcurname  Mono_pol
#define mccompcurtype  Monochromator_pol
#define mccompcurindex 6
#define mos_rms mccMono_pol_mos_rms
#define d_rms mccMono_pol_d_rms
#define mono_Q mccMono_pol_mono_Q
#define FN mccMono_pol_FN
#define FM mccMono_pol_FM
#define zwidth mccMono_pol_zwidth
#define yheight mccMono_pol_yheight
#define mosaic mccMono_pol_mosaic
#define dspread mccMono_pol_dspread
#define Q mccMono_pol_Q
#define DM mccMono_pol_DM
#define pThreshold mccMono_pol_pThreshold
#define Rup mccMono_pol_Rup
#define Rdown mccMono_pol_Rdown
#define debug mccMono_pol_debug
#line 109 "/usr/share/mcstas/2.5/optics/Monochromator_pol.comp"
{
mos_rms = MIN2RAD*mosaic/sqrt(8*log(2));

  mono_Q = Q;
  if (DM != 0)
    mono_Q = 2*PI/DM;

  DM = 2*PI/mono_Q;
  d_rms = dspread*DM/sqrt(8*log(2));

  // calculate the unit cell nuclear and magnetic structure factors
  if(debug > 0)
    printf("Rup: %f, Rdown: %f\n", Rup, Rdown);

  GetMonoPolFNFM(Rup, Rdown, &FN, &FM);

  if(debug > 0)
    printf("FN: %f, FM: %f\n", FN, FM);
}
#line 10176 "Test_Pol_TripleAxis.c"
#undef debug
#undef Rdown
#undef Rup
#undef pThreshold
#undef DM
#undef Q
#undef dspread
#undef mosaic
#undef yheight
#undef zwidth
#undef FM
#undef FN
#undef mono_Q
#undef d_rms
#undef mos_rms
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component pollambdaMonitorMono. */
  SIG_MESSAGE("pollambdaMonitorMono (Init)");
#define mccompcurname  pollambdaMonitorMono
#define mccompcurtype  PolLambda_monitor
#define mccompcurindex 7
#define xwidth mccpollambdaMonitorMono_xwidth
#define yheight mccpollambdaMonitorMono_yheight
#define nL mccpollambdaMonitorMono_nL
#define npol mccpollambdaMonitorMono_npol
#define restore_neutron mccpollambdaMonitorMono_restore_neutron
#define nowritefile mccpollambdaMonitorMono_nowritefile
#define PolL_N mccpollambdaMonitorMono_PolL_N
#define PolL_p mccpollambdaMonitorMono_PolL_p
#define PolL_p2 mccpollambdaMonitorMono_PolL_p2
#define filename mccpollambdaMonitorMono_filename
#define mx mccpollambdaMonitorMono_mx
#define my mccpollambdaMonitorMono_my
#define mz mccpollambdaMonitorMono_mz
#define Lmin mccpollambdaMonitorMono_Lmin
#define Lmax mccpollambdaMonitorMono_Lmax
#line 62 "/usr/share/mcstas/2.5/monitors/PolLambda_monitor.comp"
{
int i, j;

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

  for (i=0; i<npol; i++) {
    for (j=0; j<nL; j++) {

      PolL_N[j][i] = 0;
      PolL_p[j][i] = 0;
      PolL_p2[j][i] = 0;
    }
  }
}
#line 10256 "Test_Pol_TripleAxis.c"
#undef Lmax
#undef Lmin
#undef mz
#undef my
#undef mx
#undef filename
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef npol
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Sample_Arm. */
  SIG_MESSAGE("Sample_Arm (Init)");

  /* Initializations for component Sample_Out. */
  SIG_MESSAGE("Sample_Out (Init)");

  /* Initializations for component Sample. */
  SIG_MESSAGE("Sample (Init)");
#define mccompcurname  Sample
#define mccompcurtype  V_sample
#define mccompcurindex 10
#define VarsV mccSample_VarsV
#define radius mccSample_radius
#define thickness mccSample_thickness
#define zdepth mccSample_zdepth
#define Vc mccSample_Vc
#define sigma_abs mccSample_sigma_abs
#define sigma_inc mccSample_sigma_inc
#define radius_i mccSample_radius_i
#define radius_o mccSample_radius_o
#define h mccSample_h
#define focus_r mccSample_focus_r
#define pack mccSample_pack
#define frac mccSample_frac
#define f_QE mccSample_f_QE
#define gamma mccSample_gamma
#define target_x mccSample_target_x
#define target_y mccSample_target_y
#define target_z mccSample_target_z
#define focus_xw mccSample_focus_xw
#define focus_yh mccSample_focus_yh
#define focus_aw mccSample_focus_aw
#define focus_ah mccSample_focus_ah
#define xwidth mccSample_xwidth
#define yheight mccSample_yheight
#define zthick mccSample_zthick
#define rad_sphere mccSample_rad_sphere
#define sig_a mccSample_sig_a
#define sig_i mccSample_sig_i
#define V0 mccSample_V0
#define target_index mccSample_target_index
#define multiples mccSample_multiples
#line 121 "/usr/share/mcstas/2.5/obsolete/V_sample.comp"
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
#line 10376 "Test_Pol_TripleAxis.c"
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

  /* Initializations for component pollambdaMonitorSample. */
  SIG_MESSAGE("pollambdaMonitorSample (Init)");
#define mccompcurname  pollambdaMonitorSample
#define mccompcurtype  PolLambda_monitor
#define mccompcurindex 11
#define xwidth mccpollambdaMonitorSample_xwidth
#define yheight mccpollambdaMonitorSample_yheight
#define nL mccpollambdaMonitorSample_nL
#define npol mccpollambdaMonitorSample_npol
#define restore_neutron mccpollambdaMonitorSample_restore_neutron
#define nowritefile mccpollambdaMonitorSample_nowritefile
#define PolL_N mccpollambdaMonitorSample_PolL_N
#define PolL_p mccpollambdaMonitorSample_PolL_p
#define PolL_p2 mccpollambdaMonitorSample_PolL_p2
#define filename mccpollambdaMonitorSample_filename
#define mx mccpollambdaMonitorSample_mx
#define my mccpollambdaMonitorSample_my
#define mz mccpollambdaMonitorSample_mz
#define Lmin mccpollambdaMonitorSample_Lmin
#define Lmax mccpollambdaMonitorSample_Lmax
#line 62 "/usr/share/mcstas/2.5/monitors/PolLambda_monitor.comp"
{
int i, j;

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

  for (i=0; i<npol; i++) {
    for (j=0; j<nL; j++) {

      PolL_N[j][i] = 0;
      PolL_p[j][i] = 0;
      PolL_p2[j][i] = 0;
    }
  }
}
#line 10472 "Test_Pol_TripleAxis.c"
#undef Lmax
#undef Lmin
#undef mz
#undef my
#undef mx
#undef filename
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef npol
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component msf. */
  SIG_MESSAGE("msf (Init)");
#define mccompcurname  msf
#define mccompcurtype  Pol_constBfield
#define mccompcurindex 12
#define xwidth mccmsf_xwidth
#define yheight mccmsf_yheight
#define length mccmsf_length
#define omegaL mccmsf_omegaL
#define B mccmsf_B
#define fliplambda mccmsf_fliplambda
#define flipangle mccmsf_flipangle
#line 69 "/usr/share/mcstas/2.5/optics/Pol_constBfield.comp"
{
  omegaL = 0;

  double velocity = 0, time = 0;

  if ((xwidth<=0) || (yheight<=0) || (length<=0)) {
    fprintf(stderr, "Pol_filter: %s: Null or negative volume!\n"
	    "ERROR      (xwidth, yheight, length). Exiting\n",
	    NAME_CURRENT_COMP);
    exit(1);
  }

  // neutron Larmor frequency: omegaL=29.16 MHz/Tesla * B(Neutron Data Booklet)
  // omegaL/B = -2*mu_n/hbar = -2.0*1.913*5.051e-27/1.055e-34
  omegaL = -2 * PI * 29.16e6 * B * 1e-4; // B is in Gauss

  if(fliplambda>0) {

    velocity = K2V*2*PI/fliplambda;
    time = length/velocity;

    // omegaL*time = flipangle
    omegaL = flipangle*DEG2RAD/time;

    // omegaL = 2*PI*29.16 MHz/T*B
    B = omegaL / -2 / 29.16e6 / PI / 1e-4; // Gauss

    printf("Pol_constBfield: %s: Magnetic field set to B=%f Gauss\n",
           NAME_CURRENT_COMP, B);
  }


}
#line 10538 "Test_Pol_TripleAxis.c"
#undef flipangle
#undef fliplambda
#undef B
#undef omegaL
#undef length
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component pollambdaMonitorMSF. */
  SIG_MESSAGE("pollambdaMonitorMSF (Init)");
#define mccompcurname  pollambdaMonitorMSF
#define mccompcurtype  PolLambda_monitor
#define mccompcurindex 13
#define xwidth mccpollambdaMonitorMSF_xwidth
#define yheight mccpollambdaMonitorMSF_yheight
#define nL mccpollambdaMonitorMSF_nL
#define npol mccpollambdaMonitorMSF_npol
#define restore_neutron mccpollambdaMonitorMSF_restore_neutron
#define nowritefile mccpollambdaMonitorMSF_nowritefile
#define PolL_N mccpollambdaMonitorMSF_PolL_N
#define PolL_p mccpollambdaMonitorMSF_PolL_p
#define PolL_p2 mccpollambdaMonitorMSF_PolL_p2
#define filename mccpollambdaMonitorMSF_filename
#define mx mccpollambdaMonitorMSF_mx
#define my mccpollambdaMonitorMSF_my
#define mz mccpollambdaMonitorMSF_mz
#define Lmin mccpollambdaMonitorMSF_Lmin
#define Lmax mccpollambdaMonitorMSF_Lmax
#line 62 "/usr/share/mcstas/2.5/monitors/PolLambda_monitor.comp"
{
int i, j;

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

  for (i=0; i<npol; i++) {
    for (j=0; j<nL; j++) {

      PolL_N[j][i] = 0;
      PolL_p[j][i] = 0;
      PolL_p2[j][i] = 0;
    }
  }
}
#line 10610 "Test_Pol_TripleAxis.c"
#undef Lmax
#undef Lmin
#undef mz
#undef my
#undef mx
#undef filename
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef npol
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Ana_Arm. */
  SIG_MESSAGE("Ana_Arm (Init)");

  /* Initializations for component Ana_Out. */
  SIG_MESSAGE("Ana_Out (Init)");

  /* Initializations for component Ana. */
  SIG_MESSAGE("Ana (Init)");
#define mccompcurname  Ana
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 16
#define mos_rms_y mccAna_mos_rms_y
#define mos_rms_z mccAna_mos_rms_z
#define mos_rms_max mccAna_mos_rms_max
#define mono_Q mccAna_mono_Q
#define zmin mccAna_zmin
#define zmax mccAna_zmax
#define ymin mccAna_ymin
#define ymax mccAna_ymax
#define zwidth mccAna_zwidth
#define yheight mccAna_yheight
#define mosaich mccAna_mosaich
#define mosaicv mccAna_mosaicv
#define r0 mccAna_r0
#define Q mccAna_Q
#define DM mccAna_DM
#line 102 "/usr/share/mcstas/2.5/optics/Monochromator_flat.comp"
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
#line 10671 "Test_Pol_TripleAxis.c"
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

  /* Initializations for component Ana_pol. */
  SIG_MESSAGE("Ana_pol (Init)");
#define mccompcurname  Ana_pol
#define mccompcurtype  Monochromator_pol
#define mccompcurindex 17
#define mos_rms mccAna_pol_mos_rms
#define d_rms mccAna_pol_d_rms
#define mono_Q mccAna_pol_mono_Q
#define FN mccAna_pol_FN
#define FM mccAna_pol_FM
#define zwidth mccAna_pol_zwidth
#define yheight mccAna_pol_yheight
#define mosaic mccAna_pol_mosaic
#define dspread mccAna_pol_dspread
#define Q mccAna_pol_Q
#define DM mccAna_pol_DM
#define pThreshold mccAna_pol_pThreshold
#define Rup mccAna_pol_Rup
#define Rdown mccAna_pol_Rdown
#define debug mccAna_pol_debug
#line 109 "/usr/share/mcstas/2.5/optics/Monochromator_pol.comp"
{
mos_rms = MIN2RAD*mosaic/sqrt(8*log(2));

  mono_Q = Q;
  if (DM != 0)
    mono_Q = 2*PI/DM;

  DM = 2*PI/mono_Q;
  d_rms = dspread*DM/sqrt(8*log(2));

  // calculate the unit cell nuclear and magnetic structure factors
  if(debug > 0)
    printf("Rup: %f, Rdown: %f\n", Rup, Rdown);

  GetMonoPolFNFM(Rup, Rdown, &FN, &FM);

  if(debug > 0)
    printf("FN: %f, FM: %f\n", FN, FM);
}
#line 10731 "Test_Pol_TripleAxis.c"
#undef debug
#undef Rdown
#undef Rup
#undef pThreshold
#undef DM
#undef Q
#undef dspread
#undef mosaic
#undef yheight
#undef zwidth
#undef FM
#undef FN
#undef mono_Q
#undef d_rms
#undef mos_rms
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component pollambdaMonitorDet. */
  SIG_MESSAGE("pollambdaMonitorDet (Init)");
#define mccompcurname  pollambdaMonitorDet
#define mccompcurtype  PolLambda_monitor
#define mccompcurindex 18
#define xwidth mccpollambdaMonitorDet_xwidth
#define yheight mccpollambdaMonitorDet_yheight
#define nL mccpollambdaMonitorDet_nL
#define npol mccpollambdaMonitorDet_npol
#define restore_neutron mccpollambdaMonitorDet_restore_neutron
#define nowritefile mccpollambdaMonitorDet_nowritefile
#define PolL_N mccpollambdaMonitorDet_PolL_N
#define PolL_p mccpollambdaMonitorDet_PolL_p
#define PolL_p2 mccpollambdaMonitorDet_PolL_p2
#define filename mccpollambdaMonitorDet_filename
#define mx mccpollambdaMonitorDet_mx
#define my mccpollambdaMonitorDet_my
#define mz mccpollambdaMonitorDet_mz
#define Lmin mccpollambdaMonitorDet_Lmin
#define Lmax mccpollambdaMonitorDet_Lmax
#line 62 "/usr/share/mcstas/2.5/monitors/PolLambda_monitor.comp"
{
int i, j;

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

  for (i=0; i<npol; i++) {
    for (j=0; j<nL; j++) {

      PolL_N[j][i] = 0;
      PolL_p[j][i] = 0;
      PolL_p2[j][i] = 0;
    }
  }
}
#line 10811 "Test_Pol_TripleAxis.c"
#undef Lmax
#undef Lmin
#undef mz
#undef my
#undef mx
#undef filename
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef npol
#undef nL
#undef yheight
#undef xwidth
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
#line 10984 "Test_Pol_TripleAxis.c"
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
#line 11155 "Test_Pol_TripleAxis.c"
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

  /* TRACE Component Mono_Arm [3] */
  mccoordschange(mcposrMono_Arm, mcrotrMono_Arm,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Mono_Arm (without coords transformations) */
  mcJumpTrace_Mono_Arm:
  SIG_MESSAGE("Mono_Arm (Trace)");
  mcDEBUG_COMP("Mono_Arm")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompMono_Arm
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
#define mccompcurname  Mono_Arm
#define mccompcurtype  Arm
#define mccompcurindex 3
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompMono_Arm:
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

  /* TRACE Component Mono_Out [4] */
  mccoordschange(mcposrMono_Out, mcrotrMono_Out,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Mono_Out (without coords transformations) */
  mcJumpTrace_Mono_Out:
  SIG_MESSAGE("Mono_Out (Trace)");
  mcDEBUG_COMP("Mono_Out")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompMono_Out
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
#define mccompcurname  Mono_Out
#define mccompcurtype  Arm
#define mccompcurindex 4
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompMono_Out:
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

  /* TRACE Component Mono [5] */
  mccoordschange(mcposrMono, mcrotrMono,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Mono (without coords transformations) */
  mcJumpTrace_Mono:
  SIG_MESSAGE("Mono (Trace)");
  mcDEBUG_COMP("Mono")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompMono
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
#define mccompcurname  Mono
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 5
#define mos_rms_y mccMono_mos_rms_y
#define mos_rms_z mccMono_mos_rms_z
#define mos_rms_max mccMono_mos_rms_max
#define mono_Q mccMono_mono_Q
{   /* Declarations of Mono=Monochromator_flat() SETTING parameters. */
MCNUM zmin = mccMono_zmin;
MCNUM zmax = mccMono_zmax;
MCNUM ymin = mccMono_ymin;
MCNUM ymax = mccMono_ymax;
MCNUM zwidth = mccMono_zwidth;
MCNUM yheight = mccMono_yheight;
MCNUM mosaich = mccMono_mosaich;
MCNUM mosaicv = mccMono_mosaicv;
MCNUM r0 = mccMono_r0;
MCNUM Q = mccMono_Q;
MCNUM DM = mccMono_DM;
/* 'Mono=Monochromator_flat()' component instance has conditional execution */
if (( mcipOPTION == 0 ))

#line 118 "/usr/share/mcstas/2.5/optics/Monochromator_flat.comp"
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
#line 11622 "Test_Pol_TripleAxis.c"
}   /* End of Mono=Monochromator_flat() SETTING parameter declarations. */
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompMono:
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

  /* TRACE Component Mono_pol [6] */
  mccoordschange(mcposrMono_pol, mcrotrMono_pol,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Mono_pol (without coords transformations) */
  mcJumpTrace_Mono_pol:
  SIG_MESSAGE("Mono_pol (Trace)");
  mcDEBUG_COMP("Mono_pol")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompMono_pol
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
#define mccompcurname  Mono_pol
#define mccompcurtype  Monochromator_pol
#define mccompcurindex 6
#define mos_rms mccMono_pol_mos_rms
#define d_rms mccMono_pol_d_rms
#define mono_Q mccMono_pol_mono_Q
#define FN mccMono_pol_FN
#define FM mccMono_pol_FM
{   /* Declarations of Mono_pol=Monochromator_pol() SETTING parameters. */
MCNUM zwidth = mccMono_pol_zwidth;
MCNUM yheight = mccMono_pol_yheight;
MCNUM mosaic = mccMono_pol_mosaic;
MCNUM dspread = mccMono_pol_dspread;
MCNUM Q = mccMono_pol_Q;
MCNUM DM = mccMono_pol_DM;
MCNUM pThreshold = mccMono_pol_pThreshold;
MCNUM Rup = mccMono_pol_Rup;
MCNUM Rdown = mccMono_pol_Rdown;
int debug = mccMono_pol_debug;
/* 'Mono_pol=Monochromator_pol()' component instance has conditional execution */
if (( mcipOPTION == -1 || mcipOPTION == 1 ))

#line 130 "/usr/share/mcstas/2.5/optics/Monochromator_pol.comp"
{
  double y1, z1, t1, dt, vel;
  double sinTheta, lambdaBragg, lambda, dlambda2, sigmaLambda2, p_reflect;
  double R0; /* reflection probability based on FN and FM */
  double sx_in, sy_in, sz_in;
  int i;

  /* Propagate to crystal */
  PROP_X0;

  if (inside_rectangle(z, y, zwidth, yheight)) {/* Intersect the crystal? */

    // calculate sin(Bragg angle)
    vel = sqrt(vx*vx + vy*vy + vz*vz);
    sinTheta = abs(vx)/vel;

    // calculate lambdaBragg
    lambdaBragg = 2.0*DM*sinTheta;

    // calculate lambda of neutron
    lambda = 2*PI/(V2K*vel);


    // calculate deltalambda squared and sigmaLambda squared
    dlambda2 = (lambda-lambdaBragg)*(lambda-lambdaBragg);
    // The sigmaLambda is propagated by differentiating the Bragg
    // condition: lambda = 2*d*sinTheta
    sigmaLambda2 = 2.0*2.0 * sinTheta*sinTheta * d_rms*d_rms+
      2.0*2.0 * DM*DM * (1.0-sinTheta*sinTheta) * mos_rms*mos_rms;

    // calculate peak reflection probability
    GetMonoPolRefProb(FN, FM, sy, &R0);

    // calculate reflection probability
    p_reflect = R0*exp(-dlambda2/(2.0*sigmaLambda2));

    if(debug > 0) {
      printf("\n lambda: %f, Lambda_Bragg: %f\n", lambda, lambdaBragg);
      printf("sigmaLambda: %f, R0: %f, p_reflect: %f\n",
	     sqrt(sigmaLambda2), R0, p_reflect);
      printf("S_in:  (%f, %f, %f)\n", sx, sy, sz);
    }

    if((pThreshold>0 && p_reflect>pThreshold) || rand01()<p_reflect) {
      /* Reflect */

      // scale weight if neutron was accepted because of threshold
      if(pThreshold>0 && p_reflect>pThreshold)
	p*=p_reflect;

      vx = -vx;

      // Outgoing polarisation
      SetMonoPolRefOut(FN, FM, R0, &sx, &sy, &sz);

      if(debug > 0)
	printf("S_out: (%f, %f, %f)\n", sx, sy, sz);

      if(sx*sx+sy*sy+sz*sz>1)
        fprintf(stderr,"Pol_mirror: %s: Warning: polarisation |s| = %g > 1\n",
	      NAME_CURRENT_COMP, sx*sx+sy*sy+sz*sz); // check that polarisation is meaningfull

      SCATTER;
    } /* End MC choice to reflect or transmit neutron */
  } /* End intersect the crystal */

}
#line 11817 "Test_Pol_TripleAxis.c"
}   /* End of Mono_pol=Monochromator_pol() SETTING parameter declarations. */
#undef FM
#undef FN
#undef mono_Q
#undef d_rms
#undef mos_rms
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompMono_pol:
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

  /* TRACE Component pollambdaMonitorMono [7] */
  mccoordschange(mcposrpollambdaMonitorMono, mcrotrpollambdaMonitorMono,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component pollambdaMonitorMono (without coords transformations) */
  mcJumpTrace_pollambdaMonitorMono:
  SIG_MESSAGE("pollambdaMonitorMono (Trace)");
  mcDEBUG_COMP("pollambdaMonitorMono")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComppollambdaMonitorMono
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
#define mccompcurname  pollambdaMonitorMono
#define mccompcurtype  PolLambda_monitor
#define mccompcurindex 7
#define xwidth mccpollambdaMonitorMono_xwidth
#define yheight mccpollambdaMonitorMono_yheight
#define nL mccpollambdaMonitorMono_nL
#define npol mccpollambdaMonitorMono_npol
#define restore_neutron mccpollambdaMonitorMono_restore_neutron
#define nowritefile mccpollambdaMonitorMono_nowritefile
#define PolL_N mccpollambdaMonitorMono_PolL_N
#define PolL_p mccpollambdaMonitorMono_PolL_p
#define PolL_p2 mccpollambdaMonitorMono_PolL_p2
{   /* Declarations of pollambdaMonitorMono=PolLambda_monitor() SETTING parameters. */
char* filename = mccpollambdaMonitorMono_filename;
MCNUM mx = mccpollambdaMonitorMono_mx;
MCNUM my = mccpollambdaMonitorMono_my;
MCNUM mz = mccpollambdaMonitorMono_mz;
MCNUM Lmin = mccpollambdaMonitorMono_Lmin;
MCNUM Lmax = mccpollambdaMonitorMono_Lmax;
#line 103 "/usr/share/mcstas/2.5/monitors/PolLambda_monitor.comp"
{
  int i, j;
  double pol_proj;
  double lambda;

  PROP_Z0;
  lambda = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);

  if (inside_rectangle(x, y, xwidth, yheight) &&
      lambda > Lmin && lambda < Lmax) {

    pol_proj = scalar_prod(mx, my, mz, sx, sy, sz);
    /*protect from rounding errors introduced by trig functions*/
    if (fabs(pol_proj-1)<FLT_EPSILON)
      i = npol-1;
    else if (fabs(pol_proj+1.0)<FLT_EPSILON)
      i=0;
    else
      i = floor((pol_proj+1.0)*npol/2.0);
    j = floor((lambda - Lmin)*nL/(Lmax - Lmin));

    PolL_N[j][i]++;
    PolL_p[j][i] += p;
    PolL_p2[j][i] += p*p;
    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 11974 "Test_Pol_TripleAxis.c"
}   /* End of pollambdaMonitorMono=PolLambda_monitor() SETTING parameter declarations. */
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef npol
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppollambdaMonitorMono:
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

  /* TRACE Component Sample_Arm [8] */
  mccoordschange(mcposrSample_Arm, mcrotrSample_Arm,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Sample_Arm (without coords transformations) */
  mcJumpTrace_Sample_Arm:
  SIG_MESSAGE("Sample_Arm (Trace)");
  mcDEBUG_COMP("Sample_Arm")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompSample_Arm
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
#define mccompcurname  Sample_Arm
#define mccompcurtype  Arm
#define mccompcurindex 8
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompSample_Arm:
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

  /* TRACE Component Sample_Out [9] */
  mccoordschange(mcposrSample_Out, mcrotrSample_Out,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Sample_Out (without coords transformations) */
  mcJumpTrace_Sample_Out:
  SIG_MESSAGE("Sample_Out (Trace)");
  mcDEBUG_COMP("Sample_Out")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompSample_Out
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
#define mccompcurname  Sample_Out
#define mccompcurtype  Arm
#define mccompcurindex 9
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompSample_Out:
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

  /* TRACE Component Sample [10] */
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
#define mccompcurname  Sample
#define mccompcurtype  V_sample
#define mccompcurindex 10
#define VarsV mccSample_VarsV
{   /* Declarations of Sample=V_sample() SETTING parameters. */
MCNUM radius = mccSample_radius;
MCNUM thickness = mccSample_thickness;
MCNUM zdepth = mccSample_zdepth;
MCNUM Vc = mccSample_Vc;
MCNUM sigma_abs = mccSample_sigma_abs;
MCNUM sigma_inc = mccSample_sigma_inc;
MCNUM radius_i = mccSample_radius_i;
MCNUM radius_o = mccSample_radius_o;
MCNUM h = mccSample_h;
MCNUM focus_r = mccSample_focus_r;
MCNUM pack = mccSample_pack;
MCNUM frac = mccSample_frac;
MCNUM f_QE = mccSample_f_QE;
MCNUM gamma = mccSample_gamma;
MCNUM target_x = mccSample_target_x;
MCNUM target_y = mccSample_target_y;
MCNUM target_z = mccSample_target_z;
MCNUM focus_xw = mccSample_focus_xw;
MCNUM focus_yh = mccSample_focus_yh;
MCNUM focus_aw = mccSample_focus_aw;
MCNUM focus_ah = mccSample_focus_ah;
MCNUM xwidth = mccSample_xwidth;
MCNUM yheight = mccSample_yheight;
MCNUM zthick = mccSample_zthick;
MCNUM rad_sphere = mccSample_rad_sphere;
MCNUM sig_a = mccSample_sig_a;
MCNUM sig_i = mccSample_sig_i;
MCNUM V0 = mccSample_V0;
int target_index = mccSample_target_index;
MCNUM multiples = mccSample_multiples;
#line 180 "/usr/share/mcstas/2.5/obsolete/V_sample.comp"
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
#line 12465 "Test_Pol_TripleAxis.c"
}   /* End of Sample=V_sample() SETTING parameter declarations. */
#undef VarsV
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompSample:
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

  /* TRACE Component pollambdaMonitorSample [11] */
  mccoordschange(mcposrpollambdaMonitorSample, mcrotrpollambdaMonitorSample,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component pollambdaMonitorSample (without coords transformations) */
  mcJumpTrace_pollambdaMonitorSample:
  SIG_MESSAGE("pollambdaMonitorSample (Trace)");
  mcDEBUG_COMP("pollambdaMonitorSample")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComppollambdaMonitorSample
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
#define mccompcurname  pollambdaMonitorSample
#define mccompcurtype  PolLambda_monitor
#define mccompcurindex 11
#define xwidth mccpollambdaMonitorSample_xwidth
#define yheight mccpollambdaMonitorSample_yheight
#define nL mccpollambdaMonitorSample_nL
#define npol mccpollambdaMonitorSample_npol
#define restore_neutron mccpollambdaMonitorSample_restore_neutron
#define nowritefile mccpollambdaMonitorSample_nowritefile
#define PolL_N mccpollambdaMonitorSample_PolL_N
#define PolL_p mccpollambdaMonitorSample_PolL_p
#define PolL_p2 mccpollambdaMonitorSample_PolL_p2
{   /* Declarations of pollambdaMonitorSample=PolLambda_monitor() SETTING parameters. */
char* filename = mccpollambdaMonitorSample_filename;
MCNUM mx = mccpollambdaMonitorSample_mx;
MCNUM my = mccpollambdaMonitorSample_my;
MCNUM mz = mccpollambdaMonitorSample_mz;
MCNUM Lmin = mccpollambdaMonitorSample_Lmin;
MCNUM Lmax = mccpollambdaMonitorSample_Lmax;
#line 103 "/usr/share/mcstas/2.5/monitors/PolLambda_monitor.comp"
{
  int i, j;
  double pol_proj;
  double lambda;

  PROP_Z0;
  lambda = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);

  if (inside_rectangle(x, y, xwidth, yheight) &&
      lambda > Lmin && lambda < Lmax) {

    pol_proj = scalar_prod(mx, my, mz, sx, sy, sz);
    /*protect from rounding errors introduced by trig functions*/
    if (fabs(pol_proj-1)<FLT_EPSILON)
      i = npol-1;
    else if (fabs(pol_proj+1.0)<FLT_EPSILON)
      i=0;
    else
      i = floor((pol_proj+1.0)*npol/2.0);
    j = floor((lambda - Lmin)*nL/(Lmax - Lmin));

    PolL_N[j][i]++;
    PolL_p[j][i] += p;
    PolL_p2[j][i] += p*p;
    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 12618 "Test_Pol_TripleAxis.c"
}   /* End of pollambdaMonitorSample=PolLambda_monitor() SETTING parameter declarations. */
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef npol
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppollambdaMonitorSample:
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

  /* TRACE Component msf [12] */
  mccoordschange(mcposrmsf, mcrotrmsf,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component msf (without coords transformations) */
  mcJumpTrace_msf:
  SIG_MESSAGE("msf (Trace)");
  mcDEBUG_COMP("msf")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompmsf
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
#define mccompcurname  msf
#define mccompcurtype  Pol_constBfield
#define mccompcurindex 12
#define xwidth mccmsf_xwidth
#define yheight mccmsf_yheight
#define length mccmsf_length
#define omegaL mccmsf_omegaL
{   /* Declarations of msf=Pol_constBfield() SETTING parameters. */
MCNUM B = mccmsf_B;
MCNUM fliplambda = mccmsf_fliplambda;
MCNUM flipangle = mccmsf_flipangle;
/* 'msf=Pol_constBfield()' component instance has conditional execution */
if (( mcipOPTION == -1 ))

#line 104 "/usr/share/mcstas/2.5/optics/Pol_constBfield.comp"
{
  double deltaT, deltaTx, deltaTy, sx_in, sz_in;

  PROP_Z0;
  if (!inside_rectangle(x, y, xwidth, yheight))
    ABSORB;

  // Time spent in B-field
  deltaT = length/vz;

  // check that track goes throgh without hitting the walls
  if (!inside_rectangle(x+vx*deltaT, y+vy*deltaT, xwidth, yheight)) {

    // Propagate to the wall and absorb
    deltaTx = IntersectWall(x, vx, xwidth/2);
    deltaTy = IntersectWall(y, vy, yheight/2);

    if (deltaTx>=0 && deltaTx<deltaTy)
      deltaT = deltaTx;
    else
      deltaT = deltaTy;

    PROP_DT(deltaT);

    ABSORB;
  }

  PROP_DT(deltaT);

  // The solution to neutron precession in a constant field along the
  // z-axis is:
  // sx(deltaT) = cos(omegaL*deltaT)*sx(0) - sin(omegaL*deltaT)*sy(0)
  // sy(deltaT) = sin(omegaL*deltaT)*sx(0) + cos(omegaL*deltaT)*sy(0)
  // sz(deltaT) = sz(0)
  // see e.g. W. Gavin Willams "Polarized Neutrons", page 18.
  // In our geometry (x', y', z') where y' = z we get x' = y and z'= x

  sx_in = sx;
  sz_in = sz;

  sz = cos(omegaL*deltaT)*sz_in - sin(omegaL*deltaT)*sx_in;
  sx = sin(omegaL*deltaT)*sz_in + cos(omegaL*deltaT)*sx_in;
}
#line 12786 "Test_Pol_TripleAxis.c"
}   /* End of msf=Pol_constBfield() SETTING parameter declarations. */
#undef omegaL
#undef length
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompmsf:
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

  /* TRACE Component pollambdaMonitorMSF [13] */
  mccoordschange(mcposrpollambdaMonitorMSF, mcrotrpollambdaMonitorMSF,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component pollambdaMonitorMSF (without coords transformations) */
  mcJumpTrace_pollambdaMonitorMSF:
  SIG_MESSAGE("pollambdaMonitorMSF (Trace)");
  mcDEBUG_COMP("pollambdaMonitorMSF")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComppollambdaMonitorMSF
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
#define mccompcurname  pollambdaMonitorMSF
#define mccompcurtype  PolLambda_monitor
#define mccompcurindex 13
#define xwidth mccpollambdaMonitorMSF_xwidth
#define yheight mccpollambdaMonitorMSF_yheight
#define nL mccpollambdaMonitorMSF_nL
#define npol mccpollambdaMonitorMSF_npol
#define restore_neutron mccpollambdaMonitorMSF_restore_neutron
#define nowritefile mccpollambdaMonitorMSF_nowritefile
#define PolL_N mccpollambdaMonitorMSF_PolL_N
#define PolL_p mccpollambdaMonitorMSF_PolL_p
#define PolL_p2 mccpollambdaMonitorMSF_PolL_p2
{   /* Declarations of pollambdaMonitorMSF=PolLambda_monitor() SETTING parameters. */
char* filename = mccpollambdaMonitorMSF_filename;
MCNUM mx = mccpollambdaMonitorMSF_mx;
MCNUM my = mccpollambdaMonitorMSF_my;
MCNUM mz = mccpollambdaMonitorMSF_mz;
MCNUM Lmin = mccpollambdaMonitorMSF_Lmin;
MCNUM Lmax = mccpollambdaMonitorMSF_Lmax;
#line 103 "/usr/share/mcstas/2.5/monitors/PolLambda_monitor.comp"
{
  int i, j;
  double pol_proj;
  double lambda;

  PROP_Z0;
  lambda = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);

  if (inside_rectangle(x, y, xwidth, yheight) &&
      lambda > Lmin && lambda < Lmax) {

    pol_proj = scalar_prod(mx, my, mz, sx, sy, sz);
    /*protect from rounding errors introduced by trig functions*/
    if (fabs(pol_proj-1)<FLT_EPSILON)
      i = npol-1;
    else if (fabs(pol_proj+1.0)<FLT_EPSILON)
      i=0;
    else
      i = floor((pol_proj+1.0)*npol/2.0);
    j = floor((lambda - Lmin)*nL/(Lmax - Lmin));

    PolL_N[j][i]++;
    PolL_p[j][i] += p;
    PolL_p2[j][i] += p*p;
    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 12942 "Test_Pol_TripleAxis.c"
}   /* End of pollambdaMonitorMSF=PolLambda_monitor() SETTING parameter declarations. */
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef npol
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppollambdaMonitorMSF:
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

  /* TRACE Component Ana_Arm [14] */
  mccoordschange(mcposrAna_Arm, mcrotrAna_Arm,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Ana_Arm (without coords transformations) */
  mcJumpTrace_Ana_Arm:
  SIG_MESSAGE("Ana_Arm (Trace)");
  mcDEBUG_COMP("Ana_Arm")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompAna_Arm
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
#define mccompcurname  Ana_Arm
#define mccompcurtype  Arm
#define mccompcurindex 14
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompAna_Arm:
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

  /* TRACE Component Ana_Out [15] */
  mccoordschange(mcposrAna_Out, mcrotrAna_Out,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Ana_Out (without coords transformations) */
  mcJumpTrace_Ana_Out:
  SIG_MESSAGE("Ana_Out (Trace)");
  mcDEBUG_COMP("Ana_Out")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompAna_Out
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
#define mccompcurname  Ana_Out
#define mccompcurtype  Arm
#define mccompcurindex 15
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompAna_Out:
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

  /* TRACE Component Ana [16] */
  mccoordschange(mcposrAna, mcrotrAna,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Ana (without coords transformations) */
  mcJumpTrace_Ana:
  SIG_MESSAGE("Ana (Trace)");
  mcDEBUG_COMP("Ana")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompAna
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
#define mccompcurname  Ana
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 16
#define mos_rms_y mccAna_mos_rms_y
#define mos_rms_z mccAna_mos_rms_z
#define mos_rms_max mccAna_mos_rms_max
#define mono_Q mccAna_mono_Q
{   /* Declarations of Ana=Monochromator_flat() SETTING parameters. */
MCNUM zmin = mccAna_zmin;
MCNUM zmax = mccAna_zmax;
MCNUM ymin = mccAna_ymin;
MCNUM ymax = mccAna_ymax;
MCNUM zwidth = mccAna_zwidth;
MCNUM yheight = mccAna_yheight;
MCNUM mosaich = mccAna_mosaich;
MCNUM mosaicv = mccAna_mosaicv;
MCNUM r0 = mccAna_r0;
MCNUM Q = mccAna_Q;
MCNUM DM = mccAna_DM;
/* 'Ana=Monochromator_flat()' component instance has conditional execution */
if (( mcipOPTION == 0 ))

#line 118 "/usr/share/mcstas/2.5/optics/Monochromator_flat.comp"
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
#line 13415 "Test_Pol_TripleAxis.c"
}   /* End of Ana=Monochromator_flat() SETTING parameter declarations. */
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompAna:
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

  /* TRACE Component Ana_pol [17] */
  mccoordschange(mcposrAna_pol, mcrotrAna_pol,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Ana_pol (without coords transformations) */
  mcJumpTrace_Ana_pol:
  SIG_MESSAGE("Ana_pol (Trace)");
  mcDEBUG_COMP("Ana_pol")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompAna_pol
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
#define mccompcurname  Ana_pol
#define mccompcurtype  Monochromator_pol
#define mccompcurindex 17
#define mos_rms mccAna_pol_mos_rms
#define d_rms mccAna_pol_d_rms
#define mono_Q mccAna_pol_mono_Q
#define FN mccAna_pol_FN
#define FM mccAna_pol_FM
{   /* Declarations of Ana_pol=Monochromator_pol() SETTING parameters. */
MCNUM zwidth = mccAna_pol_zwidth;
MCNUM yheight = mccAna_pol_yheight;
MCNUM mosaic = mccAna_pol_mosaic;
MCNUM dspread = mccAna_pol_dspread;
MCNUM Q = mccAna_pol_Q;
MCNUM DM = mccAna_pol_DM;
MCNUM pThreshold = mccAna_pol_pThreshold;
MCNUM Rup = mccAna_pol_Rup;
MCNUM Rdown = mccAna_pol_Rdown;
int debug = mccAna_pol_debug;
/* 'Ana_pol=Monochromator_pol()' component instance has conditional execution */
if (( mcipOPTION == -1 || mcipOPTION == 1 ))

#line 130 "/usr/share/mcstas/2.5/optics/Monochromator_pol.comp"
{
  double y1, z1, t1, dt, vel;
  double sinTheta, lambdaBragg, lambda, dlambda2, sigmaLambda2, p_reflect;
  double R0; /* reflection probability based on FN and FM */
  double sx_in, sy_in, sz_in;
  int i;

  /* Propagate to crystal */
  PROP_X0;

  if (inside_rectangle(z, y, zwidth, yheight)) {/* Intersect the crystal? */

    // calculate sin(Bragg angle)
    vel = sqrt(vx*vx + vy*vy + vz*vz);
    sinTheta = abs(vx)/vel;

    // calculate lambdaBragg
    lambdaBragg = 2.0*DM*sinTheta;

    // calculate lambda of neutron
    lambda = 2*PI/(V2K*vel);


    // calculate deltalambda squared and sigmaLambda squared
    dlambda2 = (lambda-lambdaBragg)*(lambda-lambdaBragg);
    // The sigmaLambda is propagated by differentiating the Bragg
    // condition: lambda = 2*d*sinTheta
    sigmaLambda2 = 2.0*2.0 * sinTheta*sinTheta * d_rms*d_rms+
      2.0*2.0 * DM*DM * (1.0-sinTheta*sinTheta) * mos_rms*mos_rms;

    // calculate peak reflection probability
    GetMonoPolRefProb(FN, FM, sy, &R0);

    // calculate reflection probability
    p_reflect = R0*exp(-dlambda2/(2.0*sigmaLambda2));

    if(debug > 0) {
      printf("\n lambda: %f, Lambda_Bragg: %f\n", lambda, lambdaBragg);
      printf("sigmaLambda: %f, R0: %f, p_reflect: %f\n",
	     sqrt(sigmaLambda2), R0, p_reflect);
      printf("S_in:  (%f, %f, %f)\n", sx, sy, sz);
    }

    if((pThreshold>0 && p_reflect>pThreshold) || rand01()<p_reflect) {
      /* Reflect */

      // scale weight if neutron was accepted because of threshold
      if(pThreshold>0 && p_reflect>pThreshold)
	p*=p_reflect;

      vx = -vx;

      // Outgoing polarisation
      SetMonoPolRefOut(FN, FM, R0, &sx, &sy, &sz);

      if(debug > 0)
	printf("S_out: (%f, %f, %f)\n", sx, sy, sz);

      if(sx*sx+sy*sy+sz*sz>1)
        fprintf(stderr,"Pol_mirror: %s: Warning: polarisation |s| = %g > 1\n",
	      NAME_CURRENT_COMP, sx*sx+sy*sy+sz*sz); // check that polarisation is meaningfull

      SCATTER;
    } /* End MC choice to reflect or transmit neutron */
  } /* End intersect the crystal */

}
#line 13610 "Test_Pol_TripleAxis.c"
}   /* End of Ana_pol=Monochromator_pol() SETTING parameter declarations. */
#undef FM
#undef FN
#undef mono_Q
#undef d_rms
#undef mos_rms
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompAna_pol:
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

  /* TRACE Component pollambdaMonitorDet [18] */
  mccoordschange(mcposrpollambdaMonitorDet, mcrotrpollambdaMonitorDet,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component pollambdaMonitorDet (without coords transformations) */
  mcJumpTrace_pollambdaMonitorDet:
  SIG_MESSAGE("pollambdaMonitorDet (Trace)");
  mcDEBUG_COMP("pollambdaMonitorDet")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComppollambdaMonitorDet
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
#define mccompcurname  pollambdaMonitorDet
#define mccompcurtype  PolLambda_monitor
#define mccompcurindex 18
#define xwidth mccpollambdaMonitorDet_xwidth
#define yheight mccpollambdaMonitorDet_yheight
#define nL mccpollambdaMonitorDet_nL
#define npol mccpollambdaMonitorDet_npol
#define restore_neutron mccpollambdaMonitorDet_restore_neutron
#define nowritefile mccpollambdaMonitorDet_nowritefile
#define PolL_N mccpollambdaMonitorDet_PolL_N
#define PolL_p mccpollambdaMonitorDet_PolL_p
#define PolL_p2 mccpollambdaMonitorDet_PolL_p2
{   /* Declarations of pollambdaMonitorDet=PolLambda_monitor() SETTING parameters. */
char* filename = mccpollambdaMonitorDet_filename;
MCNUM mx = mccpollambdaMonitorDet_mx;
MCNUM my = mccpollambdaMonitorDet_my;
MCNUM mz = mccpollambdaMonitorDet_mz;
MCNUM Lmin = mccpollambdaMonitorDet_Lmin;
MCNUM Lmax = mccpollambdaMonitorDet_Lmax;
#line 103 "/usr/share/mcstas/2.5/monitors/PolLambda_monitor.comp"
{
  int i, j;
  double pol_proj;
  double lambda;

  PROP_Z0;
  lambda = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);

  if (inside_rectangle(x, y, xwidth, yheight) &&
      lambda > Lmin && lambda < Lmax) {

    pol_proj = scalar_prod(mx, my, mz, sx, sy, sz);
    /*protect from rounding errors introduced by trig functions*/
    if (fabs(pol_proj-1)<FLT_EPSILON)
      i = npol-1;
    else if (fabs(pol_proj+1.0)<FLT_EPSILON)
      i=0;
    else
      i = floor((pol_proj+1.0)*npol/2.0);
    j = floor((lambda - Lmin)*nL/(Lmax - Lmin));

    PolL_N[j][i]++;
    PolL_p[j][i] += p;
    PolL_p2[j][i] += p*p;
    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 13767 "Test_Pol_TripleAxis.c"
}   /* End of pollambdaMonitorDet=PolLambda_monitor() SETTING parameter declarations. */
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef npol
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppollambdaMonitorDet:
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
#line 13885 "Test_Pol_TripleAxis.c"
}   /* End of Origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'pollambdaMonitorMono'. */
  SIG_MESSAGE("pollambdaMonitorMono (Save)");
#define mccompcurname  pollambdaMonitorMono
#define mccompcurtype  PolLambda_monitor
#define mccompcurindex 7
#define xwidth mccpollambdaMonitorMono_xwidth
#define yheight mccpollambdaMonitorMono_yheight
#define nL mccpollambdaMonitorMono_nL
#define npol mccpollambdaMonitorMono_npol
#define restore_neutron mccpollambdaMonitorMono_restore_neutron
#define nowritefile mccpollambdaMonitorMono_nowritefile
#define PolL_N mccpollambdaMonitorMono_PolL_N
#define PolL_p mccpollambdaMonitorMono_PolL_p
#define PolL_p2 mccpollambdaMonitorMono_PolL_p2
{   /* Declarations of pollambdaMonitorMono=PolLambda_monitor() SETTING parameters. */
char* filename = mccpollambdaMonitorMono_filename;
MCNUM mx = mccpollambdaMonitorMono_mx;
MCNUM my = mccpollambdaMonitorMono_my;
MCNUM mz = mccpollambdaMonitorMono_mz;
MCNUM Lmin = mccpollambdaMonitorMono_Lmin;
MCNUM Lmax = mccpollambdaMonitorMono_Lmax;
#line 135 "/usr/share/mcstas/2.5/monitors/PolLambda_monitor.comp"
{
    if (!nowritefile) {
  DETECTOR_OUT_2D("Pol-wavelength monitor",
		  "Wavelength [AA]", "Polarisation projection",
		  Lmin, Lmax, -1.0, 1.0,
		  nL, npol,
		  &PolL_N[0][0],&PolL_p[0][0],&PolL_p2[0][0],
		  filename);
    }
}
#line 13927 "Test_Pol_TripleAxis.c"
}   /* End of pollambdaMonitorMono=PolLambda_monitor() SETTING parameter declarations. */
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef npol
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'pollambdaMonitorSample'. */
  SIG_MESSAGE("pollambdaMonitorSample (Save)");
#define mccompcurname  pollambdaMonitorSample
#define mccompcurtype  PolLambda_monitor
#define mccompcurindex 11
#define xwidth mccpollambdaMonitorSample_xwidth
#define yheight mccpollambdaMonitorSample_yheight
#define nL mccpollambdaMonitorSample_nL
#define npol mccpollambdaMonitorSample_npol
#define restore_neutron mccpollambdaMonitorSample_restore_neutron
#define nowritefile mccpollambdaMonitorSample_nowritefile
#define PolL_N mccpollambdaMonitorSample_PolL_N
#define PolL_p mccpollambdaMonitorSample_PolL_p
#define PolL_p2 mccpollambdaMonitorSample_PolL_p2
{   /* Declarations of pollambdaMonitorSample=PolLambda_monitor() SETTING parameters. */
char* filename = mccpollambdaMonitorSample_filename;
MCNUM mx = mccpollambdaMonitorSample_mx;
MCNUM my = mccpollambdaMonitorSample_my;
MCNUM mz = mccpollambdaMonitorSample_mz;
MCNUM Lmin = mccpollambdaMonitorSample_Lmin;
MCNUM Lmax = mccpollambdaMonitorSample_Lmax;
#line 135 "/usr/share/mcstas/2.5/monitors/PolLambda_monitor.comp"
{
    if (!nowritefile) {
  DETECTOR_OUT_2D("Pol-wavelength monitor",
		  "Wavelength [AA]", "Polarisation projection",
		  Lmin, Lmax, -1.0, 1.0,
		  nL, npol,
		  &PolL_N[0][0],&PolL_p[0][0],&PolL_p2[0][0],
		  filename);
    }
}
#line 13974 "Test_Pol_TripleAxis.c"
}   /* End of pollambdaMonitorSample=PolLambda_monitor() SETTING parameter declarations. */
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef npol
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'pollambdaMonitorMSF'. */
  SIG_MESSAGE("pollambdaMonitorMSF (Save)");
#define mccompcurname  pollambdaMonitorMSF
#define mccompcurtype  PolLambda_monitor
#define mccompcurindex 13
#define xwidth mccpollambdaMonitorMSF_xwidth
#define yheight mccpollambdaMonitorMSF_yheight
#define nL mccpollambdaMonitorMSF_nL
#define npol mccpollambdaMonitorMSF_npol
#define restore_neutron mccpollambdaMonitorMSF_restore_neutron
#define nowritefile mccpollambdaMonitorMSF_nowritefile
#define PolL_N mccpollambdaMonitorMSF_PolL_N
#define PolL_p mccpollambdaMonitorMSF_PolL_p
#define PolL_p2 mccpollambdaMonitorMSF_PolL_p2
{   /* Declarations of pollambdaMonitorMSF=PolLambda_monitor() SETTING parameters. */
char* filename = mccpollambdaMonitorMSF_filename;
MCNUM mx = mccpollambdaMonitorMSF_mx;
MCNUM my = mccpollambdaMonitorMSF_my;
MCNUM mz = mccpollambdaMonitorMSF_mz;
MCNUM Lmin = mccpollambdaMonitorMSF_Lmin;
MCNUM Lmax = mccpollambdaMonitorMSF_Lmax;
#line 135 "/usr/share/mcstas/2.5/monitors/PolLambda_monitor.comp"
{
    if (!nowritefile) {
  DETECTOR_OUT_2D("Pol-wavelength monitor",
		  "Wavelength [AA]", "Polarisation projection",
		  Lmin, Lmax, -1.0, 1.0,
		  nL, npol,
		  &PolL_N[0][0],&PolL_p[0][0],&PolL_p2[0][0],
		  filename);
    }
}
#line 14021 "Test_Pol_TripleAxis.c"
}   /* End of pollambdaMonitorMSF=PolLambda_monitor() SETTING parameter declarations. */
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef npol
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'pollambdaMonitorDet'. */
  SIG_MESSAGE("pollambdaMonitorDet (Save)");
#define mccompcurname  pollambdaMonitorDet
#define mccompcurtype  PolLambda_monitor
#define mccompcurindex 18
#define xwidth mccpollambdaMonitorDet_xwidth
#define yheight mccpollambdaMonitorDet_yheight
#define nL mccpollambdaMonitorDet_nL
#define npol mccpollambdaMonitorDet_npol
#define restore_neutron mccpollambdaMonitorDet_restore_neutron
#define nowritefile mccpollambdaMonitorDet_nowritefile
#define PolL_N mccpollambdaMonitorDet_PolL_N
#define PolL_p mccpollambdaMonitorDet_PolL_p
#define PolL_p2 mccpollambdaMonitorDet_PolL_p2
{   /* Declarations of pollambdaMonitorDet=PolLambda_monitor() SETTING parameters. */
char* filename = mccpollambdaMonitorDet_filename;
MCNUM mx = mccpollambdaMonitorDet_mx;
MCNUM my = mccpollambdaMonitorDet_my;
MCNUM mz = mccpollambdaMonitorDet_mz;
MCNUM Lmin = mccpollambdaMonitorDet_Lmin;
MCNUM Lmax = mccpollambdaMonitorDet_Lmax;
#line 135 "/usr/share/mcstas/2.5/monitors/PolLambda_monitor.comp"
{
    if (!nowritefile) {
  DETECTOR_OUT_2D("Pol-wavelength monitor",
		  "Wavelength [AA]", "Polarisation projection",
		  Lmin, Lmax, -1.0, 1.0,
		  nL, npol,
		  &PolL_N[0][0],&PolL_p[0][0],&PolL_p2[0][0],
		  filename);
    }
}
#line 14068 "Test_Pol_TripleAxis.c"
}   /* End of pollambdaMonitorDet=PolLambda_monitor() SETTING parameter declarations. */
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef npol
#undef nL
#undef yheight
#undef xwidth
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
#line 14117 "Test_Pol_TripleAxis.c"
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
    if (!mcNCounter[3]) fprintf(stderr, "Warning: No neutron could reach Component[3] Mono_Arm\n");
    if (mcAbsorbProp[3]) fprintf(stderr, "Warning: %g events were removed in Component[3] Mono_Arm=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[3]);
    if (!mcNCounter[4]) fprintf(stderr, "Warning: No neutron could reach Component[4] Mono_Out\n");
    if (mcAbsorbProp[4]) fprintf(stderr, "Warning: %g events were removed in Component[4] Mono_Out=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[4]);
    if (!mcNCounter[5]) fprintf(stderr, "Warning: No neutron could reach Component[5] Mono\n");
    if (mcAbsorbProp[5]) fprintf(stderr, "Warning: %g events were removed in Component[5] Mono=Monochromator_flat()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[5]);
    if (!mcNCounter[6]) fprintf(stderr, "Warning: No neutron could reach Component[6] Mono_pol\n");
    if (mcAbsorbProp[6]) fprintf(stderr, "Warning: %g events were removed in Component[6] Mono_pol=Monochromator_pol()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[6]);
    if (!mcNCounter[7]) fprintf(stderr, "Warning: No neutron could reach Component[7] pollambdaMonitorMono\n");
    if (mcAbsorbProp[7]) fprintf(stderr, "Warning: %g events were removed in Component[7] pollambdaMonitorMono=PolLambda_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[7]);
    if (!mcNCounter[8]) fprintf(stderr, "Warning: No neutron could reach Component[8] Sample_Arm\n");
    if (mcAbsorbProp[8]) fprintf(stderr, "Warning: %g events were removed in Component[8] Sample_Arm=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[8]);
    if (!mcNCounter[9]) fprintf(stderr, "Warning: No neutron could reach Component[9] Sample_Out\n");
    if (mcAbsorbProp[9]) fprintf(stderr, "Warning: %g events were removed in Component[9] Sample_Out=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[9]);
    if (!mcNCounter[10]) fprintf(stderr, "Warning: No neutron could reach Component[10] Sample\n");
    if (mcAbsorbProp[10]) fprintf(stderr, "Warning: %g events were removed in Component[10] Sample=V_sample()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[10]);
    if (!mcNCounter[11]) fprintf(stderr, "Warning: No neutron could reach Component[11] pollambdaMonitorSample\n");
    if (mcAbsorbProp[11]) fprintf(stderr, "Warning: %g events were removed in Component[11] pollambdaMonitorSample=PolLambda_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[11]);
    if (!mcNCounter[12]) fprintf(stderr, "Warning: No neutron could reach Component[12] msf\n");
    if (mcAbsorbProp[12]) fprintf(stderr, "Warning: %g events were removed in Component[12] msf=Pol_constBfield()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[12]);
    if (!mcNCounter[13]) fprintf(stderr, "Warning: No neutron could reach Component[13] pollambdaMonitorMSF\n");
    if (mcAbsorbProp[13]) fprintf(stderr, "Warning: %g events were removed in Component[13] pollambdaMonitorMSF=PolLambda_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[13]);
    if (!mcNCounter[14]) fprintf(stderr, "Warning: No neutron could reach Component[14] Ana_Arm\n");
    if (mcAbsorbProp[14]) fprintf(stderr, "Warning: %g events were removed in Component[14] Ana_Arm=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[14]);
    if (!mcNCounter[15]) fprintf(stderr, "Warning: No neutron could reach Component[15] Ana_Out\n");
    if (mcAbsorbProp[15]) fprintf(stderr, "Warning: %g events were removed in Component[15] Ana_Out=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[15]);
    if (!mcNCounter[16]) fprintf(stderr, "Warning: No neutron could reach Component[16] Ana\n");
    if (mcAbsorbProp[16]) fprintf(stderr, "Warning: %g events were removed in Component[16] Ana=Monochromator_flat()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[16]);
    if (!mcNCounter[17]) fprintf(stderr, "Warning: No neutron could reach Component[17] Ana_pol\n");
    if (mcAbsorbProp[17]) fprintf(stderr, "Warning: %g events were removed in Component[17] Ana_pol=Monochromator_pol()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[17]);
    if (!mcNCounter[18]) fprintf(stderr, "Warning: No neutron could reach Component[18] pollambdaMonitorDet\n");
    if (mcAbsorbProp[18]) fprintf(stderr, "Warning: %g events were removed in Component[18] pollambdaMonitorDet=PolLambda_monitor()\n"
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
#line 147 "/usr/share/mcstas/2.5/misc/Progress_bar.comp"
{
  
}
#line 14197 "Test_Pol_TripleAxis.c"
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
#line 14246 "Test_Pol_TripleAxis.c"
}   /* End of Source=Source_simple() SETTING parameter declarations. */
#undef srcArea
#undef square
#undef pmul
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Mono_Arm'. */
  SIG_MESSAGE("Mono_Arm (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Mono_Arm");
#define mccompcurname  Mono_Arm
#define mccompcurtype  Arm
#define mccompcurindex 3
#line 40 "/usr/share/mcstas/2.5/optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 14269 "Test_Pol_TripleAxis.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Mono_Out'. */
  SIG_MESSAGE("Mono_Out (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Mono_Out");
#define mccompcurname  Mono_Out
#define mccompcurtype  Arm
#define mccompcurindex 4
#line 40 "/usr/share/mcstas/2.5/optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 14288 "Test_Pol_TripleAxis.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Mono'. */
  SIG_MESSAGE("Mono (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Mono");
#define mccompcurname  Mono
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 5
#define mos_rms_y mccMono_mos_rms_y
#define mos_rms_z mccMono_mos_rms_z
#define mos_rms_max mccMono_mos_rms_max
#define mono_Q mccMono_mono_Q
{   /* Declarations of Mono=Monochromator_flat() SETTING parameters. */
MCNUM zmin = mccMono_zmin;
MCNUM zmax = mccMono_zmax;
MCNUM ymin = mccMono_ymin;
MCNUM ymax = mccMono_ymax;
MCNUM zwidth = mccMono_zwidth;
MCNUM yheight = mccMono_yheight;
MCNUM mosaich = mccMono_mosaich;
MCNUM mosaicv = mccMono_mosaicv;
MCNUM r0 = mccMono_r0;
MCNUM Q = mccMono_Q;
MCNUM DM = mccMono_DM;
#line 254 "/usr/share/mcstas/2.5/optics/Monochromator_flat.comp"
{
  
  multiline(5, 0.0, (double)ymin, (double)zmin,
               0.0, (double)ymax, (double)zmin,
               0.0, (double)ymax, (double)zmax,
               0.0, (double)ymin, (double)zmax,
               0.0, (double)ymin, (double)zmin);
}
#line 14324 "Test_Pol_TripleAxis.c"
}   /* End of Mono=Monochromator_flat() SETTING parameter declarations. */
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Mono_pol'. */
  SIG_MESSAGE("Mono_pol (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Mono_pol");
#define mccompcurname  Mono_pol
#define mccompcurtype  Monochromator_pol
#define mccompcurindex 6
#define mos_rms mccMono_pol_mos_rms
#define d_rms mccMono_pol_d_rms
#define mono_Q mccMono_pol_mono_Q
#define FN mccMono_pol_FN
#define FM mccMono_pol_FM
{   /* Declarations of Mono_pol=Monochromator_pol() SETTING parameters. */
MCNUM zwidth = mccMono_pol_zwidth;
MCNUM yheight = mccMono_pol_yheight;
MCNUM mosaic = mccMono_pol_mosaic;
MCNUM dspread = mccMono_pol_dspread;
MCNUM Q = mccMono_pol_Q;
MCNUM DM = mccMono_pol_DM;
MCNUM pThreshold = mccMono_pol_pThreshold;
MCNUM Rup = mccMono_pol_Rup;
MCNUM Rdown = mccMono_pol_Rdown;
int debug = mccMono_pol_debug;
#line 199 "/usr/share/mcstas/2.5/optics/Monochromator_pol.comp"
{
  
  rectangle("yz", 0, 0, 0, zwidth, yheight);
}
#line 14361 "Test_Pol_TripleAxis.c"
}   /* End of Mono_pol=Monochromator_pol() SETTING parameter declarations. */
#undef FM
#undef FN
#undef mono_Q
#undef d_rms
#undef mos_rms
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'pollambdaMonitorMono'. */
  SIG_MESSAGE("pollambdaMonitorMono (McDisplay)");
  printf("MCDISPLAY: component %s\n", "pollambdaMonitorMono");
#define mccompcurname  pollambdaMonitorMono
#define mccompcurtype  PolLambda_monitor
#define mccompcurindex 7
#define xwidth mccpollambdaMonitorMono_xwidth
#define yheight mccpollambdaMonitorMono_yheight
#define nL mccpollambdaMonitorMono_nL
#define npol mccpollambdaMonitorMono_npol
#define restore_neutron mccpollambdaMonitorMono_restore_neutron
#define nowritefile mccpollambdaMonitorMono_nowritefile
#define PolL_N mccpollambdaMonitorMono_PolL_N
#define PolL_p mccpollambdaMonitorMono_PolL_p
#define PolL_p2 mccpollambdaMonitorMono_PolL_p2
{   /* Declarations of pollambdaMonitorMono=PolLambda_monitor() SETTING parameters. */
char* filename = mccpollambdaMonitorMono_filename;
MCNUM mx = mccpollambdaMonitorMono_mx;
MCNUM my = mccpollambdaMonitorMono_my;
MCNUM mz = mccpollambdaMonitorMono_mz;
MCNUM Lmin = mccpollambdaMonitorMono_Lmin;
MCNUM Lmax = mccpollambdaMonitorMono_Lmax;
#line 147 "/usr/share/mcstas/2.5/monitors/PolLambda_monitor.comp"
{
  
  rectangle("xy", 0, 0, 0, xwidth, yheight);
}
#line 14399 "Test_Pol_TripleAxis.c"
}   /* End of pollambdaMonitorMono=PolLambda_monitor() SETTING parameter declarations. */
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef npol
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Sample_Arm'. */
  SIG_MESSAGE("Sample_Arm (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Sample_Arm");
#define mccompcurname  Sample_Arm
#define mccompcurtype  Arm
#define mccompcurindex 8
#line 40 "/usr/share/mcstas/2.5/optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 14428 "Test_Pol_TripleAxis.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Sample_Out'. */
  SIG_MESSAGE("Sample_Out (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Sample_Out");
#define mccompcurname  Sample_Out
#define mccompcurtype  Arm
#define mccompcurindex 9
#line 40 "/usr/share/mcstas/2.5/optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 14447 "Test_Pol_TripleAxis.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Sample'. */
  SIG_MESSAGE("Sample (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Sample");
#define mccompcurname  Sample
#define mccompcurtype  V_sample
#define mccompcurindex 10
#define VarsV mccSample_VarsV
{   /* Declarations of Sample=V_sample() SETTING parameters. */
MCNUM radius = mccSample_radius;
MCNUM thickness = mccSample_thickness;
MCNUM zdepth = mccSample_zdepth;
MCNUM Vc = mccSample_Vc;
MCNUM sigma_abs = mccSample_sigma_abs;
MCNUM sigma_inc = mccSample_sigma_inc;
MCNUM radius_i = mccSample_radius_i;
MCNUM radius_o = mccSample_radius_o;
MCNUM h = mccSample_h;
MCNUM focus_r = mccSample_focus_r;
MCNUM pack = mccSample_pack;
MCNUM frac = mccSample_frac;
MCNUM f_QE = mccSample_f_QE;
MCNUM gamma = mccSample_gamma;
MCNUM target_x = mccSample_target_x;
MCNUM target_y = mccSample_target_y;
MCNUM target_z = mccSample_target_z;
MCNUM focus_xw = mccSample_focus_xw;
MCNUM focus_yh = mccSample_focus_yh;
MCNUM focus_aw = mccSample_focus_aw;
MCNUM focus_ah = mccSample_focus_ah;
MCNUM xwidth = mccSample_xwidth;
MCNUM yheight = mccSample_yheight;
MCNUM zthick = mccSample_zthick;
MCNUM rad_sphere = mccSample_rad_sphere;
MCNUM sig_a = mccSample_sig_a;
MCNUM sig_i = mccSample_sig_i;
MCNUM V0 = mccSample_V0;
int target_index = mccSample_target_index;
MCNUM multiples = mccSample_multiples;
#line 320 "/usr/share/mcstas/2.5/obsolete/V_sample.comp"
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
#line 14537 "Test_Pol_TripleAxis.c"
}   /* End of Sample=V_sample() SETTING parameter declarations. */
#undef VarsV
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'pollambdaMonitorSample'. */
  SIG_MESSAGE("pollambdaMonitorSample (McDisplay)");
  printf("MCDISPLAY: component %s\n", "pollambdaMonitorSample");
#define mccompcurname  pollambdaMonitorSample
#define mccompcurtype  PolLambda_monitor
#define mccompcurindex 11
#define xwidth mccpollambdaMonitorSample_xwidth
#define yheight mccpollambdaMonitorSample_yheight
#define nL mccpollambdaMonitorSample_nL
#define npol mccpollambdaMonitorSample_npol
#define restore_neutron mccpollambdaMonitorSample_restore_neutron
#define nowritefile mccpollambdaMonitorSample_nowritefile
#define PolL_N mccpollambdaMonitorSample_PolL_N
#define PolL_p mccpollambdaMonitorSample_PolL_p
#define PolL_p2 mccpollambdaMonitorSample_PolL_p2
{   /* Declarations of pollambdaMonitorSample=PolLambda_monitor() SETTING parameters. */
char* filename = mccpollambdaMonitorSample_filename;
MCNUM mx = mccpollambdaMonitorSample_mx;
MCNUM my = mccpollambdaMonitorSample_my;
MCNUM mz = mccpollambdaMonitorSample_mz;
MCNUM Lmin = mccpollambdaMonitorSample_Lmin;
MCNUM Lmax = mccpollambdaMonitorSample_Lmax;
#line 147 "/usr/share/mcstas/2.5/monitors/PolLambda_monitor.comp"
{
  
  rectangle("xy", 0, 0, 0, xwidth, yheight);
}
#line 14571 "Test_Pol_TripleAxis.c"
}   /* End of pollambdaMonitorSample=PolLambda_monitor() SETTING parameter declarations. */
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef npol
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'msf'. */
  SIG_MESSAGE("msf (McDisplay)");
  printf("MCDISPLAY: component %s\n", "msf");
#define mccompcurname  msf
#define mccompcurtype  Pol_constBfield
#define mccompcurindex 12
#define xwidth mccmsf_xwidth
#define yheight mccmsf_yheight
#define length mccmsf_length
#define omegaL mccmsf_omegaL
{   /* Declarations of msf=Pol_constBfield() SETTING parameters. */
MCNUM B = mccmsf_B;
MCNUM fliplambda = mccmsf_fliplambda;
MCNUM flipangle = mccmsf_flipangle;
#line 149 "/usr/share/mcstas/2.5/optics/Pol_constBfield.comp"
{
  

  box(0, 0, length/2, xwidth, yheight, length);
}
#line 14606 "Test_Pol_TripleAxis.c"
}   /* End of msf=Pol_constBfield() SETTING parameter declarations. */
#undef omegaL
#undef length
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'pollambdaMonitorMSF'. */
  SIG_MESSAGE("pollambdaMonitorMSF (McDisplay)");
  printf("MCDISPLAY: component %s\n", "pollambdaMonitorMSF");
#define mccompcurname  pollambdaMonitorMSF
#define mccompcurtype  PolLambda_monitor
#define mccompcurindex 13
#define xwidth mccpollambdaMonitorMSF_xwidth
#define yheight mccpollambdaMonitorMSF_yheight
#define nL mccpollambdaMonitorMSF_nL
#define npol mccpollambdaMonitorMSF_npol
#define restore_neutron mccpollambdaMonitorMSF_restore_neutron
#define nowritefile mccpollambdaMonitorMSF_nowritefile
#define PolL_N mccpollambdaMonitorMSF_PolL_N
#define PolL_p mccpollambdaMonitorMSF_PolL_p
#define PolL_p2 mccpollambdaMonitorMSF_PolL_p2
{   /* Declarations of pollambdaMonitorMSF=PolLambda_monitor() SETTING parameters. */
char* filename = mccpollambdaMonitorMSF_filename;
MCNUM mx = mccpollambdaMonitorMSF_mx;
MCNUM my = mccpollambdaMonitorMSF_my;
MCNUM mz = mccpollambdaMonitorMSF_mz;
MCNUM Lmin = mccpollambdaMonitorMSF_Lmin;
MCNUM Lmax = mccpollambdaMonitorMSF_Lmax;
#line 147 "/usr/share/mcstas/2.5/monitors/PolLambda_monitor.comp"
{
  
  rectangle("xy", 0, 0, 0, xwidth, yheight);
}
#line 14643 "Test_Pol_TripleAxis.c"
}   /* End of pollambdaMonitorMSF=PolLambda_monitor() SETTING parameter declarations. */
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef npol
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Ana_Arm'. */
  SIG_MESSAGE("Ana_Arm (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Ana_Arm");
#define mccompcurname  Ana_Arm
#define mccompcurtype  Arm
#define mccompcurindex 14
#line 40 "/usr/share/mcstas/2.5/optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 14672 "Test_Pol_TripleAxis.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Ana_Out'. */
  SIG_MESSAGE("Ana_Out (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Ana_Out");
#define mccompcurname  Ana_Out
#define mccompcurtype  Arm
#define mccompcurindex 15
#line 40 "/usr/share/mcstas/2.5/optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 14691 "Test_Pol_TripleAxis.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Ana'. */
  SIG_MESSAGE("Ana (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Ana");
#define mccompcurname  Ana
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 16
#define mos_rms_y mccAna_mos_rms_y
#define mos_rms_z mccAna_mos_rms_z
#define mos_rms_max mccAna_mos_rms_max
#define mono_Q mccAna_mono_Q
{   /* Declarations of Ana=Monochromator_flat() SETTING parameters. */
MCNUM zmin = mccAna_zmin;
MCNUM zmax = mccAna_zmax;
MCNUM ymin = mccAna_ymin;
MCNUM ymax = mccAna_ymax;
MCNUM zwidth = mccAna_zwidth;
MCNUM yheight = mccAna_yheight;
MCNUM mosaich = mccAna_mosaich;
MCNUM mosaicv = mccAna_mosaicv;
MCNUM r0 = mccAna_r0;
MCNUM Q = mccAna_Q;
MCNUM DM = mccAna_DM;
#line 254 "/usr/share/mcstas/2.5/optics/Monochromator_flat.comp"
{
  
  multiline(5, 0.0, (double)ymin, (double)zmin,
               0.0, (double)ymax, (double)zmin,
               0.0, (double)ymax, (double)zmax,
               0.0, (double)ymin, (double)zmax,
               0.0, (double)ymin, (double)zmin);
}
#line 14727 "Test_Pol_TripleAxis.c"
}   /* End of Ana=Monochromator_flat() SETTING parameter declarations. */
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Ana_pol'. */
  SIG_MESSAGE("Ana_pol (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Ana_pol");
#define mccompcurname  Ana_pol
#define mccompcurtype  Monochromator_pol
#define mccompcurindex 17
#define mos_rms mccAna_pol_mos_rms
#define d_rms mccAna_pol_d_rms
#define mono_Q mccAna_pol_mono_Q
#define FN mccAna_pol_FN
#define FM mccAna_pol_FM
{   /* Declarations of Ana_pol=Monochromator_pol() SETTING parameters. */
MCNUM zwidth = mccAna_pol_zwidth;
MCNUM yheight = mccAna_pol_yheight;
MCNUM mosaic = mccAna_pol_mosaic;
MCNUM dspread = mccAna_pol_dspread;
MCNUM Q = mccAna_pol_Q;
MCNUM DM = mccAna_pol_DM;
MCNUM pThreshold = mccAna_pol_pThreshold;
MCNUM Rup = mccAna_pol_Rup;
MCNUM Rdown = mccAna_pol_Rdown;
int debug = mccAna_pol_debug;
#line 199 "/usr/share/mcstas/2.5/optics/Monochromator_pol.comp"
{
  
  rectangle("yz", 0, 0, 0, zwidth, yheight);
}
#line 14764 "Test_Pol_TripleAxis.c"
}   /* End of Ana_pol=Monochromator_pol() SETTING parameter declarations. */
#undef FM
#undef FN
#undef mono_Q
#undef d_rms
#undef mos_rms
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'pollambdaMonitorDet'. */
  SIG_MESSAGE("pollambdaMonitorDet (McDisplay)");
  printf("MCDISPLAY: component %s\n", "pollambdaMonitorDet");
#define mccompcurname  pollambdaMonitorDet
#define mccompcurtype  PolLambda_monitor
#define mccompcurindex 18
#define xwidth mccpollambdaMonitorDet_xwidth
#define yheight mccpollambdaMonitorDet_yheight
#define nL mccpollambdaMonitorDet_nL
#define npol mccpollambdaMonitorDet_npol
#define restore_neutron mccpollambdaMonitorDet_restore_neutron
#define nowritefile mccpollambdaMonitorDet_nowritefile
#define PolL_N mccpollambdaMonitorDet_PolL_N
#define PolL_p mccpollambdaMonitorDet_PolL_p
#define PolL_p2 mccpollambdaMonitorDet_PolL_p2
{   /* Declarations of pollambdaMonitorDet=PolLambda_monitor() SETTING parameters. */
char* filename = mccpollambdaMonitorDet_filename;
MCNUM mx = mccpollambdaMonitorDet_mx;
MCNUM my = mccpollambdaMonitorDet_my;
MCNUM mz = mccpollambdaMonitorDet_mz;
MCNUM Lmin = mccpollambdaMonitorDet_Lmin;
MCNUM Lmax = mccpollambdaMonitorDet_Lmax;
#line 147 "/usr/share/mcstas/2.5/monitors/PolLambda_monitor.comp"
{
  
  rectangle("xy", 0, 0, 0, xwidth, yheight);
}
#line 14802 "Test_Pol_TripleAxis.c"
}   /* End of pollambdaMonitorDet=PolLambda_monitor() SETTING parameter declarations. */
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef npol
#undef nL
#undef yheight
#undef xwidth
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
/* end of generated C code Test_Pol_TripleAxis.c */
