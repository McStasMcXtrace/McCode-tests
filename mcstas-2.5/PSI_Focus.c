/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: PSI_Focus.instr (PSI_Focus)
 * Date:       Fri Oct 11 20:08:59 2019
 * File:       PSI_Focus.c
 * Compile:    cc -o PSI_Focus.out PSI_Focus.c 
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

#line 695 "PSI_Focus.c"

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

#line 928 "PSI_Focus.c"

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

#line 4947 "PSI_Focus.c"

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

#line 5307 "PSI_Focus.c"
#ifdef MC_TRACE_ENABLED
int mctraceenabled = 1;
#else
int mctraceenabled = 0;
#endif
#define MCSTAS "/usr/share/mcstas/2.5/"
int mcdefaultmain = 1;
char mcinstrument_name[] = "PSI_Focus";
char mcinstrument_source[] = "PSI_Focus.instr";
char *mcinstrument_exe=NULL; /* will be set to argv[0] in main */
int main(int argc, char *argv[]){return mccode_main(argc, argv);}
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);

/* Shared user declarations for all components 'Source_gen'. */
#line 140 "/usr/share/mcstas/2.5/sources/Source_gen.comp"
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

#line 6803 "PSI_Focus.c"

/* Shared user declarations for all components 'Guide'. */
#line 63 "/usr/share/mcstas/2.5/optics/Guide.comp"

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

#line 6956 "PSI_Focus.c"

/* Shared user declarations for all components 'Bender'. */
#line 102 "/usr/share/mcstas/2.5/optics/Bender.comp"

#line 6961 "PSI_Focus.c"

/* Shared user declarations for all components 'Guide_channeled'. */
#line 76 "/usr/share/mcstas/2.5/optics/Guide_channeled.comp"

#line 6966 "PSI_Focus.c"

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

#line 9901 "PSI_Focus.c"

/* Shared user declarations for all components 'Monochromator_2foc'. */
#line 83 "/usr/share/mcstas/2.5/contrib/Monochromator_2foc.comp"

#line 9906 "PSI_Focus.c"

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

#line 10086 "PSI_Focus.c"

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
#line 10103 "PSI_Focus.c"

/* Instrument parameters. */
MCNUM mciplambda;
MCNUM mcipchopp_ratio;
MCNUM mcipDET;

#define mcNUMIPAR 3
int mcnumipar = 3;
struct mcinputtable_struct mcinputtable[mcNUMIPAR+1] = {
  "lambda", &mciplambda, instr_type_double, "3.4", 
  "chopp_ratio", &mcipchopp_ratio, instr_type_double, "1", 
  "DET", &mcipDET, instr_type_double, "-69.9", 
  NULL, NULL, instr_type_double, ""
};

/* User declarations from instrument definition. */
#define mccompcurname  PSI_Focus
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaPSI_Focus coords_set(0,0,0)
#define lambda mciplambda
#define chopp_ratio mcipchopp_ratio
#define DET mcipDET
#line 31 "PSI_Focus.instr"
 double lambda;
 double chopp_ratio;
 double DET;
 double PHM, TTM;
 double RV_2, RH_2;
 double FO_PHA;
 double DISC_SPEED, FERMI_SPEED;
 /* mono d-spacing */
 double monodd = 3.355;
 /* time constant [mus/A/m] */
 double c = 252.78;
 /* distance sample - detector */
 double dsd = 2.5;
 /* distance guide - monochromator */
 double dgm = 2.996;
 /* distance fermi-chopper - sample */
 double dfcs = 0.4997;
 /* distance disk chopper - monochromator */
 double ddcm = 2.944;
 /* distance  monochromator - fermi-chopper */
 double mfc = 1.002;
 /* distance  monochromator - sample */
 double dms = 1.4997;
 /* setting range of sources */
 double dL = 0.7;
 double LMIN, LMAX;
 double EMIN, EMAX;

#line 10156 "PSI_Focus.c"
#undef DET
#undef chopp_ratio
#undef lambda
#undef mcposaPSI_Focus
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname

/* neutron state table at each component input (local coords) */
/* [x, y, z, vx, vy, vz, t, sx, sy, sz, p] */
MCNUM mccomp_storein[11*46];
/* Components position table (absolute and relative coords) */
Coords mccomp_posa[46];
Coords mccomp_posr[46];
/* Counter for each comp to check for inactive ones */
MCNUM  mcNCounter[46];
MCNUM  mcPCounter[46];
MCNUM  mcP2Counter[46];
#define mcNUMCOMP 45 /* number of components */
/* Counter for PROP ABSORB */
MCNUM  mcAbsorbProp[46];
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

/* Setting parameters for component 'csource' [2]. */
char mcccsource_flux_file[16384];
char mcccsource_xdiv_file[16384];
char mcccsource_ydiv_file[16384];
MCNUM mcccsource_radius;
MCNUM mcccsource_dist;
MCNUM mcccsource_focus_xw;
MCNUM mcccsource_focus_yh;
MCNUM mcccsource_focus_aw;
MCNUM mcccsource_focus_ah;
MCNUM mcccsource_E0;
MCNUM mcccsource_dE;
MCNUM mcccsource_lambda0;
MCNUM mcccsource_dlambda;
MCNUM mcccsource_I1;
MCNUM mcccsource_yheight;
MCNUM mcccsource_xwidth;
MCNUM mcccsource_verbose;
MCNUM mcccsource_T1;
MCNUM mcccsource_flux_file_perAA;
MCNUM mcccsource_flux_file_log;
MCNUM mcccsource_Lmin;
MCNUM mcccsource_Lmax;
MCNUM mcccsource_Emin;
MCNUM mcccsource_Emax;
MCNUM mcccsource_T2;
MCNUM mcccsource_I2;
MCNUM mcccsource_T3;
MCNUM mcccsource_I3;
MCNUM mcccsource_zdepth;
int mcccsource_target_index;

/* Setting parameters for component 'guide1' [3]. */
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

/* Setting parameters for component 'guide2' [4]. */
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

/* Setting parameters for component 'bunker' [5]. */
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

/* Setting parameters for component 'guide3' [6]. */
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

/* Definition parameters for component 'lambdaGuideExit' [7]. */
#define mcclambdaGuideExit_nL 60
/* Setting parameters for component 'lambdaGuideExit' [7]. */
char mcclambdaGuideExit_filename[16384];
MCNUM mcclambdaGuideExit_xmin;
MCNUM mcclambdaGuideExit_xmax;
MCNUM mcclambdaGuideExit_ymin;
MCNUM mcclambdaGuideExit_ymax;
MCNUM mcclambdaGuideExit_xwidth;
MCNUM mcclambdaGuideExit_yheight;
MCNUM mcclambdaGuideExit_Lmin;
MCNUM mcclambdaGuideExit_Lmax;
MCNUM mcclambdaGuideExit_restore_neutron;
int mcclambdaGuideExit_nowritefile;

/* Definition parameters for component 'DivMonGuideExit' [8]. */
#define mccDivMonGuideExit_nh 60
#define mccDivMonGuideExit_nv 60
/* Setting parameters for component 'DivMonGuideExit' [8]. */
char mccDivMonGuideExit_filename[16384];
MCNUM mccDivMonGuideExit_xmin;
MCNUM mccDivMonGuideExit_xmax;
MCNUM mccDivMonGuideExit_ymin;
MCNUM mccDivMonGuideExit_ymax;
MCNUM mccDivMonGuideExit_xwidth;
MCNUM mccDivMonGuideExit_yheight;
MCNUM mccDivMonGuideExit_maxdiv_h;
MCNUM mccDivMonGuideExit_maxdiv_v;
MCNUM mccDivMonGuideExit_restore_neutron;
MCNUM mccDivMonGuideExit_nx;
MCNUM mccDivMonGuideExit_ny;
MCNUM mccDivMonGuideExit_nz;
int mccDivMonGuideExit_nowritefile;

/* Definition parameters for component 'PSDGuideExit' [9]. */
#define mccPSDGuideExit_nx 60
#define mccPSDGuideExit_ny 60
/* Setting parameters for component 'PSDGuideExit' [9]. */
char mccPSDGuideExit_filename[16384];
MCNUM mccPSDGuideExit_xmin;
MCNUM mccPSDGuideExit_xmax;
MCNUM mccPSDGuideExit_ymin;
MCNUM mccPSDGuideExit_ymax;
MCNUM mccPSDGuideExit_xwidth;
MCNUM mccPSDGuideExit_yheight;
MCNUM mccPSDGuideExit_restore_neutron;
int mccPSDGuideExit_nowritefile;

/* Setting parameters for component 'FOCUSguide' [10]. */
MCNUM mccFOCUSguide_w1;
MCNUM mccFOCUSguide_h1;
MCNUM mccFOCUSguide_w2;
MCNUM mccFOCUSguide_h2;
MCNUM mccFOCUSguide_l;
MCNUM mccFOCUSguide_R0;
MCNUM mccFOCUSguide_Qc;
MCNUM mccFOCUSguide_alpha;
MCNUM mccFOCUSguide_m;
MCNUM mccFOCUSguide_nslit;
MCNUM mccFOCUSguide_d;
MCNUM mccFOCUSguide_Qcx;
MCNUM mccFOCUSguide_Qcy;
MCNUM mccFOCUSguide_alphax;
MCNUM mccFOCUSguide_alphay;
MCNUM mccFOCUSguide_W;
MCNUM mccFOCUSguide_mx;
MCNUM mccFOCUSguide_my;
MCNUM mccFOCUSguide_nu;
MCNUM mccFOCUSguide_phase;

/* Setting parameters for component 'FirstChopper' [11]. */
MCNUM mccFirstChopper_theta_0;
MCNUM mccFirstChopper_radius;
MCNUM mccFirstChopper_yheight;
MCNUM mccFirstChopper_nu;
MCNUM mccFirstChopper_nslit;
MCNUM mccFirstChopper_jitter;
MCNUM mccFirstChopper_delay;
MCNUM mccFirstChopper_isfirst;
MCNUM mccFirstChopper_n_pulse;
MCNUM mccFirstChopper_abs_out;
MCNUM mccFirstChopper_phase;
MCNUM mccFirstChopper_xwidth;
MCNUM mccFirstChopper_verbose;

/* Definition parameters for component 'DISCTOF' [12]. */
#define mccDISCTOF_user1 FLT_MAX
#define mccDISCTOF_user2 FLT_MAX
#define mccDISCTOF_user3 FLT_MAX
/* Setting parameters for component 'DISCTOF' [12]. */
MCNUM mccDISCTOF_xwidth;
MCNUM mccDISCTOF_yheight;
MCNUM mccDISCTOF_zdepth;
MCNUM mccDISCTOF_xmin;
MCNUM mccDISCTOF_xmax;
MCNUM mccDISCTOF_ymin;
MCNUM mccDISCTOF_ymax;
MCNUM mccDISCTOF_zmin;
MCNUM mccDISCTOF_zmax;
MCNUM mccDISCTOF_bins;
MCNUM mccDISCTOF_min;
MCNUM mccDISCTOF_max;
MCNUM mccDISCTOF_restore_neutron;
MCNUM mccDISCTOF_radius;
char mccDISCTOF_options[16384];
char mccDISCTOF_filename[16384];
char mccDISCTOF_geometry[16384];
char mccDISCTOF_username1[16384];
char mccDISCTOF_username2[16384];
char mccDISCTOF_username3[16384];
int mccDISCTOF_nowritefile;

/* Definition parameters for component 'PSDmon1Chopper' [13]. */
#define mccPSDmon1Chopper_nx 60
#define mccPSDmon1Chopper_ny 60
/* Setting parameters for component 'PSDmon1Chopper' [13]. */
char mccPSDmon1Chopper_filename[16384];
MCNUM mccPSDmon1Chopper_xmin;
MCNUM mccPSDmon1Chopper_xmax;
MCNUM mccPSDmon1Chopper_ymin;
MCNUM mccPSDmon1Chopper_ymax;
MCNUM mccPSDmon1Chopper_xwidth;
MCNUM mccPSDmon1Chopper_yheight;
MCNUM mccPSDmon1Chopper_restore_neutron;
int mccPSDmon1Chopper_nowritefile;

/* Setting parameters for component 'VacuumTube1entry' [14]. */
MCNUM mccVacuumTube1entry_xmin;
MCNUM mccVacuumTube1entry_xmax;
MCNUM mccVacuumTube1entry_ymin;
MCNUM mccVacuumTube1entry_ymax;
MCNUM mccVacuumTube1entry_radius;
MCNUM mccVacuumTube1entry_xwidth;
MCNUM mccVacuumTube1entry_yheight;

/* Setting parameters for component 'VacuumTube1exit' [15]. */
MCNUM mccVacuumTube1exit_xmin;
MCNUM mccVacuumTube1exit_xmax;
MCNUM mccVacuumTube1exit_ymin;
MCNUM mccVacuumTube1exit_ymax;
MCNUM mccVacuumTube1exit_radius;
MCNUM mccVacuumTube1exit_xwidth;
MCNUM mccVacuumTube1exit_yheight;

/* Setting parameters for component 'VacuumTube2entry' [16]. */
MCNUM mccVacuumTube2entry_xmin;
MCNUM mccVacuumTube2entry_xmax;
MCNUM mccVacuumTube2entry_ymin;
MCNUM mccVacuumTube2entry_ymax;
MCNUM mccVacuumTube2entry_radius;
MCNUM mccVacuumTube2entry_xwidth;
MCNUM mccVacuumTube2entry_yheight;

/* Setting parameters for component 'VacuumTube2exit' [17]. */
MCNUM mccVacuumTube2exit_xmin;
MCNUM mccVacuumTube2exit_xmax;
MCNUM mccVacuumTube2exit_ymin;
MCNUM mccVacuumTube2exit_ymax;
MCNUM mccVacuumTube2exit_radius;
MCNUM mccVacuumTube2exit_xwidth;
MCNUM mccVacuumTube2exit_yheight;

/* Setting parameters for component 'VacuumTube3entry' [18]. */
MCNUM mccVacuumTube3entry_xmin;
MCNUM mccVacuumTube3entry_xmax;
MCNUM mccVacuumTube3entry_ymin;
MCNUM mccVacuumTube3entry_ymax;
MCNUM mccVacuumTube3entry_radius;
MCNUM mccVacuumTube3entry_xwidth;
MCNUM mccVacuumTube3entry_yheight;

/* Setting parameters for component 'VacuumTube3exit' [19]. */
MCNUM mccVacuumTube3exit_xmin;
MCNUM mccVacuumTube3exit_xmax;
MCNUM mccVacuumTube3exit_ymin;
MCNUM mccVacuumTube3exit_ymax;
MCNUM mccVacuumTube3exit_radius;
MCNUM mccVacuumTube3exit_xwidth;
MCNUM mccVacuumTube3exit_yheight;

/* Definition parameters for component 'PSDmonMono' [20]. */
#define mccPSDmonMono_nx 60
#define mccPSDmonMono_ny 60
/* Setting parameters for component 'PSDmonMono' [20]. */
char mccPSDmonMono_filename[16384];
MCNUM mccPSDmonMono_xmin;
MCNUM mccPSDmonMono_xmax;
MCNUM mccPSDmonMono_ymin;
MCNUM mccPSDmonMono_ymax;
MCNUM mccPSDmonMono_xwidth;
MCNUM mccPSDmonMono_yheight;
MCNUM mccPSDmonMono_restore_neutron;
int mccPSDmonMono_nowritefile;

/* Definition parameters for component 'MONOTOF' [21]. */
#define mccMONOTOF_user1 FLT_MAX
#define mccMONOTOF_user2 FLT_MAX
#define mccMONOTOF_user3 FLT_MAX
/* Setting parameters for component 'MONOTOF' [21]. */
MCNUM mccMONOTOF_xwidth;
MCNUM mccMONOTOF_yheight;
MCNUM mccMONOTOF_zdepth;
MCNUM mccMONOTOF_xmin;
MCNUM mccMONOTOF_xmax;
MCNUM mccMONOTOF_ymin;
MCNUM mccMONOTOF_ymax;
MCNUM mccMONOTOF_zmin;
MCNUM mccMONOTOF_zmax;
MCNUM mccMONOTOF_bins;
MCNUM mccMONOTOF_min;
MCNUM mccMONOTOF_max;
MCNUM mccMONOTOF_restore_neutron;
MCNUM mccMONOTOF_radius;
char mccMONOTOF_options[16384];
char mccMONOTOF_filename[16384];
char mccMONOTOF_geometry[16384];
char mccMONOTOF_username1[16384];
char mccMONOTOF_username2[16384];
char mccMONOTOF_username3[16384];
int mccMONOTOF_nowritefile;

/* Definition parameters for component 'DivMonMono' [22]. */
#define mccDivMonMono_nh 60
#define mccDivMonMono_nv 60
/* Setting parameters for component 'DivMonMono' [22]. */
char mccDivMonMono_filename[16384];
MCNUM mccDivMonMono_xmin;
MCNUM mccDivMonMono_xmax;
MCNUM mccDivMonMono_ymin;
MCNUM mccDivMonMono_ymax;
MCNUM mccDivMonMono_xwidth;
MCNUM mccDivMonMono_yheight;
MCNUM mccDivMonMono_maxdiv_h;
MCNUM mccDivMonMono_maxdiv_v;
MCNUM mccDivMonMono_restore_neutron;
MCNUM mccDivMonMono_nx;
MCNUM mccDivMonMono_ny;
MCNUM mccDivMonMono_nz;
int mccDivMonMono_nowritefile;

/* Setting parameters for component 'mono' [24]. */
char mccmono_reflect[16384];
MCNUM mccmono_zwidth;
MCNUM mccmono_yheight;
MCNUM mccmono_gap;
MCNUM mccmono_NH;
MCNUM mccmono_NV;
MCNUM mccmono_mosaich;
MCNUM mccmono_mosaicv;
MCNUM mccmono_r0;
MCNUM mccmono_Q;
MCNUM mccmono_RV;
MCNUM mccmono_RH;
MCNUM mccmono_DM;
MCNUM mccmono_mosaic;
MCNUM mccmono_width;
MCNUM mccmono_height;
MCNUM mccmono_verbose;

/* Definition parameters for component 'FERMITOF_before' [26]. */
#define mccFERMITOF_before_user1 FLT_MAX
#define mccFERMITOF_before_user2 FLT_MAX
#define mccFERMITOF_before_user3 FLT_MAX
/* Setting parameters for component 'FERMITOF_before' [26]. */
MCNUM mccFERMITOF_before_xwidth;
MCNUM mccFERMITOF_before_yheight;
MCNUM mccFERMITOF_before_zdepth;
MCNUM mccFERMITOF_before_xmin;
MCNUM mccFERMITOF_before_xmax;
MCNUM mccFERMITOF_before_ymin;
MCNUM mccFERMITOF_before_ymax;
MCNUM mccFERMITOF_before_zmin;
MCNUM mccFERMITOF_before_zmax;
MCNUM mccFERMITOF_before_bins;
MCNUM mccFERMITOF_before_min;
MCNUM mccFERMITOF_before_max;
MCNUM mccFERMITOF_before_restore_neutron;
MCNUM mccFERMITOF_before_radius;
char mccFERMITOF_before_options[16384];
char mccFERMITOF_before_filename[16384];
char mccFERMITOF_before_geometry[16384];
char mccFERMITOF_before_username1[16384];
char mccFERMITOF_before_username2[16384];
char mccFERMITOF_before_username3[16384];
int mccFERMITOF_before_nowritefile;

/* Definition parameters for component 'lambdaFermi' [27]. */
#define mcclambdaFermi_nL 60
/* Setting parameters for component 'lambdaFermi' [27]. */
char mcclambdaFermi_filename[16384];
MCNUM mcclambdaFermi_xmin;
MCNUM mcclambdaFermi_xmax;
MCNUM mcclambdaFermi_ymin;
MCNUM mcclambdaFermi_ymax;
MCNUM mcclambdaFermi_xwidth;
MCNUM mcclambdaFermi_yheight;
MCNUM mcclambdaFermi_Lmin;
MCNUM mcclambdaFermi_Lmax;
MCNUM mcclambdaFermi_restore_neutron;
int mcclambdaFermi_nowritefile;

/* Definition parameters for component 'EMON_Fermi' [28]. */
#define mccEMON_Fermi_nE 60
/* Setting parameters for component 'EMON_Fermi' [28]. */
char mccEMON_Fermi_filename[16384];
MCNUM mccEMON_Fermi_xmin;
MCNUM mccEMON_Fermi_xmax;
MCNUM mccEMON_Fermi_ymin;
MCNUM mccEMON_Fermi_ymax;
MCNUM mccEMON_Fermi_xwidth;
MCNUM mccEMON_Fermi_yheight;
MCNUM mccEMON_Fermi_Emin;
MCNUM mccEMON_Fermi_Emax;
MCNUM mccEMON_Fermi_restore_neutron;
int mccEMON_Fermi_nowritefile;

/* Definition parameters for component 'DivMonfermi1' [29]. */
#define mccDivMonfermi1_nh 30
#define mccDivMonfermi1_nv 30
/* Setting parameters for component 'DivMonfermi1' [29]. */
char mccDivMonfermi1_filename[16384];
MCNUM mccDivMonfermi1_xmin;
MCNUM mccDivMonfermi1_xmax;
MCNUM mccDivMonfermi1_ymin;
MCNUM mccDivMonfermi1_ymax;
MCNUM mccDivMonfermi1_xwidth;
MCNUM mccDivMonfermi1_yheight;
MCNUM mccDivMonfermi1_maxdiv_h;
MCNUM mccDivMonfermi1_maxdiv_v;
MCNUM mccDivMonfermi1_restore_neutron;
MCNUM mccDivMonfermi1_nx;
MCNUM mccDivMonfermi1_ny;
MCNUM mccDivMonfermi1_nz;
int mccDivMonfermi1_nowritefile;

/* Definition parameters for component 'PSD_Fermi1' [30]. */
#define mccPSD_Fermi1_nx 30
#define mccPSD_Fermi1_ny 30
/* Setting parameters for component 'PSD_Fermi1' [30]. */
char mccPSD_Fermi1_filename[16384];
MCNUM mccPSD_Fermi1_xmin;
MCNUM mccPSD_Fermi1_xmax;
MCNUM mccPSD_Fermi1_ymin;
MCNUM mccPSD_Fermi1_ymax;
MCNUM mccPSD_Fermi1_xwidth;
MCNUM mccPSD_Fermi1_yheight;
MCNUM mccPSD_Fermi1_restore_neutron;
int mccPSD_Fermi1_nowritefile;

/* Setting parameters for component 'FoChopper' [31]. */
MCNUM mccFoChopper_phase;
MCNUM mccFoChopper_radius;
MCNUM mccFoChopper_nu;
MCNUM mccFoChopper_w;
MCNUM mccFoChopper_nslit;
MCNUM mccFoChopper_R0;
MCNUM mccFoChopper_Qc;
MCNUM mccFoChopper_alpha;
MCNUM mccFoChopper_m;
MCNUM mccFoChopper_W;
MCNUM mccFoChopper_length;
MCNUM mccFoChopper_eff;
MCNUM mccFoChopper_zero_time;
MCNUM mccFoChopper_xwidth;
MCNUM mccFoChopper_verbose;
MCNUM mccFoChopper_yheight;
MCNUM mccFoChopper_curvature;
MCNUM mccFoChopper_delay;

/* Definition parameters for component 'PSD_Fermi2' [32]. */
#define mccPSD_Fermi2_nx 30
#define mccPSD_Fermi2_ny 30
/* Setting parameters for component 'PSD_Fermi2' [32]. */
char mccPSD_Fermi2_filename[16384];
MCNUM mccPSD_Fermi2_xmin;
MCNUM mccPSD_Fermi2_xmax;
MCNUM mccPSD_Fermi2_ymin;
MCNUM mccPSD_Fermi2_ymax;
MCNUM mccPSD_Fermi2_xwidth;
MCNUM mccPSD_Fermi2_yheight;
MCNUM mccPSD_Fermi2_restore_neutron;
int mccPSD_Fermi2_nowritefile;

/* Definition parameters for component 'DivMonfermi2' [33]. */
#define mccDivMonfermi2_nh 30
#define mccDivMonfermi2_nv 30
/* Setting parameters for component 'DivMonfermi2' [33]. */
char mccDivMonfermi2_filename[16384];
MCNUM mccDivMonfermi2_xmin;
MCNUM mccDivMonfermi2_xmax;
MCNUM mccDivMonfermi2_ymin;
MCNUM mccDivMonfermi2_ymax;
MCNUM mccDivMonfermi2_xwidth;
MCNUM mccDivMonfermi2_yheight;
MCNUM mccDivMonfermi2_maxdiv_h;
MCNUM mccDivMonfermi2_maxdiv_v;
MCNUM mccDivMonfermi2_restore_neutron;
MCNUM mccDivMonfermi2_nx;
MCNUM mccDivMonfermi2_ny;
MCNUM mccDivMonfermi2_nz;
int mccDivMonfermi2_nowritefile;

/* Definition parameters for component 'FERMITOF1' [34]. */
#define mccFERMITOF1_user1 FLT_MAX
#define mccFERMITOF1_user2 FLT_MAX
#define mccFERMITOF1_user3 FLT_MAX
/* Setting parameters for component 'FERMITOF1' [34]. */
MCNUM mccFERMITOF1_xwidth;
MCNUM mccFERMITOF1_yheight;
MCNUM mccFERMITOF1_zdepth;
MCNUM mccFERMITOF1_xmin;
MCNUM mccFERMITOF1_xmax;
MCNUM mccFERMITOF1_ymin;
MCNUM mccFERMITOF1_ymax;
MCNUM mccFERMITOF1_zmin;
MCNUM mccFERMITOF1_zmax;
MCNUM mccFERMITOF1_bins;
MCNUM mccFERMITOF1_min;
MCNUM mccFERMITOF1_max;
MCNUM mccFERMITOF1_restore_neutron;
MCNUM mccFERMITOF1_radius;
char mccFERMITOF1_options[16384];
char mccFERMITOF1_filename[16384];
char mccFERMITOF1_geometry[16384];
char mccFERMITOF1_username1[16384];
char mccFERMITOF1_username2[16384];
char mccFERMITOF1_username3[16384];
int mccFERMITOF1_nowritefile;

/* Setting parameters for component 'SAMPLE_SLIT' [35]. */
MCNUM mccSAMPLE_SLIT_xmin;
MCNUM mccSAMPLE_SLIT_xmax;
MCNUM mccSAMPLE_SLIT_ymin;
MCNUM mccSAMPLE_SLIT_ymax;
MCNUM mccSAMPLE_SLIT_radius;
MCNUM mccSAMPLE_SLIT_xwidth;
MCNUM mccSAMPLE_SLIT_yheight;

/* Definition parameters for component 'FERMITOF2' [36]. */
#define mccFERMITOF2_user1 FLT_MAX
#define mccFERMITOF2_user2 FLT_MAX
#define mccFERMITOF2_user3 FLT_MAX
/* Setting parameters for component 'FERMITOF2' [36]. */
MCNUM mccFERMITOF2_xwidth;
MCNUM mccFERMITOF2_yheight;
MCNUM mccFERMITOF2_zdepth;
MCNUM mccFERMITOF2_xmin;
MCNUM mccFERMITOF2_xmax;
MCNUM mccFERMITOF2_ymin;
MCNUM mccFERMITOF2_ymax;
MCNUM mccFERMITOF2_zmin;
MCNUM mccFERMITOF2_zmax;
MCNUM mccFERMITOF2_bins;
MCNUM mccFERMITOF2_min;
MCNUM mccFERMITOF2_max;
MCNUM mccFERMITOF2_restore_neutron;
MCNUM mccFERMITOF2_radius;
char mccFERMITOF2_options[16384];
char mccFERMITOF2_filename[16384];
char mccFERMITOF2_geometry[16384];
char mccFERMITOF2_username1[16384];
char mccFERMITOF2_username2[16384];
char mccFERMITOF2_username3[16384];
int mccFERMITOF2_nowritefile;

/* Definition parameters for component 'PSD_SAMPLE' [37]. */
#define mccPSD_SAMPLE_nx 30
#define mccPSD_SAMPLE_ny 30
/* Setting parameters for component 'PSD_SAMPLE' [37]. */
char mccPSD_SAMPLE_filename[16384];
MCNUM mccPSD_SAMPLE_xmin;
MCNUM mccPSD_SAMPLE_xmax;
MCNUM mccPSD_SAMPLE_ymin;
MCNUM mccPSD_SAMPLE_ymax;
MCNUM mccPSD_SAMPLE_xwidth;
MCNUM mccPSD_SAMPLE_yheight;
MCNUM mccPSD_SAMPLE_restore_neutron;
int mccPSD_SAMPLE_nowritefile;

/* Definition parameters for component 'DivMon_Sample' [38]. */
#define mccDivMon_Sample_nh 30
#define mccDivMon_Sample_nv 30
/* Setting parameters for component 'DivMon_Sample' [38]. */
char mccDivMon_Sample_filename[16384];
MCNUM mccDivMon_Sample_xmin;
MCNUM mccDivMon_Sample_xmax;
MCNUM mccDivMon_Sample_ymin;
MCNUM mccDivMon_Sample_ymax;
MCNUM mccDivMon_Sample_xwidth;
MCNUM mccDivMon_Sample_yheight;
MCNUM mccDivMon_Sample_maxdiv_h;
MCNUM mccDivMon_Sample_maxdiv_v;
MCNUM mccDivMon_Sample_restore_neutron;
MCNUM mccDivMon_Sample_nx;
MCNUM mccDivMon_Sample_ny;
MCNUM mccDivMon_Sample_nz;
int mccDivMon_Sample_nowritefile;

/* Definition parameters for component 'EMON_SAMPLE' [39]. */
#define mccEMON_SAMPLE_nE 60
/* Setting parameters for component 'EMON_SAMPLE' [39]. */
char mccEMON_SAMPLE_filename[16384];
MCNUM mccEMON_SAMPLE_xmin;
MCNUM mccEMON_SAMPLE_xmax;
MCNUM mccEMON_SAMPLE_ymin;
MCNUM mccEMON_SAMPLE_ymax;
MCNUM mccEMON_SAMPLE_xwidth;
MCNUM mccEMON_SAMPLE_yheight;
MCNUM mccEMON_SAMPLE_Emin;
MCNUM mccEMON_SAMPLE_Emax;
MCNUM mccEMON_SAMPLE_restore_neutron;
int mccEMON_SAMPLE_nowritefile;

/* Setting parameters for component 'Sample' [41]. */
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

/* Definition parameters for component 'TOF_Det' [42]. */
#define mccTOF_Det_user1 FLT_MAX
#define mccTOF_Det_user2 FLT_MAX
#define mccTOF_Det_user3 FLT_MAX
/* Setting parameters for component 'TOF_Det' [42]. */
MCNUM mccTOF_Det_xwidth;
MCNUM mccTOF_Det_yheight;
MCNUM mccTOF_Det_zdepth;
MCNUM mccTOF_Det_xmin;
MCNUM mccTOF_Det_xmax;
MCNUM mccTOF_Det_ymin;
MCNUM mccTOF_Det_ymax;
MCNUM mccTOF_Det_zmin;
MCNUM mccTOF_Det_zmax;
MCNUM mccTOF_Det_bins;
MCNUM mccTOF_Det_min;
MCNUM mccTOF_Det_max;
MCNUM mccTOF_Det_restore_neutron;
MCNUM mccTOF_Det_radius;
char mccTOF_Det_options[16384];
char mccTOF_Det_filename[16384];
char mccTOF_Det_geometry[16384];
char mccTOF_Det_username1[16384];
char mccTOF_Det_username2[16384];
char mccTOF_Det_username3[16384];
int mccTOF_Det_nowritefile;

/* Definition parameters for component 'FoDet' [43]. */
#define mccFoDet_user1 FLT_MAX
#define mccFoDet_user2 FLT_MAX
#define mccFoDet_user3 FLT_MAX
/* Setting parameters for component 'FoDet' [43]. */
MCNUM mccFoDet_xwidth;
MCNUM mccFoDet_yheight;
MCNUM mccFoDet_zdepth;
MCNUM mccFoDet_xmin;
MCNUM mccFoDet_xmax;
MCNUM mccFoDet_ymin;
MCNUM mccFoDet_ymax;
MCNUM mccFoDet_zmin;
MCNUM mccFoDet_zmax;
MCNUM mccFoDet_bins;
MCNUM mccFoDet_min;
MCNUM mccFoDet_max;
MCNUM mccFoDet_restore_neutron;
MCNUM mccFoDet_radius;
char mccFoDet_options[16384];
char mccFoDet_filename[16384];
char mccFoDet_geometry[16384];
char mccFoDet_username1[16384];
char mccFoDet_username2[16384];
char mccFoDet_username3[16384];
int mccFoDet_nowritefile;

/* Definition parameters for component 'EMON_DET' [44]. */
#define mccEMON_DET_nE 80
/* Setting parameters for component 'EMON_DET' [44]. */
char mccEMON_DET_filename[16384];
MCNUM mccEMON_DET_xmin;
MCNUM mccEMON_DET_xmax;
MCNUM mccEMON_DET_ymin;
MCNUM mccEMON_DET_ymax;
MCNUM mccEMON_DET_xwidth;
MCNUM mccEMON_DET_yheight;
MCNUM mccEMON_DET_Emin;
MCNUM mccEMON_DET_Emax;
MCNUM mccEMON_DET_restore_neutron;
int mccEMON_DET_nowritefile;

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
#line 10924 "PSI_Focus.c"
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

/* User declarations for component 'csource' [2]. */
#define mccompcurname  csource
#define mccompcurtype  Source_gen
#define mccompcurindex 2
#define p_in mcccsource_p_in
#define lambda1 mcccsource_lambda1
#define lambda2 mcccsource_lambda2
#define lambda3 mcccsource_lambda3
#define pTable mcccsource_pTable
#define pTable_x mcccsource_pTable_x
#define pTable_y mcccsource_pTable_y
#define pTable_xmin mcccsource_pTable_xmin
#define pTable_xmax mcccsource_pTable_xmax
#define pTable_xsum mcccsource_pTable_xsum
#define pTable_ymin mcccsource_pTable_ymin
#define pTable_ymax mcccsource_pTable_ymax
#define pTable_ysum mcccsource_pTable_ysum
#define pTable_dxmin mcccsource_pTable_dxmin
#define pTable_dxmax mcccsource_pTable_dxmax
#define pTable_dymin mcccsource_pTable_dymin
#define pTable_dymax mcccsource_pTable_dymax
#define flux_file mcccsource_flux_file
#define xdiv_file mcccsource_xdiv_file
#define ydiv_file mcccsource_ydiv_file
#define radius mcccsource_radius
#define dist mcccsource_dist
#define focus_xw mcccsource_focus_xw
#define focus_yh mcccsource_focus_yh
#define focus_aw mcccsource_focus_aw
#define focus_ah mcccsource_focus_ah
#define E0 mcccsource_E0
#define dE mcccsource_dE
#define lambda0 mcccsource_lambda0
#define dlambda mcccsource_dlambda
#define I1 mcccsource_I1
#define yheight mcccsource_yheight
#define xwidth mcccsource_xwidth
#define verbose mcccsource_verbose
#define T1 mcccsource_T1
#define flux_file_perAA mcccsource_flux_file_perAA
#define flux_file_log mcccsource_flux_file_log
#define Lmin mcccsource_Lmin
#define Lmax mcccsource_Lmax
#define Emin mcccsource_Emin
#define Emax mcccsource_Emax
#define T2 mcccsource_T2
#define I2 mcccsource_I2
#define T3 mcccsource_T3
#define I3 mcccsource_I3
#define zdepth mcccsource_zdepth
#define target_index mcccsource_target_index
#line 184 "/usr/share/mcstas/2.5/sources/Source_gen.comp"

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

#line 11008 "PSI_Focus.c"
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

/* User declarations for component 'guide1' [3]. */
#define mccompcurname  guide1
#define mccompcurtype  Guide
#define mccompcurindex 3
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
#line 70 "/usr/share/mcstas/2.5/optics/Guide.comp"
t_Table pTable;
#line 11078 "PSI_Focus.c"
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

/* User declarations for component 'guide2' [4]. */
#define mccompcurname  guide2
#define mccompcurtype  Bender
#define mccompcurindex 4
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
#line 108 "/usr/share/mcstas/2.5/optics/Bender.comp"
double bk, mWin;
#line 11125 "PSI_Focus.c"
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

/* User declarations for component 'bunker' [5]. */
#define mccompcurname  bunker
#define mccompcurtype  Guide
#define mccompcurindex 5
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
#line 70 "/usr/share/mcstas/2.5/optics/Guide.comp"
t_Table pTable;
#line 11172 "PSI_Focus.c"
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

/* User declarations for component 'guide3' [6]. */
#define mccompcurname  guide3
#define mccompcurtype  Guide
#define mccompcurindex 6
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
#line 70 "/usr/share/mcstas/2.5/optics/Guide.comp"
t_Table pTable;
#line 11207 "PSI_Focus.c"
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

/* User declarations for component 'lambdaGuideExit' [7]. */
#define mccompcurname  lambdaGuideExit
#define mccompcurtype  L_monitor
#define mccompcurindex 7
#define nL mcclambdaGuideExit_nL
#define L_N mcclambdaGuideExit_L_N
#define L_p mcclambdaGuideExit_L_p
#define L_p2 mcclambdaGuideExit_L_p2
#define filename mcclambdaGuideExit_filename
#define xmin mcclambdaGuideExit_xmin
#define xmax mcclambdaGuideExit_xmax
#define ymin mcclambdaGuideExit_ymin
#define ymax mcclambdaGuideExit_ymax
#define xwidth mcclambdaGuideExit_xwidth
#define yheight mcclambdaGuideExit_yheight
#define Lmin mcclambdaGuideExit_Lmin
#define Lmax mcclambdaGuideExit_Lmax
#define restore_neutron mcclambdaGuideExit_restore_neutron
#define nowritefile mcclambdaGuideExit_nowritefile
#line 57 "/usr/share/mcstas/2.5/monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 11246 "PSI_Focus.c"
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

/* User declarations for component 'DivMonGuideExit' [8]. */
#define mccompcurname  DivMonGuideExit
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 8
#define nh mccDivMonGuideExit_nh
#define nv mccDivMonGuideExit_nv
#define Div_N mccDivMonGuideExit_Div_N
#define Div_p mccDivMonGuideExit_Div_p
#define Div_p2 mccDivMonGuideExit_Div_p2
#define filename mccDivMonGuideExit_filename
#define xmin mccDivMonGuideExit_xmin
#define xmax mccDivMonGuideExit_xmax
#define ymin mccDivMonGuideExit_ymin
#define ymax mccDivMonGuideExit_ymax
#define xwidth mccDivMonGuideExit_xwidth
#define yheight mccDivMonGuideExit_yheight
#define maxdiv_h mccDivMonGuideExit_maxdiv_h
#define maxdiv_v mccDivMonGuideExit_maxdiv_v
#define restore_neutron mccDivMonGuideExit_restore_neutron
#define nx mccDivMonGuideExit_nx
#define ny mccDivMonGuideExit_ny
#define nz mccDivMonGuideExit_nz
#define nowritefile mccDivMonGuideExit_nowritefile
#line 61 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
double Div_N[nh][nv];
double Div_p[nh][nv];
double Div_p2[nh][nv];
#line 11293 "PSI_Focus.c"
#undef nowritefile
#undef nz
#undef ny
#undef nx
#undef restore_neutron
#undef maxdiv_v
#undef maxdiv_h
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
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'PSDGuideExit' [9]. */
#define mccompcurname  PSDGuideExit
#define mccompcurtype  PSD_monitor
#define mccompcurindex 9
#define nx mccPSDGuideExit_nx
#define ny mccPSDGuideExit_ny
#define PSD_N mccPSDGuideExit_PSD_N
#define PSD_p mccPSDGuideExit_PSD_p
#define PSD_p2 mccPSDGuideExit_PSD_p2
#define filename mccPSDGuideExit_filename
#define xmin mccPSDGuideExit_xmin
#define xmax mccPSDGuideExit_xmax
#define ymin mccPSDGuideExit_ymin
#define ymax mccPSDGuideExit_ymax
#define xwidth mccPSDGuideExit_xwidth
#define yheight mccPSDGuideExit_yheight
#define restore_neutron mccPSDGuideExit_restore_neutron
#define nowritefile mccPSDGuideExit_nowritefile
#line 56 "/usr/share/mcstas/2.5/monitors/PSD_monitor.comp"
double PSD_N[nx][ny];
double PSD_p[nx][ny];
double PSD_p2[nx][ny];
#line 11339 "PSI_Focus.c"
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

/* User declarations for component 'FOCUSguide' [10]. */
#define mccompcurname  FOCUSguide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 10
#define w1c mccFOCUSguide_w1c
#define w2c mccFOCUSguide_w2c
#define ww mccFOCUSguide_ww
#define hh mccFOCUSguide_hh
#define whalf mccFOCUSguide_whalf
#define hhalf mccFOCUSguide_hhalf
#define lwhalf mccFOCUSguide_lwhalf
#define lhhalf mccFOCUSguide_lhhalf
#define w1 mccFOCUSguide_w1
#define h1 mccFOCUSguide_h1
#define w2 mccFOCUSguide_w2
#define h2 mccFOCUSguide_h2
#define l mccFOCUSguide_l
#define R0 mccFOCUSguide_R0
#define Qc mccFOCUSguide_Qc
#define alpha mccFOCUSguide_alpha
#define m mccFOCUSguide_m
#define nslit mccFOCUSguide_nslit
#define d mccFOCUSguide_d
#define Qcx mccFOCUSguide_Qcx
#define Qcy mccFOCUSguide_Qcy
#define alphax mccFOCUSguide_alphax
#define alphay mccFOCUSguide_alphay
#define W mccFOCUSguide_W
#define mx mccFOCUSguide_mx
#define my mccFOCUSguide_my
#define nu mccFOCUSguide_nu
#define phase mccFOCUSguide_phase
#line 81 "/usr/share/mcstas/2.5/optics/Guide_channeled.comp"
double w1c;
double w2c;
double ww, hh;
double whalf, hhalf;
double lwhalf, lhhalf;
#line 11396 "PSI_Focus.c"
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

/* User declarations for component 'FirstChopper' [11]. */
#define mccompcurname  FirstChopper
#define mccompcurtype  DiskChopper
#define mccompcurindex 11
#define Tg mccFirstChopper_Tg
#define To mccFirstChopper_To
#define delta_y mccFirstChopper_delta_y
#define height mccFirstChopper_height
#define omega mccFirstChopper_omega
#define theta_0 mccFirstChopper_theta_0
#define radius mccFirstChopper_radius
#define yheight mccFirstChopper_yheight
#define nu mccFirstChopper_nu
#define nslit mccFirstChopper_nslit
#define jitter mccFirstChopper_jitter
#define delay mccFirstChopper_delay
#define isfirst mccFirstChopper_isfirst
#define n_pulse mccFirstChopper_n_pulse
#define abs_out mccFirstChopper_abs_out
#define phase mccFirstChopper_phase
#define xwidth mccFirstChopper_xwidth
#define verbose mccFirstChopper_verbose
#line 66 "/usr/share/mcstas/2.5/optics/DiskChopper.comp"
double Tg,To,delta_y,height,omega;
#line 11453 "PSI_Focus.c"
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

/* User declarations for component 'DISCTOF' [12]. */
#define mccompcurname  DISCTOF
#define mccompcurtype  Monitor_nD
#define mccompcurindex 12
#define user1 mccDISCTOF_user1
#define user2 mccDISCTOF_user2
#define user3 mccDISCTOF_user3
#define DEFS mccDISCTOF_DEFS
#define Vars mccDISCTOF_Vars
#define detector mccDISCTOF_detector
#define offdata mccDISCTOF_offdata
#define xwidth mccDISCTOF_xwidth
#define yheight mccDISCTOF_yheight
#define zdepth mccDISCTOF_zdepth
#define xmin mccDISCTOF_xmin
#define xmax mccDISCTOF_xmax
#define ymin mccDISCTOF_ymin
#define ymax mccDISCTOF_ymax
#define zmin mccDISCTOF_zmin
#define zmax mccDISCTOF_zmax
#define bins mccDISCTOF_bins
#define min mccDISCTOF_min
#define max mccDISCTOF_max
#define restore_neutron mccDISCTOF_restore_neutron
#define radius mccDISCTOF_radius
#define options mccDISCTOF_options
#define filename mccDISCTOF_filename
#define geometry mccDISCTOF_geometry
#define username1 mccDISCTOF_username1
#define username2 mccDISCTOF_username2
#define username3 mccDISCTOF_username3
#define nowritefile mccDISCTOF_nowritefile
#line 222 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
  MonitornD_Defines_type DEFS;
  MonitornD_Variables_type Vars;
  MCDETECTOR detector;
  off_struct offdata;
#line 11513 "PSI_Focus.c"
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

/* User declarations for component 'PSDmon1Chopper' [13]. */
#define mccompcurname  PSDmon1Chopper
#define mccompcurtype  PSD_monitor
#define mccompcurindex 13
#define nx mccPSDmon1Chopper_nx
#define ny mccPSDmon1Chopper_ny
#define PSD_N mccPSDmon1Chopper_PSD_N
#define PSD_p mccPSDmon1Chopper_PSD_p
#define PSD_p2 mccPSDmon1Chopper_PSD_p2
#define filename mccPSDmon1Chopper_filename
#define xmin mccPSDmon1Chopper_xmin
#define xmax mccPSDmon1Chopper_xmax
#define ymin mccPSDmon1Chopper_ymin
#define ymax mccPSDmon1Chopper_ymax
#define xwidth mccPSDmon1Chopper_xwidth
#define yheight mccPSDmon1Chopper_yheight
#define restore_neutron mccPSDmon1Chopper_restore_neutron
#define nowritefile mccPSDmon1Chopper_nowritefile
#line 56 "/usr/share/mcstas/2.5/monitors/PSD_monitor.comp"
double PSD_N[nx][ny];
double PSD_p[nx][ny];
double PSD_p2[nx][ny];
#line 11568 "PSI_Focus.c"
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

/* User declarations for component 'VacuumTube1entry' [14]. */
#define mccompcurname  VacuumTube1entry
#define mccompcurtype  Slit
#define mccompcurindex 14
#define xmin mccVacuumTube1entry_xmin
#define xmax mccVacuumTube1entry_xmax
#define ymin mccVacuumTube1entry_ymin
#define ymax mccVacuumTube1entry_ymax
#define radius mccVacuumTube1entry_radius
#define xwidth mccVacuumTube1entry_xwidth
#define yheight mccVacuumTube1entry_yheight
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

/* User declarations for component 'VacuumTube1exit' [15]. */
#define mccompcurname  VacuumTube1exit
#define mccompcurtype  Slit
#define mccompcurindex 15
#define xmin mccVacuumTube1exit_xmin
#define xmax mccVacuumTube1exit_xmax
#define ymin mccVacuumTube1exit_ymin
#define ymax mccVacuumTube1exit_ymax
#define radius mccVacuumTube1exit_radius
#define xwidth mccVacuumTube1exit_xwidth
#define yheight mccVacuumTube1exit_yheight
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

/* User declarations for component 'VacuumTube2entry' [16]. */
#define mccompcurname  VacuumTube2entry
#define mccompcurtype  Slit
#define mccompcurindex 16
#define xmin mccVacuumTube2entry_xmin
#define xmax mccVacuumTube2entry_xmax
#define ymin mccVacuumTube2entry_ymin
#define ymax mccVacuumTube2entry_ymax
#define radius mccVacuumTube2entry_radius
#define xwidth mccVacuumTube2entry_xwidth
#define yheight mccVacuumTube2entry_yheight
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

/* User declarations for component 'VacuumTube2exit' [17]. */
#define mccompcurname  VacuumTube2exit
#define mccompcurtype  Slit
#define mccompcurindex 17
#define xmin mccVacuumTube2exit_xmin
#define xmax mccVacuumTube2exit_xmax
#define ymin mccVacuumTube2exit_ymin
#define ymax mccVacuumTube2exit_ymax
#define radius mccVacuumTube2exit_radius
#define xwidth mccVacuumTube2exit_xwidth
#define yheight mccVacuumTube2exit_yheight
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

/* User declarations for component 'VacuumTube3entry' [18]. */
#define mccompcurname  VacuumTube3entry
#define mccompcurtype  Slit
#define mccompcurindex 18
#define xmin mccVacuumTube3entry_xmin
#define xmax mccVacuumTube3entry_xmax
#define ymin mccVacuumTube3entry_ymin
#define ymax mccVacuumTube3entry_ymax
#define radius mccVacuumTube3entry_radius
#define xwidth mccVacuumTube3entry_xwidth
#define yheight mccVacuumTube3entry_yheight
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

/* User declarations for component 'VacuumTube3exit' [19]. */
#define mccompcurname  VacuumTube3exit
#define mccompcurtype  Slit
#define mccompcurindex 19
#define xmin mccVacuumTube3exit_xmin
#define xmax mccVacuumTube3exit_xmax
#define ymin mccVacuumTube3exit_ymin
#define ymax mccVacuumTube3exit_ymax
#define radius mccVacuumTube3exit_radius
#define xwidth mccVacuumTube3exit_xwidth
#define yheight mccVacuumTube3exit_yheight
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

/* User declarations for component 'PSDmonMono' [20]. */
#define mccompcurname  PSDmonMono
#define mccompcurtype  PSD_monitor
#define mccompcurindex 20
#define nx mccPSDmonMono_nx
#define ny mccPSDmonMono_ny
#define PSD_N mccPSDmonMono_PSD_N
#define PSD_p mccPSDmonMono_PSD_p
#define PSD_p2 mccPSDmonMono_PSD_p2
#define filename mccPSDmonMono_filename
#define xmin mccPSDmonMono_xmin
#define xmax mccPSDmonMono_xmax
#define ymin mccPSDmonMono_ymin
#define ymax mccPSDmonMono_ymax
#define xwidth mccPSDmonMono_xwidth
#define yheight mccPSDmonMono_yheight
#define restore_neutron mccPSDmonMono_restore_neutron
#define nowritefile mccPSDmonMono_nowritefile
#line 56 "/usr/share/mcstas/2.5/monitors/PSD_monitor.comp"
double PSD_N[nx][ny];
double PSD_p[nx][ny];
double PSD_p2[nx][ny];
#line 11741 "PSI_Focus.c"
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

/* User declarations for component 'MONOTOF' [21]. */
#define mccompcurname  MONOTOF
#define mccompcurtype  Monitor_nD
#define mccompcurindex 21
#define user1 mccMONOTOF_user1
#define user2 mccMONOTOF_user2
#define user3 mccMONOTOF_user3
#define DEFS mccMONOTOF_DEFS
#define Vars mccMONOTOF_Vars
#define detector mccMONOTOF_detector
#define offdata mccMONOTOF_offdata
#define xwidth mccMONOTOF_xwidth
#define yheight mccMONOTOF_yheight
#define zdepth mccMONOTOF_zdepth
#define xmin mccMONOTOF_xmin
#define xmax mccMONOTOF_xmax
#define ymin mccMONOTOF_ymin
#define ymax mccMONOTOF_ymax
#define zmin mccMONOTOF_zmin
#define zmax mccMONOTOF_zmax
#define bins mccMONOTOF_bins
#define min mccMONOTOF_min
#define max mccMONOTOF_max
#define restore_neutron mccMONOTOF_restore_neutron
#define radius mccMONOTOF_radius
#define options mccMONOTOF_options
#define filename mccMONOTOF_filename
#define geometry mccMONOTOF_geometry
#define username1 mccMONOTOF_username1
#define username2 mccMONOTOF_username2
#define username3 mccMONOTOF_username3
#define nowritefile mccMONOTOF_nowritefile
#line 222 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
  MonitornD_Defines_type DEFS;
  MonitornD_Variables_type Vars;
  MCDETECTOR detector;
  off_struct offdata;
#line 11797 "PSI_Focus.c"
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

/* User declarations for component 'DivMonMono' [22]. */
#define mccompcurname  DivMonMono
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 22
#define nh mccDivMonMono_nh
#define nv mccDivMonMono_nv
#define Div_N mccDivMonMono_Div_N
#define Div_p mccDivMonMono_Div_p
#define Div_p2 mccDivMonMono_Div_p2
#define filename mccDivMonMono_filename
#define xmin mccDivMonMono_xmin
#define xmax mccDivMonMono_xmax
#define ymin mccDivMonMono_ymin
#define ymax mccDivMonMono_ymax
#define xwidth mccDivMonMono_xwidth
#define yheight mccDivMonMono_yheight
#define maxdiv_h mccDivMonMono_maxdiv_h
#define maxdiv_v mccDivMonMono_maxdiv_v
#define restore_neutron mccDivMonMono_restore_neutron
#define nx mccDivMonMono_nx
#define ny mccDivMonMono_ny
#define nz mccDivMonMono_nz
#define nowritefile mccDivMonMono_nowritefile
#line 61 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
double Div_N[nh][nv];
double Div_p[nh][nv];
double Div_p2[nh][nv];
#line 11857 "PSI_Focus.c"
#undef nowritefile
#undef nz
#undef ny
#undef nx
#undef restore_neutron
#undef maxdiv_v
#undef maxdiv_h
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
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'focus_mono' [23]. */
#define mccompcurname  focus_mono
#define mccompcurtype  Arm
#define mccompcurindex 23
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'mono' [24]. */
#define mccompcurname  mono
#define mccompcurtype  Monochromator_2foc
#define mccompcurindex 24
#define mos_y mccmono_mos_y
#define mos_z mccmono_mos_z
#define mono_Q mccmono_mono_Q
#define SlabWidth mccmono_SlabWidth
#define SlabHeight mccmono_SlabHeight
#define rTable mccmono_rTable
#define reflect mccmono_reflect
#define zwidth mccmono_zwidth
#define yheight mccmono_yheight
#define gap mccmono_gap
#define NH mccmono_NH
#define NV mccmono_NV
#define mosaich mccmono_mosaich
#define mosaicv mccmono_mosaicv
#define r0 mccmono_r0
#define Q mccmono_Q
#define RV mccmono_RV
#define RH mccmono_RH
#define DM mccmono_DM
#define mosaic mccmono_mosaic
#define width mccmono_width
#define height mccmono_height
#define verbose mccmono_verbose
#line 88 "/usr/share/mcstas/2.5/contrib/Monochromator_2foc.comp"
#ifndef DIV_CUTOFF
#define DIV_CUTOFF 2            /* ~ 10^-5 cutoff. */
#endif
double mos_y; /* mosaic, in radians */
double mos_z;
double mono_Q;
double SlabWidth, SlabHeight;
t_Table rTable;
#line 11925 "PSI_Focus.c"
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

/* User declarations for component 'a2' [25]. */
#define mccompcurname  a2
#define mccompcurtype  Arm
#define mccompcurindex 25
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'FERMITOF_before' [26]. */
#define mccompcurname  FERMITOF_before
#define mccompcurtype  Monitor_nD
#define mccompcurindex 26
#define user1 mccFERMITOF_before_user1
#define user2 mccFERMITOF_before_user2
#define user3 mccFERMITOF_before_user3
#define DEFS mccFERMITOF_before_DEFS
#define Vars mccFERMITOF_before_Vars
#define detector mccFERMITOF_before_detector
#define offdata mccFERMITOF_before_offdata
#define xwidth mccFERMITOF_before_xwidth
#define yheight mccFERMITOF_before_yheight
#define zdepth mccFERMITOF_before_zdepth
#define xmin mccFERMITOF_before_xmin
#define xmax mccFERMITOF_before_xmax
#define ymin mccFERMITOF_before_ymin
#define ymax mccFERMITOF_before_ymax
#define zmin mccFERMITOF_before_zmin
#define zmax mccFERMITOF_before_zmax
#define bins mccFERMITOF_before_bins
#define min mccFERMITOF_before_min
#define max mccFERMITOF_before_max
#define restore_neutron mccFERMITOF_before_restore_neutron
#define radius mccFERMITOF_before_radius
#define options mccFERMITOF_before_options
#define filename mccFERMITOF_before_filename
#define geometry mccFERMITOF_before_geometry
#define username1 mccFERMITOF_before_username1
#define username2 mccFERMITOF_before_username2
#define username3 mccFERMITOF_before_username3
#define nowritefile mccFERMITOF_before_nowritefile
#line 222 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
  MonitornD_Defines_type DEFS;
  MonitornD_Variables_type Vars;
  MCDETECTOR detector;
  off_struct offdata;
#line 11998 "PSI_Focus.c"
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

/* User declarations for component 'lambdaFermi' [27]. */
#define mccompcurname  lambdaFermi
#define mccompcurtype  L_monitor
#define mccompcurindex 27
#define nL mcclambdaFermi_nL
#define L_N mcclambdaFermi_L_N
#define L_p mcclambdaFermi_L_p
#define L_p2 mcclambdaFermi_L_p2
#define filename mcclambdaFermi_filename
#define xmin mcclambdaFermi_xmin
#define xmax mcclambdaFermi_xmax
#define ymin mcclambdaFermi_ymin
#define ymax mcclambdaFermi_ymax
#define xwidth mcclambdaFermi_xwidth
#define yheight mcclambdaFermi_yheight
#define Lmin mcclambdaFermi_Lmin
#define Lmax mcclambdaFermi_Lmax
#define restore_neutron mcclambdaFermi_restore_neutron
#define nowritefile mcclambdaFermi_nowritefile
#line 57 "/usr/share/mcstas/2.5/monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 12053 "PSI_Focus.c"
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

/* User declarations for component 'EMON_Fermi' [28]. */
#define mccompcurname  EMON_Fermi
#define mccompcurtype  E_monitor
#define mccompcurindex 28
#define nE mccEMON_Fermi_nE
#define E_N mccEMON_Fermi_E_N
#define E_p mccEMON_Fermi_E_p
#define E_p2 mccEMON_Fermi_E_p2
#define S_p mccEMON_Fermi_S_p
#define S_pE mccEMON_Fermi_S_pE
#define S_pE2 mccEMON_Fermi_S_pE2
#define filename mccEMON_Fermi_filename
#define xmin mccEMON_Fermi_xmin
#define xmax mccEMON_Fermi_xmax
#define ymin mccEMON_Fermi_ymin
#define ymax mccEMON_Fermi_ymax
#define xwidth mccEMON_Fermi_xwidth
#define yheight mccEMON_Fermi_yheight
#define Emin mccEMON_Fermi_Emin
#define Emax mccEMON_Fermi_Emax
#define restore_neutron mccEMON_Fermi_restore_neutron
#define nowritefile mccEMON_Fermi_nowritefile
#line 60 "/usr/share/mcstas/2.5/monitors/E_monitor.comp"
double E_N[nE];
double E_p[nE], E_p2[nE];
double S_p, S_pE, S_pE2;
#line 12099 "PSI_Focus.c"
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

/* User declarations for component 'DivMonfermi1' [29]. */
#define mccompcurname  DivMonfermi1
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 29
#define nh mccDivMonfermi1_nh
#define nv mccDivMonfermi1_nv
#define Div_N mccDivMonfermi1_Div_N
#define Div_p mccDivMonfermi1_Div_p
#define Div_p2 mccDivMonfermi1_Div_p2
#define filename mccDivMonfermi1_filename
#define xmin mccDivMonfermi1_xmin
#define xmax mccDivMonfermi1_xmax
#define ymin mccDivMonfermi1_ymin
#define ymax mccDivMonfermi1_ymax
#define xwidth mccDivMonfermi1_xwidth
#define yheight mccDivMonfermi1_yheight
#define maxdiv_h mccDivMonfermi1_maxdiv_h
#define maxdiv_v mccDivMonfermi1_maxdiv_v
#define restore_neutron mccDivMonfermi1_restore_neutron
#define nx mccDivMonfermi1_nx
#define ny mccDivMonfermi1_ny
#define nz mccDivMonfermi1_nz
#define nowritefile mccDivMonfermi1_nowritefile
#line 61 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
double Div_N[nh][nv];
double Div_p[nh][nv];
double Div_p2[nh][nv];
#line 12149 "PSI_Focus.c"
#undef nowritefile
#undef nz
#undef ny
#undef nx
#undef restore_neutron
#undef maxdiv_v
#undef maxdiv_h
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
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'PSD_Fermi1' [30]. */
#define mccompcurname  PSD_Fermi1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 30
#define nx mccPSD_Fermi1_nx
#define ny mccPSD_Fermi1_ny
#define PSD_N mccPSD_Fermi1_PSD_N
#define PSD_p mccPSD_Fermi1_PSD_p
#define PSD_p2 mccPSD_Fermi1_PSD_p2
#define filename mccPSD_Fermi1_filename
#define xmin mccPSD_Fermi1_xmin
#define xmax mccPSD_Fermi1_xmax
#define ymin mccPSD_Fermi1_ymin
#define ymax mccPSD_Fermi1_ymax
#define xwidth mccPSD_Fermi1_xwidth
#define yheight mccPSD_Fermi1_yheight
#define restore_neutron mccPSD_Fermi1_restore_neutron
#define nowritefile mccPSD_Fermi1_nowritefile
#line 56 "/usr/share/mcstas/2.5/monitors/PSD_monitor.comp"
double PSD_N[nx][ny];
double PSD_p[nx][ny];
double PSD_p2[nx][ny];
#line 12195 "PSI_Focus.c"
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

/* User declarations for component 'FoChopper' [31]. */
#define mccompcurname  FoChopper
#define mccompcurtype  FermiChopper
#define mccompcurindex 31
#define FCVars mccFoChopper_FCVars
#define phase mccFoChopper_phase
#define radius mccFoChopper_radius
#define nu mccFoChopper_nu
#define w mccFoChopper_w
#define nslit mccFoChopper_nslit
#define R0 mccFoChopper_R0
#define Qc mccFoChopper_Qc
#define alpha mccFoChopper_alpha
#define m mccFoChopper_m
#define W mccFoChopper_W
#define length mccFoChopper_length
#define eff mccFoChopper_eff
#define zero_time mccFoChopper_zero_time
#define xwidth mccFoChopper_xwidth
#define verbose mccFoChopper_verbose
#define yheight mccFoChopper_yheight
#define curvature mccFoChopper_curvature
#define delay mccFoChopper_delay
#line 278 "/usr/share/mcstas/2.5/optics/FermiChopper.comp"
        struct FermiChopper_struct FCVars;
#line 12239 "PSI_Focus.c"
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

/* User declarations for component 'PSD_Fermi2' [32]. */
#define mccompcurname  PSD_Fermi2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 32
#define nx mccPSD_Fermi2_nx
#define ny mccPSD_Fermi2_ny
#define PSD_N mccPSD_Fermi2_PSD_N
#define PSD_p mccPSD_Fermi2_PSD_p
#define PSD_p2 mccPSD_Fermi2_PSD_p2
#define filename mccPSD_Fermi2_filename
#define xmin mccPSD_Fermi2_xmin
#define xmax mccPSD_Fermi2_xmax
#define ymin mccPSD_Fermi2_ymin
#define ymax mccPSD_Fermi2_ymax
#define xwidth mccPSD_Fermi2_xwidth
#define yheight mccPSD_Fermi2_yheight
#define restore_neutron mccPSD_Fermi2_restore_neutron
#define nowritefile mccPSD_Fermi2_nowritefile
#line 56 "/usr/share/mcstas/2.5/monitors/PSD_monitor.comp"
double PSD_N[nx][ny];
double PSD_p[nx][ny];
double PSD_p2[nx][ny];
#line 12285 "PSI_Focus.c"
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

/* User declarations for component 'DivMonfermi2' [33]. */
#define mccompcurname  DivMonfermi2
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 33
#define nh mccDivMonfermi2_nh
#define nv mccDivMonfermi2_nv
#define Div_N mccDivMonfermi2_Div_N
#define Div_p mccDivMonfermi2_Div_p
#define Div_p2 mccDivMonfermi2_Div_p2
#define filename mccDivMonfermi2_filename
#define xmin mccDivMonfermi2_xmin
#define xmax mccDivMonfermi2_xmax
#define ymin mccDivMonfermi2_ymin
#define ymax mccDivMonfermi2_ymax
#define xwidth mccDivMonfermi2_xwidth
#define yheight mccDivMonfermi2_yheight
#define maxdiv_h mccDivMonfermi2_maxdiv_h
#define maxdiv_v mccDivMonfermi2_maxdiv_v
#define restore_neutron mccDivMonfermi2_restore_neutron
#define nx mccDivMonfermi2_nx
#define ny mccDivMonfermi2_ny
#define nz mccDivMonfermi2_nz
#define nowritefile mccDivMonfermi2_nowritefile
#line 61 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
double Div_N[nh][nv];
double Div_p[nh][nv];
double Div_p2[nh][nv];
#line 12331 "PSI_Focus.c"
#undef nowritefile
#undef nz
#undef ny
#undef nx
#undef restore_neutron
#undef maxdiv_v
#undef maxdiv_h
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
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'FERMITOF1' [34]. */
#define mccompcurname  FERMITOF1
#define mccompcurtype  Monitor_nD
#define mccompcurindex 34
#define user1 mccFERMITOF1_user1
#define user2 mccFERMITOF1_user2
#define user3 mccFERMITOF1_user3
#define DEFS mccFERMITOF1_DEFS
#define Vars mccFERMITOF1_Vars
#define detector mccFERMITOF1_detector
#define offdata mccFERMITOF1_offdata
#define xwidth mccFERMITOF1_xwidth
#define yheight mccFERMITOF1_yheight
#define zdepth mccFERMITOF1_zdepth
#define xmin mccFERMITOF1_xmin
#define xmax mccFERMITOF1_xmax
#define ymin mccFERMITOF1_ymin
#define ymax mccFERMITOF1_ymax
#define zmin mccFERMITOF1_zmin
#define zmax mccFERMITOF1_zmax
#define bins mccFERMITOF1_bins
#define min mccFERMITOF1_min
#define max mccFERMITOF1_max
#define restore_neutron mccFERMITOF1_restore_neutron
#define radius mccFERMITOF1_radius
#define options mccFERMITOF1_options
#define filename mccFERMITOF1_filename
#define geometry mccFERMITOF1_geometry
#define username1 mccFERMITOF1_username1
#define username2 mccFERMITOF1_username2
#define username3 mccFERMITOF1_username3
#define nowritefile mccFERMITOF1_nowritefile
#line 222 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
  MonitornD_Defines_type DEFS;
  MonitornD_Variables_type Vars;
  MCDETECTOR detector;
  off_struct offdata;
#line 12392 "PSI_Focus.c"
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

/* User declarations for component 'SAMPLE_SLIT' [35]. */
#define mccompcurname  SAMPLE_SLIT
#define mccompcurtype  Slit
#define mccompcurindex 35
#define xmin mccSAMPLE_SLIT_xmin
#define xmax mccSAMPLE_SLIT_xmax
#define ymin mccSAMPLE_SLIT_ymin
#define ymax mccSAMPLE_SLIT_ymax
#define radius mccSAMPLE_SLIT_radius
#define xwidth mccSAMPLE_SLIT_xwidth
#define yheight mccSAMPLE_SLIT_yheight
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

/* User declarations for component 'FERMITOF2' [36]. */
#define mccompcurname  FERMITOF2
#define mccompcurtype  Monitor_nD
#define mccompcurindex 36
#define user1 mccFERMITOF2_user1
#define user2 mccFERMITOF2_user2
#define user3 mccFERMITOF2_user3
#define DEFS mccFERMITOF2_DEFS
#define Vars mccFERMITOF2_Vars
#define detector mccFERMITOF2_detector
#define offdata mccFERMITOF2_offdata
#define xwidth mccFERMITOF2_xwidth
#define yheight mccFERMITOF2_yheight
#define zdepth mccFERMITOF2_zdepth
#define xmin mccFERMITOF2_xmin
#define xmax mccFERMITOF2_xmax
#define ymin mccFERMITOF2_ymin
#define ymax mccFERMITOF2_ymax
#define zmin mccFERMITOF2_zmin
#define zmax mccFERMITOF2_zmax
#define bins mccFERMITOF2_bins
#define min mccFERMITOF2_min
#define max mccFERMITOF2_max
#define restore_neutron mccFERMITOF2_restore_neutron
#define radius mccFERMITOF2_radius
#define options mccFERMITOF2_options
#define filename mccFERMITOF2_filename
#define geometry mccFERMITOF2_geometry
#define username1 mccFERMITOF2_username1
#define username2 mccFERMITOF2_username2
#define username3 mccFERMITOF2_username3
#define nowritefile mccFERMITOF2_nowritefile
#line 222 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
  MonitornD_Defines_type DEFS;
  MonitornD_Variables_type Vars;
  MCDETECTOR detector;
  off_struct offdata;
#line 12484 "PSI_Focus.c"
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

/* User declarations for component 'PSD_SAMPLE' [37]. */
#define mccompcurname  PSD_SAMPLE
#define mccompcurtype  PSD_monitor
#define mccompcurindex 37
#define nx mccPSD_SAMPLE_nx
#define ny mccPSD_SAMPLE_ny
#define PSD_N mccPSD_SAMPLE_PSD_N
#define PSD_p mccPSD_SAMPLE_PSD_p
#define PSD_p2 mccPSD_SAMPLE_PSD_p2
#define filename mccPSD_SAMPLE_filename
#define xmin mccPSD_SAMPLE_xmin
#define xmax mccPSD_SAMPLE_xmax
#define ymin mccPSD_SAMPLE_ymin
#define ymax mccPSD_SAMPLE_ymax
#define xwidth mccPSD_SAMPLE_xwidth
#define yheight mccPSD_SAMPLE_yheight
#define restore_neutron mccPSD_SAMPLE_restore_neutron
#define nowritefile mccPSD_SAMPLE_nowritefile
#line 56 "/usr/share/mcstas/2.5/monitors/PSD_monitor.comp"
double PSD_N[nx][ny];
double PSD_p[nx][ny];
double PSD_p2[nx][ny];
#line 12539 "PSI_Focus.c"
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

/* User declarations for component 'DivMon_Sample' [38]. */
#define mccompcurname  DivMon_Sample
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 38
#define nh mccDivMon_Sample_nh
#define nv mccDivMon_Sample_nv
#define Div_N mccDivMon_Sample_Div_N
#define Div_p mccDivMon_Sample_Div_p
#define Div_p2 mccDivMon_Sample_Div_p2
#define filename mccDivMon_Sample_filename
#define xmin mccDivMon_Sample_xmin
#define xmax mccDivMon_Sample_xmax
#define ymin mccDivMon_Sample_ymin
#define ymax mccDivMon_Sample_ymax
#define xwidth mccDivMon_Sample_xwidth
#define yheight mccDivMon_Sample_yheight
#define maxdiv_h mccDivMon_Sample_maxdiv_h
#define maxdiv_v mccDivMon_Sample_maxdiv_v
#define restore_neutron mccDivMon_Sample_restore_neutron
#define nx mccDivMon_Sample_nx
#define ny mccDivMon_Sample_ny
#define nz mccDivMon_Sample_nz
#define nowritefile mccDivMon_Sample_nowritefile
#line 61 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
double Div_N[nh][nv];
double Div_p[nh][nv];
double Div_p2[nh][nv];
#line 12585 "PSI_Focus.c"
#undef nowritefile
#undef nz
#undef ny
#undef nx
#undef restore_neutron
#undef maxdiv_v
#undef maxdiv_h
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
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'EMON_SAMPLE' [39]. */
#define mccompcurname  EMON_SAMPLE
#define mccompcurtype  E_monitor
#define mccompcurindex 39
#define nE mccEMON_SAMPLE_nE
#define E_N mccEMON_SAMPLE_E_N
#define E_p mccEMON_SAMPLE_E_p
#define E_p2 mccEMON_SAMPLE_E_p2
#define S_p mccEMON_SAMPLE_S_p
#define S_pE mccEMON_SAMPLE_S_pE
#define S_pE2 mccEMON_SAMPLE_S_pE2
#define filename mccEMON_SAMPLE_filename
#define xmin mccEMON_SAMPLE_xmin
#define xmax mccEMON_SAMPLE_xmax
#define ymin mccEMON_SAMPLE_ymin
#define ymax mccEMON_SAMPLE_ymax
#define xwidth mccEMON_SAMPLE_xwidth
#define yheight mccEMON_SAMPLE_yheight
#define Emin mccEMON_SAMPLE_Emin
#define Emax mccEMON_SAMPLE_Emax
#define restore_neutron mccEMON_SAMPLE_restore_neutron
#define nowritefile mccEMON_SAMPLE_nowritefile
#line 60 "/usr/share/mcstas/2.5/monitors/E_monitor.comp"
double E_N[nE];
double E_p[nE], E_p2[nE];
double S_p, S_pE, S_pE2;
#line 12635 "PSI_Focus.c"
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

/* User declarations for component 'a3' [40]. */
#define mccompcurname  a3
#define mccompcurtype  Arm
#define mccompcurindex 40
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Sample' [41]. */
#define mccompcurname  Sample
#define mccompcurtype  V_sample
#define mccompcurindex 41
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
#line 12703 "PSI_Focus.c"
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

/* User declarations for component 'TOF_Det' [42]. */
#define mccompcurname  TOF_Det
#define mccompcurtype  Monitor_nD
#define mccompcurindex 42
#define user1 mccTOF_Det_user1
#define user2 mccTOF_Det_user2
#define user3 mccTOF_Det_user3
#define DEFS mccTOF_Det_DEFS
#define Vars mccTOF_Det_Vars
#define detector mccTOF_Det_detector
#define offdata mccTOF_Det_offdata
#define xwidth mccTOF_Det_xwidth
#define yheight mccTOF_Det_yheight
#define zdepth mccTOF_Det_zdepth
#define xmin mccTOF_Det_xmin
#define xmax mccTOF_Det_xmax
#define ymin mccTOF_Det_ymin
#define ymax mccTOF_Det_ymax
#define zmin mccTOF_Det_zmin
#define zmax mccTOF_Det_zmax
#define bins mccTOF_Det_bins
#define min mccTOF_Det_min
#define max mccTOF_Det_max
#define restore_neutron mccTOF_Det_restore_neutron
#define radius mccTOF_Det_radius
#define options mccTOF_Det_options
#define filename mccTOF_Det_filename
#define geometry mccTOF_Det_geometry
#define username1 mccTOF_Det_username1
#define username2 mccTOF_Det_username2
#define username3 mccTOF_Det_username3
#define nowritefile mccTOF_Det_nowritefile
#line 222 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
  MonitornD_Defines_type DEFS;
  MonitornD_Variables_type Vars;
  MCDETECTOR detector;
  off_struct offdata;
#line 12776 "PSI_Focus.c"
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

/* User declarations for component 'FoDet' [43]. */
#define mccompcurname  FoDet
#define mccompcurtype  Monitor_nD
#define mccompcurindex 43
#define user1 mccFoDet_user1
#define user2 mccFoDet_user2
#define user3 mccFoDet_user3
#define DEFS mccFoDet_DEFS
#define Vars mccFoDet_Vars
#define detector mccFoDet_detector
#define offdata mccFoDet_offdata
#define xwidth mccFoDet_xwidth
#define yheight mccFoDet_yheight
#define zdepth mccFoDet_zdepth
#define xmin mccFoDet_xmin
#define xmax mccFoDet_xmax
#define ymin mccFoDet_ymin
#define ymax mccFoDet_ymax
#define zmin mccFoDet_zmin
#define zmax mccFoDet_zmax
#define bins mccFoDet_bins
#define min mccFoDet_min
#define max mccFoDet_max
#define restore_neutron mccFoDet_restore_neutron
#define radius mccFoDet_radius
#define options mccFoDet_options
#define filename mccFoDet_filename
#define geometry mccFoDet_geometry
#define username1 mccFoDet_username1
#define username2 mccFoDet_username2
#define username3 mccFoDet_username3
#define nowritefile mccFoDet_nowritefile
#line 222 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
  MonitornD_Defines_type DEFS;
  MonitornD_Variables_type Vars;
  MCDETECTOR detector;
  off_struct offdata;
#line 12846 "PSI_Focus.c"
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

/* User declarations for component 'EMON_DET' [44]. */
#define mccompcurname  EMON_DET
#define mccompcurtype  E_monitor
#define mccompcurindex 44
#define nE mccEMON_DET_nE
#define E_N mccEMON_DET_E_N
#define E_p mccEMON_DET_E_p
#define E_p2 mccEMON_DET_E_p2
#define S_p mccEMON_DET_S_p
#define S_pE mccEMON_DET_S_pE
#define S_pE2 mccEMON_DET_S_pE2
#define filename mccEMON_DET_filename
#define xmin mccEMON_DET_xmin
#define xmax mccEMON_DET_xmax
#define ymin mccEMON_DET_ymin
#define ymax mccEMON_DET_ymax
#define xwidth mccEMON_DET_xwidth
#define yheight mccEMON_DET_yheight
#define Emin mccEMON_DET_Emin
#define Emax mccEMON_DET_Emax
#define restore_neutron mccEMON_DET_restore_neutron
#define nowritefile mccEMON_DET_nowritefile
#line 60 "/usr/share/mcstas/2.5/monitors/E_monitor.comp"
double E_N[nE];
double E_p[nE], E_p2[nE];
double S_p, S_pE, S_pE2;
#line 12905 "PSI_Focus.c"
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
Coords mcposacsource, mcposrcsource;
Rotation mcrotacsource, mcrotrcsource;
Coords mcposaguide1, mcposrguide1;
Rotation mcrotaguide1, mcrotrguide1;
Coords mcposaguide2, mcposrguide2;
Rotation mcrotaguide2, mcrotrguide2;
Coords mcposabunker, mcposrbunker;
Rotation mcrotabunker, mcrotrbunker;
Coords mcposaguide3, mcposrguide3;
Rotation mcrotaguide3, mcrotrguide3;
Coords mcposalambdaGuideExit, mcposrlambdaGuideExit;
Rotation mcrotalambdaGuideExit, mcrotrlambdaGuideExit;
Coords mcposaDivMonGuideExit, mcposrDivMonGuideExit;
Rotation mcrotaDivMonGuideExit, mcrotrDivMonGuideExit;
Coords mcposaPSDGuideExit, mcposrPSDGuideExit;
Rotation mcrotaPSDGuideExit, mcrotrPSDGuideExit;
Coords mcposaFOCUSguide, mcposrFOCUSguide;
Rotation mcrotaFOCUSguide, mcrotrFOCUSguide;
Coords mcposaFirstChopper, mcposrFirstChopper;
Rotation mcrotaFirstChopper, mcrotrFirstChopper;
Coords mcposaDISCTOF, mcposrDISCTOF;
Rotation mcrotaDISCTOF, mcrotrDISCTOF;
Coords mcposaPSDmon1Chopper, mcposrPSDmon1Chopper;
Rotation mcrotaPSDmon1Chopper, mcrotrPSDmon1Chopper;
Coords mcposaVacuumTube1entry, mcposrVacuumTube1entry;
Rotation mcrotaVacuumTube1entry, mcrotrVacuumTube1entry;
Coords mcposaVacuumTube1exit, mcposrVacuumTube1exit;
Rotation mcrotaVacuumTube1exit, mcrotrVacuumTube1exit;
Coords mcposaVacuumTube2entry, mcposrVacuumTube2entry;
Rotation mcrotaVacuumTube2entry, mcrotrVacuumTube2entry;
Coords mcposaVacuumTube2exit, mcposrVacuumTube2exit;
Rotation mcrotaVacuumTube2exit, mcrotrVacuumTube2exit;
Coords mcposaVacuumTube3entry, mcposrVacuumTube3entry;
Rotation mcrotaVacuumTube3entry, mcrotrVacuumTube3entry;
Coords mcposaVacuumTube3exit, mcposrVacuumTube3exit;
Rotation mcrotaVacuumTube3exit, mcrotrVacuumTube3exit;
Coords mcposaPSDmonMono, mcposrPSDmonMono;
Rotation mcrotaPSDmonMono, mcrotrPSDmonMono;
Coords mcposaMONOTOF, mcposrMONOTOF;
Rotation mcrotaMONOTOF, mcrotrMONOTOF;
Coords mcposaDivMonMono, mcposrDivMonMono;
Rotation mcrotaDivMonMono, mcrotrDivMonMono;
Coords mcposafocus_mono, mcposrfocus_mono;
Rotation mcrotafocus_mono, mcrotrfocus_mono;
Coords mcposamono, mcposrmono;
Rotation mcrotamono, mcrotrmono;
Coords mcposaa2, mcposra2;
Rotation mcrotaa2, mcrotra2;
Coords mcposaFERMITOF_before, mcposrFERMITOF_before;
Rotation mcrotaFERMITOF_before, mcrotrFERMITOF_before;
Coords mcposalambdaFermi, mcposrlambdaFermi;
Rotation mcrotalambdaFermi, mcrotrlambdaFermi;
Coords mcposaEMON_Fermi, mcposrEMON_Fermi;
Rotation mcrotaEMON_Fermi, mcrotrEMON_Fermi;
Coords mcposaDivMonfermi1, mcposrDivMonfermi1;
Rotation mcrotaDivMonfermi1, mcrotrDivMonfermi1;
Coords mcposaPSD_Fermi1, mcposrPSD_Fermi1;
Rotation mcrotaPSD_Fermi1, mcrotrPSD_Fermi1;
Coords mcposaFoChopper, mcposrFoChopper;
Rotation mcrotaFoChopper, mcrotrFoChopper;
Coords mcposaPSD_Fermi2, mcposrPSD_Fermi2;
Rotation mcrotaPSD_Fermi2, mcrotrPSD_Fermi2;
Coords mcposaDivMonfermi2, mcposrDivMonfermi2;
Rotation mcrotaDivMonfermi2, mcrotrDivMonfermi2;
Coords mcposaFERMITOF1, mcposrFERMITOF1;
Rotation mcrotaFERMITOF1, mcrotrFERMITOF1;
Coords mcposaSAMPLE_SLIT, mcposrSAMPLE_SLIT;
Rotation mcrotaSAMPLE_SLIT, mcrotrSAMPLE_SLIT;
Coords mcposaFERMITOF2, mcposrFERMITOF2;
Rotation mcrotaFERMITOF2, mcrotrFERMITOF2;
Coords mcposaPSD_SAMPLE, mcposrPSD_SAMPLE;
Rotation mcrotaPSD_SAMPLE, mcrotrPSD_SAMPLE;
Coords mcposaDivMon_Sample, mcposrDivMon_Sample;
Rotation mcrotaDivMon_Sample, mcrotrDivMon_Sample;
Coords mcposaEMON_SAMPLE, mcposrEMON_SAMPLE;
Rotation mcrotaEMON_SAMPLE, mcrotrEMON_SAMPLE;
Coords mcposaa3, mcposra3;
Rotation mcrotaa3, mcrotra3;
Coords mcposaSample, mcposrSample;
Rotation mcrotaSample, mcrotrSample;
Coords mcposaTOF_Det, mcposrTOF_Det;
Rotation mcrotaTOF_Det, mcrotrTOF_Det;
Coords mcposaFoDet, mcposrFoDet;
Rotation mcrotaFoDet, mcrotrFoDet;
Coords mcposaEMON_DET, mcposrEMON_DET;
Rotation mcrotaEMON_DET, mcrotrEMON_DET;

MCNUM mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz, mcnt, mcnsx, mcnsy, mcnsz, mcnp;

/* end declare */

void mcinit(void) {
#define mccompcurname  PSI_Focus
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaPSI_Focus coords_set(0,0,0)
#define lambda mciplambda
#define chopp_ratio mcipchopp_ratio
#define DET mcipDET
#line 62 "PSI_Focus.instr"
{

 /* setting theta and 2Theta */
 PHM = asin(lambda/2/monodd)*180/PI;
 TTM = 2*PHM;

 /* setting speed of Fermi chopper, Disk chopper */
 FERMI_SPEED=1.0e6/(PI*lambda*c*(1./tan(PHM*PI/180))*(dsd*pow((1),1.5)+dfcs)*(1.0-(dms/dgm)));
 DISC_SPEED = FERMI_SPEED/chopp_ratio;

 /* calculate the phase between Fermi and Disk chopper */
 FO_PHA = (c*lambda*(ddcm+mfc))/1e6;

 /* setting monochromator curvature */
 RV_2 = 2*dms*(sin(DEG2RAD*PHM));
 RH_2 = 2*dms/(sin(DEG2RAD*PHM));

 /* setting E and L monitors */
 LMIN = lambda - dL;
 LMAX = lambda + dL;
 EMIN = 81.81/(LMAX*LMAX);
 EMAX = 81.81/(LMIN*LMIN);

}
#line 13054 "PSI_Focus.c"
#undef DET
#undef chopp_ratio
#undef lambda
#undef mcposaPSI_Focus
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
#line 39 "PSI_Focus.instr"
  if("NULL") strncpy(mcca1_profile, "NULL" ? "NULL" : "", 16384); else mcca1_profile[0]='\0';
#line 39 "PSI_Focus.instr"
  mcca1_percent = 10;
#line 39 "PSI_Focus.instr"
  mcca1_flag_save = 0;
#line 39 "PSI_Focus.instr"
  mcca1_minutes = 0;
#line 13083 "PSI_Focus.c"

  SIG_MESSAGE("a1 (Init:Place/Rotate)");
  rot_set_rotation(mcrotaa1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 13090 "PSI_Focus.c"
  rot_copy(mcrotra1, mcrotaa1);
  mcposaa1 = coords_set(
#line 90 "PSI_Focus.instr"
    0,
#line 90 "PSI_Focus.instr"
    0,
#line 90 "PSI_Focus.instr"
    0);
#line 13099 "PSI_Focus.c"
  mctc1 = coords_neg(mcposaa1);
  mcposra1 = rot_apply(mcrotaa1, mctc1);
  mcDEBUG_COMPONENT("a1", mcposaa1, mcrotaa1)
  mccomp_posa[1] = mcposaa1;
  mccomp_posr[1] = mcposra1;
  mcNCounter[1]  = mcPCounter[1] = mcP2Counter[1] = 0;
  mcAbsorbProp[1]= 0;
    /* Component csource. */
  /* Setting parameters for component csource. */
  SIG_MESSAGE("csource (Init:SetPar)");
#line 129 "PSI_Focus.instr"
  if("NULL") strncpy(mcccsource_flux_file, "NULL" ? "NULL" : "", 16384); else mcccsource_flux_file[0]='\0';
#line 129 "PSI_Focus.instr"
  if("NULL") strncpy(mcccsource_xdiv_file, "NULL" ? "NULL" : "", 16384); else mcccsource_xdiv_file[0]='\0';
#line 129 "PSI_Focus.instr"
  if("NULL") strncpy(mcccsource_ydiv_file, "NULL" ? "NULL" : "", 16384); else mcccsource_ydiv_file[0]='\0';
#line 130 "PSI_Focus.instr"
  mcccsource_radius = 0.0;
#line 94 "PSI_Focus.instr"
  mcccsource_dist = 1.77;
#line 94 "PSI_Focus.instr"
  mcccsource_focus_xw = 0.05;
#line 94 "PSI_Focus.instr"
  mcccsource_focus_yh = 0.12;
#line 130 "PSI_Focus.instr"
  mcccsource_focus_aw = 0;
#line 130 "PSI_Focus.instr"
  mcccsource_focus_ah = 0;
#line 131 "PSI_Focus.instr"
  mcccsource_E0 = 0;
#line 131 "PSI_Focus.instr"
  mcccsource_dE = 0;
#line 95 "PSI_Focus.instr"
  mcccsource_lambda0 = 3.4;
#line 95 "PSI_Focus.instr"
  mcccsource_dlambda = dL;
#line 95 "PSI_Focus.instr"
  mcccsource_I1 = 8.5E11;
#line 94 "PSI_Focus.instr"
  mcccsource_yheight = 0.135;
#line 94 "PSI_Focus.instr"
  mcccsource_xwidth = 0.08;
#line 132 "PSI_Focus.instr"
  mcccsource_verbose = 0;
#line 95 "PSI_Focus.instr"
  mcccsource_T1 = 296.16;
#line 133 "PSI_Focus.instr"
  mcccsource_flux_file_perAA = 0;
#line 133 "PSI_Focus.instr"
  mcccsource_flux_file_log = 0;
#line 134 "PSI_Focus.instr"
  mcccsource_Lmin = 0;
#line 134 "PSI_Focus.instr"
  mcccsource_Lmax = 0;
#line 134 "PSI_Focus.instr"
  mcccsource_Emin = 0;
#line 134 "PSI_Focus.instr"
  mcccsource_Emax = 0;
#line 96 "PSI_Focus.instr"
  mcccsource_T2 = 40.68;
#line 96 "PSI_Focus.instr"
  mcccsource_I2 = 5.2E11;
#line 134 "PSI_Focus.instr"
  mcccsource_T3 = 0;
#line 134 "PSI_Focus.instr"
  mcccsource_I3 = 0;
#line 134 "PSI_Focus.instr"
  mcccsource_zdepth = 0;
#line 134 "PSI_Focus.instr"
  mcccsource_target_index = + 1;
#line 13170 "PSI_Focus.c"

  SIG_MESSAGE("csource (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 97 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 97 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 97 "PSI_Focus.instr"
    (0)*DEG2RAD);
#line 13180 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa1, mcrotacsource);
  rot_transpose(mcrotaa1, mctr1);
  rot_mul(mcrotacsource, mctr1, mcrotrcsource);
  mctc1 = coords_set(
#line 97 "PSI_Focus.instr"
    0,
#line 97 "PSI_Focus.instr"
    0,
#line 97 "PSI_Focus.instr"
    0);
#line 13191 "PSI_Focus.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposacsource = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaa1, mcposacsource);
  mcposrcsource = rot_apply(mcrotacsource, mctc1);
  mcDEBUG_COMPONENT("csource", mcposacsource, mcrotacsource)
  mccomp_posa[2] = mcposacsource;
  mccomp_posr[2] = mcposrcsource;
  mcNCounter[2]  = mcPCounter[2] = mcP2Counter[2] = 0;
  mcAbsorbProp[2]= 0;
    /* Component guide1. */
  /* Setting parameters for component guide1. */
  SIG_MESSAGE("guide1 (Init:SetPar)");
#line 58 "PSI_Focus.instr"
  if(0) strncpy(mccguide1_reflect, 0 ? 0 : "", 16384); else mccguide1_reflect[0]='\0';
#line 100 "PSI_Focus.instr"
  mccguide1_w1 = 0.05;
#line 100 "PSI_Focus.instr"
  mccguide1_h1 = 0.12;
#line 100 "PSI_Focus.instr"
  mccguide1_w2 = 0.05;
#line 100 "PSI_Focus.instr"
  mccguide1_h2 = 0.12;
#line 101 "PSI_Focus.instr"
  mccguide1_l = 4.66;
#line 101 "PSI_Focus.instr"
  mccguide1_R0 = 0.995;
#line 101 "PSI_Focus.instr"
  mccguide1_Qc = 0.0217;
#line 101 "PSI_Focus.instr"
  mccguide1_alpha = 5.76;
#line 102 "PSI_Focus.instr"
  mccguide1_m = 2.0;
#line 102 "PSI_Focus.instr"
  mccguide1_W = 0.0033;
#line 13227 "PSI_Focus.c"

  SIG_MESSAGE("guide1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 103 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 103 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 103 "PSI_Focus.instr"
    (0)*DEG2RAD);
#line 13237 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa1, mcrotaguide1);
  rot_transpose(mcrotacsource, mctr1);
  rot_mul(mcrotaguide1, mctr1, mcrotrguide1);
  mctc1 = coords_set(
#line 103 "PSI_Focus.instr"
    0,
#line 103 "PSI_Focus.instr"
    0,
#line 103 "PSI_Focus.instr"
    1.5);
#line 13248 "PSI_Focus.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaguide1 = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposacsource, mcposaguide1);
  mcposrguide1 = rot_apply(mcrotaguide1, mctc1);
  mcDEBUG_COMPONENT("guide1", mcposaguide1, mcrotaguide1)
  mccomp_posa[3] = mcposaguide1;
  mccomp_posr[3] = mcposrguide1;
  mcNCounter[3]  = mcPCounter[3] = mcP2Counter[3] = 0;
  mcAbsorbProp[3]= 0;
    /* Component guide2. */
  /* Setting parameters for component guide2. */
  SIG_MESSAGE("guide2 (Init:SetPar)");
#line 107 "PSI_Focus.instr"
  mccguide2_w = 0.05;
#line 107 "PSI_Focus.instr"
  mccguide2_h = 0.12;
#line 107 "PSI_Focus.instr"
  mccguide2_r = 1445;
#line 94 "PSI_Focus.instr"
  mccguide2_Win = 0.04;
#line 94 "PSI_Focus.instr"
  mccguide2_k = 1;
#line 94 "PSI_Focus.instr"
  mccguide2_d = 0.001;
#line 110 "PSI_Focus.instr"
  mccguide2_l = 24.5;
#line 107 "PSI_Focus.instr"
  mccguide2_R0a = 0.995;
#line 107 "PSI_Focus.instr"
  mccguide2_Qca = 0.0217;
#line 108 "PSI_Focus.instr"
  mccguide2_alphaa = 5.76;
#line 108 "PSI_Focus.instr"
  mccguide2_ma = 2;
#line 108 "PSI_Focus.instr"
  mccguide2_Wa = 0.0033;
#line 108 "PSI_Focus.instr"
  mccguide2_R0i = 0.995;
#line 108 "PSI_Focus.instr"
  mccguide2_Qci = 0.0217;
#line 109 "PSI_Focus.instr"
  mccguide2_alphai = 5.76;
#line 109 "PSI_Focus.instr"
  mccguide2_mi = 2;
#line 109 "PSI_Focus.instr"
  mccguide2_Wi = 0.0033;
#line 109 "PSI_Focus.instr"
  mccguide2_R0s = 0.995;
#line 109 "PSI_Focus.instr"
  mccguide2_Qcs = 0.0217;
#line 110 "PSI_Focus.instr"
  mccguide2_alphas = 5.76;
#line 110 "PSI_Focus.instr"
  mccguide2_ms = 2;
#line 110 "PSI_Focus.instr"
  mccguide2_Ws = 0.0033;
#line 13306 "PSI_Focus.c"

  SIG_MESSAGE("guide2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 111 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 111 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 111 "PSI_Focus.instr"
    (0)*DEG2RAD);
#line 13316 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa1, mcrotaguide2);
  rot_transpose(mcrotaguide1, mctr1);
  rot_mul(mcrotaguide2, mctr1, mcrotrguide2);
  mctc1 = coords_set(
#line 111 "PSI_Focus.instr"
    0,
#line 111 "PSI_Focus.instr"
    0,
#line 111 "PSI_Focus.instr"
    6.16);
#line 13327 "PSI_Focus.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaguide2 = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaguide1, mcposaguide2);
  mcposrguide2 = rot_apply(mcrotaguide2, mctc1);
  mcDEBUG_COMPONENT("guide2", mcposaguide2, mcrotaguide2)
  mccomp_posa[4] = mcposaguide2;
  mccomp_posr[4] = mcposrguide2;
  mcNCounter[4]  = mcPCounter[4] = mcP2Counter[4] = 0;
  mcAbsorbProp[4]= 0;
    /* Component bunker. */
  /* Setting parameters for component bunker. */
  SIG_MESSAGE("bunker (Init:SetPar)");
#line 58 "PSI_Focus.instr"
  if(0) strncpy(mccbunker_reflect, 0 ? 0 : "", 16384); else mccbunker_reflect[0]='\0';
#line 114 "PSI_Focus.instr"
  mccbunker_w1 = 0.05;
#line 114 "PSI_Focus.instr"
  mccbunker_h1 = .12;
#line 114 "PSI_Focus.instr"
  mccbunker_w2 = 0.05;
#line 114 "PSI_Focus.instr"
  mccbunker_h2 = .12;
#line 115 "PSI_Focus.instr"
  mccbunker_l = 3.0;
#line 115 "PSI_Focus.instr"
  mccbunker_R0 = 0.995;
#line 115 "PSI_Focus.instr"
  mccbunker_Qc = 0.0217;
#line 115 "PSI_Focus.instr"
  mccbunker_alpha = 5.76;
#line 115 "PSI_Focus.instr"
  mccbunker_m = 2.0;
#line 115 "PSI_Focus.instr"
  mccbunker_W = 0.0033;
#line 13363 "PSI_Focus.c"

  SIG_MESSAGE("bunker (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 116 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 116 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 116 "PSI_Focus.instr"
    (0)*DEG2RAD);
#line 13373 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa1, mcrotabunker);
  rot_transpose(mcrotaguide2, mctr1);
  rot_mul(mcrotabunker, mctr1, mcrotrbunker);
  mctc1 = coords_set(
#line 116 "PSI_Focus.instr"
    0,
#line 116 "PSI_Focus.instr"
    0,
#line 116 "PSI_Focus.instr"
    30.66);
#line 13384 "PSI_Focus.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposabunker = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaguide2, mcposabunker);
  mcposrbunker = rot_apply(mcrotabunker, mctc1);
  mcDEBUG_COMPONENT("bunker", mcposabunker, mcrotabunker)
  mccomp_posa[5] = mcposabunker;
  mccomp_posr[5] = mcposrbunker;
  mcNCounter[5]  = mcPCounter[5] = mcP2Counter[5] = 0;
  mcAbsorbProp[5]= 0;
    /* Component guide3. */
  /* Setting parameters for component guide3. */
  SIG_MESSAGE("guide3 (Init:SetPar)");
#line 58 "PSI_Focus.instr"
  if(0) strncpy(mccguide3_reflect, 0 ? 0 : "", 16384); else mccguide3_reflect[0]='\0';
#line 119 "PSI_Focus.instr"
  mccguide3_w1 = 0.05;
#line 119 "PSI_Focus.instr"
  mccguide3_h1 = .12;
#line 119 "PSI_Focus.instr"
  mccguide3_w2 = 0.05;
#line 119 "PSI_Focus.instr"
  mccguide3_h2 = .12;
#line 120 "PSI_Focus.instr"
  mccguide3_l = 32.95;
#line 120 "PSI_Focus.instr"
  mccguide3_R0 = 0.995;
#line 120 "PSI_Focus.instr"
  mccguide3_Qc = 0.0217;
#line 120 "PSI_Focus.instr"
  mccguide3_alpha = 5.76;
#line 120 "PSI_Focus.instr"
  mccguide3_m = 2.0;
#line 120 "PSI_Focus.instr"
  mccguide3_W = 0.0033;
#line 13420 "PSI_Focus.c"

  SIG_MESSAGE("guide3 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 121 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 121 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 121 "PSI_Focus.instr"
    (0)*DEG2RAD);
#line 13430 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa1, mcrotaguide3);
  rot_transpose(mcrotabunker, mctr1);
  rot_mul(mcrotaguide3, mctr1, mcrotrguide3);
  mctc1 = coords_set(
#line 121 "PSI_Focus.instr"
    0,
#line 121 "PSI_Focus.instr"
    0,
#line 121 "PSI_Focus.instr"
    33.66);
#line 13441 "PSI_Focus.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaguide3 = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposabunker, mcposaguide3);
  mcposrguide3 = rot_apply(mcrotaguide3, mctc1);
  mcDEBUG_COMPONENT("guide3", mcposaguide3, mcrotaguide3)
  mccomp_posa[6] = mcposaguide3;
  mccomp_posr[6] = mcposrguide3;
  mcNCounter[6]  = mcPCounter[6] = mcP2Counter[6] = 0;
  mcAbsorbProp[6]= 0;
    /* Component lambdaGuideExit. */
  /* Setting parameters for component lambdaGuideExit. */
  SIG_MESSAGE("lambdaGuideExit (Init:SetPar)");
#line 128 "PSI_Focus.instr"
  if("lambdaguide.dat") strncpy(mcclambdaGuideExit_filename, "lambdaguide.dat" ? "lambdaguide.dat" : "", 16384); else mcclambdaGuideExit_filename[0]='\0';
#line 50 "PSI_Focus.instr"
  mcclambdaGuideExit_xmin = -0.05;
#line 50 "PSI_Focus.instr"
  mcclambdaGuideExit_xmax = 0.05;
#line 50 "PSI_Focus.instr"
  mcclambdaGuideExit_ymin = -0.05;
#line 50 "PSI_Focus.instr"
  mcclambdaGuideExit_ymax = 0.05;
#line 126 "PSI_Focus.instr"
  mcclambdaGuideExit_xwidth = 0.05;
#line 126 "PSI_Focus.instr"
  mcclambdaGuideExit_yheight = 0.12;
#line 127 "PSI_Focus.instr"
  mcclambdaGuideExit_Lmin = LMIN;
#line 127 "PSI_Focus.instr"
  mcclambdaGuideExit_Lmax = LMAX;
#line 51 "PSI_Focus.instr"
  mcclambdaGuideExit_restore_neutron = 0;
#line 51 "PSI_Focus.instr"
  mcclambdaGuideExit_nowritefile = 0;
#line 13477 "PSI_Focus.c"

  SIG_MESSAGE("lambdaGuideExit (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 129 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 129 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 129 "PSI_Focus.instr"
    (0)*DEG2RAD);
#line 13487 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa1, mcrotalambdaGuideExit);
  rot_transpose(mcrotaguide3, mctr1);
  rot_mul(mcrotalambdaGuideExit, mctr1, mcrotrlambdaGuideExit);
  mctc1 = coords_set(
#line 129 "PSI_Focus.instr"
    0,
#line 129 "PSI_Focus.instr"
    0,
#line 129 "PSI_Focus.instr"
    66.611);
#line 13498 "PSI_Focus.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposalambdaGuideExit = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaguide3, mcposalambdaGuideExit);
  mcposrlambdaGuideExit = rot_apply(mcrotalambdaGuideExit, mctc1);
  mcDEBUG_COMPONENT("lambdaGuideExit", mcposalambdaGuideExit, mcrotalambdaGuideExit)
  mccomp_posa[7] = mcposalambdaGuideExit;
  mccomp_posr[7] = mcposrlambdaGuideExit;
  mcNCounter[7]  = mcPCounter[7] = mcP2Counter[7] = 0;
  mcAbsorbProp[7]= 0;
    /* Component DivMonGuideExit. */
  /* Setting parameters for component DivMonGuideExit. */
  SIG_MESSAGE("DivMonGuideExit (Init:SetPar)");
#line 132 "PSI_Focus.instr"
  if("divguide.dat") strncpy(mccDivMonGuideExit_filename, "divguide.dat" ? "divguide.dat" : "", 16384); else mccDivMonGuideExit_filename[0]='\0';
#line 55 "PSI_Focus.instr"
  mccDivMonGuideExit_xmin = -0.05;
#line 55 "PSI_Focus.instr"
  mccDivMonGuideExit_xmax = 0.05;
#line 55 "PSI_Focus.instr"
  mccDivMonGuideExit_ymin = -0.05;
#line 55 "PSI_Focus.instr"
  mccDivMonGuideExit_ymax = 0.05;
#line 133 "PSI_Focus.instr"
  mccDivMonGuideExit_xwidth = 0.05;
#line 133 "PSI_Focus.instr"
  mccDivMonGuideExit_yheight = 0.12;
#line 134 "PSI_Focus.instr"
  mccDivMonGuideExit_maxdiv_h = 2;
#line 134 "PSI_Focus.instr"
  mccDivMonGuideExit_maxdiv_v = 2;
#line 56 "PSI_Focus.instr"
  mccDivMonGuideExit_restore_neutron = 0;
#line 56 "PSI_Focus.instr"
  mccDivMonGuideExit_nx = 0;
#line 56 "PSI_Focus.instr"
  mccDivMonGuideExit_ny = 0;
#line 56 "PSI_Focus.instr"
  mccDivMonGuideExit_nz = 1;
#line 56 "PSI_Focus.instr"
  mccDivMonGuideExit_nowritefile = 0;
#line 13540 "PSI_Focus.c"

  SIG_MESSAGE("DivMonGuideExit (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 135 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 135 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 135 "PSI_Focus.instr"
    (0)*DEG2RAD);
#line 13550 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa1, mcrotaDivMonGuideExit);
  rot_transpose(mcrotalambdaGuideExit, mctr1);
  rot_mul(mcrotaDivMonGuideExit, mctr1, mcrotrDivMonGuideExit);
  mctc1 = coords_set(
#line 135 "PSI_Focus.instr"
    0,
#line 135 "PSI_Focus.instr"
    0,
#line 135 "PSI_Focus.instr"
    66.612);
#line 13561 "PSI_Focus.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaDivMonGuideExit = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposalambdaGuideExit, mcposaDivMonGuideExit);
  mcposrDivMonGuideExit = rot_apply(mcrotaDivMonGuideExit, mctc1);
  mcDEBUG_COMPONENT("DivMonGuideExit", mcposaDivMonGuideExit, mcrotaDivMonGuideExit)
  mccomp_posa[8] = mcposaDivMonGuideExit;
  mccomp_posr[8] = mcposrDivMonGuideExit;
  mcNCounter[8]  = mcPCounter[8] = mcP2Counter[8] = 0;
  mcAbsorbProp[8]= 0;
    /* Component PSDGuideExit. */
  /* Setting parameters for component PSDGuideExit. */
  SIG_MESSAGE("PSDGuideExit (Init:SetPar)");
#line 139 "PSI_Focus.instr"
  if("psdguide.dat") strncpy(mccPSDGuideExit_filename, "psdguide.dat" ? "psdguide.dat" : "", 16384); else mccPSDGuideExit_filename[0]='\0';
#line 50 "PSI_Focus.instr"
  mccPSDGuideExit_xmin = -0.05;
#line 50 "PSI_Focus.instr"
  mccPSDGuideExit_xmax = 0.05;
#line 50 "PSI_Focus.instr"
  mccPSDGuideExit_ymin = -0.05;
#line 50 "PSI_Focus.instr"
  mccPSDGuideExit_ymax = 0.05;
#line 138 "PSI_Focus.instr"
  mccPSDGuideExit_xwidth = 0.05;
#line 138 "PSI_Focus.instr"
  mccPSDGuideExit_yheight = 0.12;
#line 50 "PSI_Focus.instr"
  mccPSDGuideExit_restore_neutron = 0;
#line 50 "PSI_Focus.instr"
  mccPSDGuideExit_nowritefile = 0;
#line 13593 "PSI_Focus.c"

  SIG_MESSAGE("PSDGuideExit (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 140 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 140 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 140 "PSI_Focus.instr"
    (0)*DEG2RAD);
#line 13603 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa1, mcrotaPSDGuideExit);
  rot_transpose(mcrotaDivMonGuideExit, mctr1);
  rot_mul(mcrotaPSDGuideExit, mctr1, mcrotrPSDGuideExit);
  mctc1 = coords_set(
#line 140 "PSI_Focus.instr"
    0,
#line 140 "PSI_Focus.instr"
    0,
#line 140 "PSI_Focus.instr"
    66.613);
#line 13614 "PSI_Focus.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPSDGuideExit = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaDivMonGuideExit, mcposaPSDGuideExit);
  mcposrPSDGuideExit = rot_apply(mcrotaPSDGuideExit, mctc1);
  mcDEBUG_COMPONENT("PSDGuideExit", mcposaPSDGuideExit, mcrotaPSDGuideExit)
  mccomp_posa[9] = mcposaPSDGuideExit;
  mccomp_posr[9] = mcposrPSDGuideExit;
  mcNCounter[9]  = mcPCounter[9] = mcP2Counter[9] = 0;
  mcAbsorbProp[9]= 0;
    /* Component FOCUSguide. */
  /* Setting parameters for component FOCUSguide. */
  SIG_MESSAGE("FOCUSguide (Init:SetPar)");
#line 145 "PSI_Focus.instr"
  mccFOCUSguide_w1 = 0.05;
#line 145 "PSI_Focus.instr"
  mccFOCUSguide_h1 = 0.12;
#line 146 "PSI_Focus.instr"
  mccFOCUSguide_w2 = 0.05;
#line 146 "PSI_Focus.instr"
  mccFOCUSguide_h2 = 0.095;
#line 147 "PSI_Focus.instr"
  mccFOCUSguide_l = 2.95;
#line 147 "PSI_Focus.instr"
  mccFOCUSguide_R0 = 0.995;
#line 71 "PSI_Focus.instr"
  mccFOCUSguide_Qc = 0;
#line 71 "PSI_Focus.instr"
  mccFOCUSguide_alpha = 0;
#line 71 "PSI_Focus.instr"
  mccFOCUSguide_m = 0;
#line 149 "PSI_Focus.instr"
  mccFOCUSguide_nslit = 1;
#line 149 "PSI_Focus.instr"
  mccFOCUSguide_d = 0.002;
#line 147 "PSI_Focus.instr"
  mccFOCUSguide_Qcx = 0.0217;
#line 147 "PSI_Focus.instr"
  mccFOCUSguide_Qcy = 0.0217;
#line 148 "PSI_Focus.instr"
  mccFOCUSguide_alphax = 5.76;
#line 148 "PSI_Focus.instr"
  mccFOCUSguide_alphay = 5.64;
#line 148 "PSI_Focus.instr"
  mccFOCUSguide_W = 0.0033;
#line 149 "PSI_Focus.instr"
  mccFOCUSguide_mx = 2.4;
#line 149 "PSI_Focus.instr"
  mccFOCUSguide_my = 3;
#line 72 "PSI_Focus.instr"
  mccFOCUSguide_nu = 0;
#line 72 "PSI_Focus.instr"
  mccFOCUSguide_phase = 0;
#line 13668 "PSI_Focus.c"

  SIG_MESSAGE("FOCUSguide (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 150 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 150 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 150 "PSI_Focus.instr"
    (0)*DEG2RAD);
#line 13678 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa1, mcrotaFOCUSguide);
  rot_transpose(mcrotaPSDGuideExit, mctr1);
  rot_mul(mcrotaFOCUSguide, mctr1, mcrotrFOCUSguide);
  mctc1 = coords_set(
#line 150 "PSI_Focus.instr"
    0,
#line 150 "PSI_Focus.instr"
    0,
#line 150 "PSI_Focus.instr"
    66.67);
#line 13689 "PSI_Focus.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaFOCUSguide = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaPSDGuideExit, mcposaFOCUSguide);
  mcposrFOCUSguide = rot_apply(mcrotaFOCUSguide, mctc1);
  mcDEBUG_COMPONENT("FOCUSguide", mcposaFOCUSguide, mcrotaFOCUSguide)
  mccomp_posa[10] = mcposaFOCUSguide;
  mccomp_posr[10] = mcposrFOCUSguide;
  mcNCounter[10]  = mcPCounter[10] = mcP2Counter[10] = 0;
  mcAbsorbProp[10]= 0;
    /* Component FirstChopper. */
  /* Setting parameters for component FirstChopper. */
  SIG_MESSAGE("FirstChopper (Init:SetPar)");
#line 154 "PSI_Focus.instr"
  mccFirstChopper_theta_0 = 0;
#line 155 "PSI_Focus.instr"
  mccFirstChopper_radius = 0.70;
#line 155 "PSI_Focus.instr"
  mccFirstChopper_yheight = 0.70;
#line 155 "PSI_Focus.instr"
  mccFirstChopper_nu = DISC_SPEED / 2 / PI;
#line 155 "PSI_Focus.instr"
  mccFirstChopper_nslit = 2;
#line 59 "PSI_Focus.instr"
  mccFirstChopper_jitter = 0;
#line 155 "PSI_Focus.instr"
  mccFirstChopper_delay = FO_PHA;
#line 155 "PSI_Focus.instr"
  mccFirstChopper_isfirst = 1;
#line 59 "PSI_Focus.instr"
  mccFirstChopper_n_pulse = 1;
#line 154 "PSI_Focus.instr"
  mccFirstChopper_abs_out = 0;
#line 59 "PSI_Focus.instr"
  mccFirstChopper_phase = 0;
#line 155 "PSI_Focus.instr"
  mccFirstChopper_xwidth = 0.07;
#line 154 "PSI_Focus.instr"
  mccFirstChopper_verbose = 1;
#line 13729 "PSI_Focus.c"

  SIG_MESSAGE("FirstChopper (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 13736 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa1, mcrotaFirstChopper);
  rot_transpose(mcrotaFOCUSguide, mctr1);
  rot_mul(mcrotaFirstChopper, mctr1, mcrotrFirstChopper);
  mctc1 = coords_set(
#line 156 "PSI_Focus.instr"
    0,
#line 156 "PSI_Focus.instr"
    0,
#line 156 "PSI_Focus.instr"
    69.674);
#line 13747 "PSI_Focus.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaFirstChopper = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaFOCUSguide, mcposaFirstChopper);
  mcposrFirstChopper = rot_apply(mcrotaFirstChopper, mctc1);
  mcDEBUG_COMPONENT("FirstChopper", mcposaFirstChopper, mcrotaFirstChopper)
  mccomp_posa[11] = mcposaFirstChopper;
  mccomp_posr[11] = mcposrFirstChopper;
  mcNCounter[11]  = mcPCounter[11] = mcP2Counter[11] = 0;
  mcAbsorbProp[11]= 0;
    /* Component DISCTOF. */
  /* Setting parameters for component DISCTOF. */
  SIG_MESSAGE("DISCTOF (Init:SetPar)");
#line 201 "PSI_Focus.instr"
  mccDISCTOF_xwidth = 0;
#line 201 "PSI_Focus.instr"
  mccDISCTOF_yheight = 0;
#line 201 "PSI_Focus.instr"
  mccDISCTOF_zdepth = 0;
#line 159 "PSI_Focus.instr"
  mccDISCTOF_xmin = -0.05;
#line 160 "PSI_Focus.instr"
  mccDISCTOF_xmax = 0.05;
#line 160 "PSI_Focus.instr"
  mccDISCTOF_ymin = -0.1;
#line 160 "PSI_Focus.instr"
  mccDISCTOF_ymax = 0.1;
#line 202 "PSI_Focus.instr"
  mccDISCTOF_zmin = 0;
#line 202 "PSI_Focus.instr"
  mccDISCTOF_zmax = 0;
#line 159 "PSI_Focus.instr"
  mccDISCTOF_bins = 30;
#line 203 "PSI_Focus.instr"
  mccDISCTOF_min = -1e40;
#line 203 "PSI_Focus.instr"
  mccDISCTOF_max = 1e40;
#line 203 "PSI_Focus.instr"
  mccDISCTOF_restore_neutron = 0;
#line 203 "PSI_Focus.instr"
  mccDISCTOF_radius = 0;
#line 160 "PSI_Focus.instr"
  if("auto, time") strncpy(mccDISCTOF_options, "auto, time" ? "auto, time" : "", 16384); else mccDISCTOF_options[0]='\0';
#line 159 "PSI_Focus.instr"
  if("DISC_TOF.dat") strncpy(mccDISCTOF_filename, "DISC_TOF.dat" ? "DISC_TOF.dat" : "", 16384); else mccDISCTOF_filename[0]='\0';
#line 204 "PSI_Focus.instr"
  if("NULL") strncpy(mccDISCTOF_geometry, "NULL" ? "NULL" : "", 16384); else mccDISCTOF_geometry[0]='\0';
#line 205 "PSI_Focus.instr"
  if("NULL") strncpy(mccDISCTOF_username1, "NULL" ? "NULL" : "", 16384); else mccDISCTOF_username1[0]='\0';
#line 205 "PSI_Focus.instr"
  if("NULL") strncpy(mccDISCTOF_username2, "NULL" ? "NULL" : "", 16384); else mccDISCTOF_username2[0]='\0';
#line 205 "PSI_Focus.instr"
  if("NULL") strncpy(mccDISCTOF_username3, "NULL" ? "NULL" : "", 16384); else mccDISCTOF_username3[0]='\0';
#line 206 "PSI_Focus.instr"
  mccDISCTOF_nowritefile = 0;
#line 13803 "PSI_Focus.c"

  SIG_MESSAGE("DISCTOF (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 13810 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa1, mcrotaDISCTOF);
  rot_transpose(mcrotaFirstChopper, mctr1);
  rot_mul(mcrotaDISCTOF, mctr1, mcrotrDISCTOF);
  mctc1 = coords_set(
#line 161 "PSI_Focus.instr"
    0,
#line 161 "PSI_Focus.instr"
    0,
#line 161 "PSI_Focus.instr"
    69.718);
#line 13821 "PSI_Focus.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaDISCTOF = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaFirstChopper, mcposaDISCTOF);
  mcposrDISCTOF = rot_apply(mcrotaDISCTOF, mctc1);
  mcDEBUG_COMPONENT("DISCTOF", mcposaDISCTOF, mcrotaDISCTOF)
  mccomp_posa[12] = mcposaDISCTOF;
  mccomp_posr[12] = mcposrDISCTOF;
  mcNCounter[12]  = mcPCounter[12] = mcP2Counter[12] = 0;
  mcAbsorbProp[12]= 0;
    /* Component PSDmon1Chopper. */
  /* Setting parameters for component PSDmon1Chopper. */
  SIG_MESSAGE("PSDmon1Chopper (Init:SetPar)");
#line 165 "PSI_Focus.instr"
  if("psdchopper.dat") strncpy(mccPSDmon1Chopper_filename, "psdchopper.dat" ? "psdchopper.dat" : "", 16384); else mccPSDmon1Chopper_filename[0]='\0';
#line 50 "PSI_Focus.instr"
  mccPSDmon1Chopper_xmin = -0.05;
#line 50 "PSI_Focus.instr"
  mccPSDmon1Chopper_xmax = 0.05;
#line 50 "PSI_Focus.instr"
  mccPSDmon1Chopper_ymin = -0.05;
#line 50 "PSI_Focus.instr"
  mccPSDmon1Chopper_ymax = 0.05;
#line 164 "PSI_Focus.instr"
  mccPSDmon1Chopper_xwidth = 0.1;
#line 164 "PSI_Focus.instr"
  mccPSDmon1Chopper_yheight = 0.2;
#line 50 "PSI_Focus.instr"
  mccPSDmon1Chopper_restore_neutron = 0;
#line 50 "PSI_Focus.instr"
  mccPSDmon1Chopper_nowritefile = 0;
#line 13853 "PSI_Focus.c"

  SIG_MESSAGE("PSDmon1Chopper (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 166 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 166 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 166 "PSI_Focus.instr"
    (0)*DEG2RAD);
#line 13863 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa1, mcrotaPSDmon1Chopper);
  rot_transpose(mcrotaDISCTOF, mctr1);
  rot_mul(mcrotaPSDmon1Chopper, mctr1, mcrotrPSDmon1Chopper);
  mctc1 = coords_set(
#line 166 "PSI_Focus.instr"
    0,
#line 166 "PSI_Focus.instr"
    0,
#line 166 "PSI_Focus.instr"
    69.719);
#line 13874 "PSI_Focus.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPSDmon1Chopper = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaDISCTOF, mcposaPSDmon1Chopper);
  mcposrPSDmon1Chopper = rot_apply(mcrotaPSDmon1Chopper, mctc1);
  mcDEBUG_COMPONENT("PSDmon1Chopper", mcposaPSDmon1Chopper, mcrotaPSDmon1Chopper)
  mccomp_posa[13] = mcposaPSDmon1Chopper;
  mccomp_posr[13] = mcposrPSDmon1Chopper;
  mcNCounter[13]  = mcPCounter[13] = mcP2Counter[13] = 0;
  mcAbsorbProp[13]= 0;
    /* Component VacuumTube1entry. */
  /* Setting parameters for component VacuumTube1entry. */
  SIG_MESSAGE("VacuumTube1entry (Init:SetPar)");
#line 46 "PSI_Focus.instr"
  mccVacuumTube1entry_xmin = 0;
#line 46 "PSI_Focus.instr"
  mccVacuumTube1entry_xmax = 0;
#line 46 "PSI_Focus.instr"
  mccVacuumTube1entry_ymin = 0;
#line 46 "PSI_Focus.instr"
  mccVacuumTube1entry_ymax = 0;
#line 182 "PSI_Focus.instr"
  mccVacuumTube1entry_radius = 0.085;
#line 46 "PSI_Focus.instr"
  mccVacuumTube1entry_xwidth = 0;
#line 46 "PSI_Focus.instr"
  mccVacuumTube1entry_yheight = 0;
#line 13902 "PSI_Focus.c"

  SIG_MESSAGE("VacuumTube1entry (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 13909 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa1, mcrotaVacuumTube1entry);
  rot_transpose(mcrotaPSDmon1Chopper, mctr1);
  rot_mul(mcrotaVacuumTube1entry, mctr1, mcrotrVacuumTube1entry);
  mctc1 = coords_set(
#line 183 "PSI_Focus.instr"
    0,
#line 183 "PSI_Focus.instr"
    0,
#line 183 "PSI_Focus.instr"
    70.072);
#line 13920 "PSI_Focus.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaVacuumTube1entry = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaPSDmon1Chopper, mcposaVacuumTube1entry);
  mcposrVacuumTube1entry = rot_apply(mcrotaVacuumTube1entry, mctc1);
  mcDEBUG_COMPONENT("VacuumTube1entry", mcposaVacuumTube1entry, mcrotaVacuumTube1entry)
  mccomp_posa[14] = mcposaVacuumTube1entry;
  mccomp_posr[14] = mcposrVacuumTube1entry;
  mcNCounter[14]  = mcPCounter[14] = mcP2Counter[14] = 0;
  mcAbsorbProp[14]= 0;
    /* Component VacuumTube1exit. */
  /* Setting parameters for component VacuumTube1exit. */
  SIG_MESSAGE("VacuumTube1exit (Init:SetPar)");
#line 46 "PSI_Focus.instr"
  mccVacuumTube1exit_xmin = 0;
#line 46 "PSI_Focus.instr"
  mccVacuumTube1exit_xmax = 0;
#line 46 "PSI_Focus.instr"
  mccVacuumTube1exit_ymin = 0;
#line 46 "PSI_Focus.instr"
  mccVacuumTube1exit_ymax = 0;
#line 185 "PSI_Focus.instr"
  mccVacuumTube1exit_radius = 0.085;
#line 46 "PSI_Focus.instr"
  mccVacuumTube1exit_xwidth = 0;
#line 46 "PSI_Focus.instr"
  mccVacuumTube1exit_yheight = 0;
#line 13948 "PSI_Focus.c"

  SIG_MESSAGE("VacuumTube1exit (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 13955 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa1, mcrotaVacuumTube1exit);
  rot_transpose(mcrotaVacuumTube1entry, mctr1);
  rot_mul(mcrotaVacuumTube1exit, mctr1, mcrotrVacuumTube1exit);
  mctc1 = coords_set(
#line 186 "PSI_Focus.instr"
    0,
#line 186 "PSI_Focus.instr"
    0,
#line 186 "PSI_Focus.instr"
    70.607);
#line 13966 "PSI_Focus.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaVacuumTube1exit = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaVacuumTube1entry, mcposaVacuumTube1exit);
  mcposrVacuumTube1exit = rot_apply(mcrotaVacuumTube1exit, mctc1);
  mcDEBUG_COMPONENT("VacuumTube1exit", mcposaVacuumTube1exit, mcrotaVacuumTube1exit)
  mccomp_posa[15] = mcposaVacuumTube1exit;
  mccomp_posr[15] = mcposrVacuumTube1exit;
  mcNCounter[15]  = mcPCounter[15] = mcP2Counter[15] = 0;
  mcAbsorbProp[15]= 0;
    /* Component VacuumTube2entry. */
  /* Setting parameters for component VacuumTube2entry. */
  SIG_MESSAGE("VacuumTube2entry (Init:SetPar)");
#line 46 "PSI_Focus.instr"
  mccVacuumTube2entry_xmin = 0;
#line 46 "PSI_Focus.instr"
  mccVacuumTube2entry_xmax = 0;
#line 46 "PSI_Focus.instr"
  mccVacuumTube2entry_ymin = 0;
#line 46 "PSI_Focus.instr"
  mccVacuumTube2entry_ymax = 0;
#line 190 "PSI_Focus.instr"
  mccVacuumTube2entry_radius = 0.1;
#line 46 "PSI_Focus.instr"
  mccVacuumTube2entry_xwidth = 0;
#line 46 "PSI_Focus.instr"
  mccVacuumTube2entry_yheight = 0;
#line 13994 "PSI_Focus.c"

  SIG_MESSAGE("VacuumTube2entry (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 14001 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa1, mcrotaVacuumTube2entry);
  rot_transpose(mcrotaVacuumTube1exit, mctr1);
  rot_mul(mcrotaVacuumTube2entry, mctr1, mcrotrVacuumTube2entry);
  mctc1 = coords_set(
#line 191 "PSI_Focus.instr"
    0,
#line 191 "PSI_Focus.instr"
    0,
#line 191 "PSI_Focus.instr"
    70.607);
#line 14012 "PSI_Focus.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaVacuumTube2entry = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaVacuumTube1exit, mcposaVacuumTube2entry);
  mcposrVacuumTube2entry = rot_apply(mcrotaVacuumTube2entry, mctc1);
  mcDEBUG_COMPONENT("VacuumTube2entry", mcposaVacuumTube2entry, mcrotaVacuumTube2entry)
  mccomp_posa[16] = mcposaVacuumTube2entry;
  mccomp_posr[16] = mcposrVacuumTube2entry;
  mcNCounter[16]  = mcPCounter[16] = mcP2Counter[16] = 0;
  mcAbsorbProp[16]= 0;
    /* Component VacuumTube2exit. */
  /* Setting parameters for component VacuumTube2exit. */
  SIG_MESSAGE("VacuumTube2exit (Init:SetPar)");
#line 46 "PSI_Focus.instr"
  mccVacuumTube2exit_xmin = 0;
#line 46 "PSI_Focus.instr"
  mccVacuumTube2exit_xmax = 0;
#line 46 "PSI_Focus.instr"
  mccVacuumTube2exit_ymin = 0;
#line 46 "PSI_Focus.instr"
  mccVacuumTube2exit_ymax = 0;
#line 193 "PSI_Focus.instr"
  mccVacuumTube2exit_radius = 0.1;
#line 46 "PSI_Focus.instr"
  mccVacuumTube2exit_xwidth = 0;
#line 46 "PSI_Focus.instr"
  mccVacuumTube2exit_yheight = 0;
#line 14040 "PSI_Focus.c"

  SIG_MESSAGE("VacuumTube2exit (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 14047 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa1, mcrotaVacuumTube2exit);
  rot_transpose(mcrotaVacuumTube2entry, mctr1);
  rot_mul(mcrotaVacuumTube2exit, mctr1, mcrotrVacuumTube2exit);
  mctc1 = coords_set(
#line 194 "PSI_Focus.instr"
    0,
#line 194 "PSI_Focus.instr"
    0,
#line 194 "PSI_Focus.instr"
    71.402);
#line 14058 "PSI_Focus.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaVacuumTube2exit = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaVacuumTube2entry, mcposaVacuumTube2exit);
  mcposrVacuumTube2exit = rot_apply(mcrotaVacuumTube2exit, mctc1);
  mcDEBUG_COMPONENT("VacuumTube2exit", mcposaVacuumTube2exit, mcrotaVacuumTube2exit)
  mccomp_posa[17] = mcposaVacuumTube2exit;
  mccomp_posr[17] = mcposrVacuumTube2exit;
  mcNCounter[17]  = mcPCounter[17] = mcP2Counter[17] = 0;
  mcAbsorbProp[17]= 0;
    /* Component VacuumTube3entry. */
  /* Setting parameters for component VacuumTube3entry. */
  SIG_MESSAGE("VacuumTube3entry (Init:SetPar)");
#line 46 "PSI_Focus.instr"
  mccVacuumTube3entry_xmin = 0;
#line 46 "PSI_Focus.instr"
  mccVacuumTube3entry_xmax = 0;
#line 46 "PSI_Focus.instr"
  mccVacuumTube3entry_ymin = 0;
#line 46 "PSI_Focus.instr"
  mccVacuumTube3entry_ymax = 0;
#line 198 "PSI_Focus.instr"
  mccVacuumTube3entry_radius = 0.1335;
#line 46 "PSI_Focus.instr"
  mccVacuumTube3entry_xwidth = 0;
#line 46 "PSI_Focus.instr"
  mccVacuumTube3entry_yheight = 0;
#line 14086 "PSI_Focus.c"

  SIG_MESSAGE("VacuumTube3entry (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 14093 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa1, mcrotaVacuumTube3entry);
  rot_transpose(mcrotaVacuumTube2exit, mctr1);
  rot_mul(mcrotaVacuumTube3entry, mctr1, mcrotrVacuumTube3entry);
  mctc1 = coords_set(
#line 199 "PSI_Focus.instr"
    0,
#line 199 "PSI_Focus.instr"
    0,
#line 199 "PSI_Focus.instr"
    71.402);
#line 14104 "PSI_Focus.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaVacuumTube3entry = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaVacuumTube2exit, mcposaVacuumTube3entry);
  mcposrVacuumTube3entry = rot_apply(mcrotaVacuumTube3entry, mctc1);
  mcDEBUG_COMPONENT("VacuumTube3entry", mcposaVacuumTube3entry, mcrotaVacuumTube3entry)
  mccomp_posa[18] = mcposaVacuumTube3entry;
  mccomp_posr[18] = mcposrVacuumTube3entry;
  mcNCounter[18]  = mcPCounter[18] = mcP2Counter[18] = 0;
  mcAbsorbProp[18]= 0;
    /* Component VacuumTube3exit. */
  /* Setting parameters for component VacuumTube3exit. */
  SIG_MESSAGE("VacuumTube3exit (Init:SetPar)");
#line 46 "PSI_Focus.instr"
  mccVacuumTube3exit_xmin = 0;
#line 46 "PSI_Focus.instr"
  mccVacuumTube3exit_xmax = 0;
#line 46 "PSI_Focus.instr"
  mccVacuumTube3exit_ymin = 0;
#line 46 "PSI_Focus.instr"
  mccVacuumTube3exit_ymax = 0;
#line 201 "PSI_Focus.instr"
  mccVacuumTube3exit_radius = 0.1335;
#line 46 "PSI_Focus.instr"
  mccVacuumTube3exit_xwidth = 0;
#line 46 "PSI_Focus.instr"
  mccVacuumTube3exit_yheight = 0;
#line 14132 "PSI_Focus.c"

  SIG_MESSAGE("VacuumTube3exit (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 14139 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa1, mcrotaVacuumTube3exit);
  rot_transpose(mcrotaVacuumTube3entry, mctr1);
  rot_mul(mcrotaVacuumTube3exit, mctr1, mcrotrVacuumTube3exit);
  mctc1 = coords_set(
#line 202 "PSI_Focus.instr"
    0,
#line 202 "PSI_Focus.instr"
    0,
#line 202 "PSI_Focus.instr"
    72.007);
#line 14150 "PSI_Focus.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaVacuumTube3exit = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaVacuumTube3entry, mcposaVacuumTube3exit);
  mcposrVacuumTube3exit = rot_apply(mcrotaVacuumTube3exit, mctc1);
  mcDEBUG_COMPONENT("VacuumTube3exit", mcposaVacuumTube3exit, mcrotaVacuumTube3exit)
  mccomp_posa[19] = mcposaVacuumTube3exit;
  mccomp_posr[19] = mcposrVacuumTube3exit;
  mcNCounter[19]  = mcPCounter[19] = mcP2Counter[19] = 0;
  mcAbsorbProp[19]= 0;
    /* Component PSDmonMono. */
  /* Setting parameters for component PSDmonMono. */
  SIG_MESSAGE("PSDmonMono (Init:SetPar)");
#line 208 "PSI_Focus.instr"
  if("psdmono.dat") strncpy(mccPSDmonMono_filename, "psdmono.dat" ? "psdmono.dat" : "", 16384); else mccPSDmonMono_filename[0]='\0';
#line 50 "PSI_Focus.instr"
  mccPSDmonMono_xmin = -0.05;
#line 50 "PSI_Focus.instr"
  mccPSDmonMono_xmax = 0.05;
#line 50 "PSI_Focus.instr"
  mccPSDmonMono_ymin = -0.05;
#line 50 "PSI_Focus.instr"
  mccPSDmonMono_ymax = 0.05;
#line 207 "PSI_Focus.instr"
  mccPSDmonMono_xwidth = 0.1;
#line 207 "PSI_Focus.instr"
  mccPSDmonMono_yheight = 0.2;
#line 50 "PSI_Focus.instr"
  mccPSDmonMono_restore_neutron = 0;
#line 50 "PSI_Focus.instr"
  mccPSDmonMono_nowritefile = 0;
#line 14182 "PSI_Focus.c"

  SIG_MESSAGE("PSDmonMono (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 209 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 209 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 209 "PSI_Focus.instr"
    (0)*DEG2RAD);
#line 14192 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa1, mcrotaPSDmonMono);
  rot_transpose(mcrotaVacuumTube3exit, mctr1);
  rot_mul(mcrotaPSDmonMono, mctr1, mcrotrPSDmonMono);
  mctc1 = coords_set(
#line 209 "PSI_Focus.instr"
    0,
#line 209 "PSI_Focus.instr"
    0,
#line 209 "PSI_Focus.instr"
    72.4);
#line 14203 "PSI_Focus.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPSDmonMono = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaVacuumTube3exit, mcposaPSDmonMono);
  mcposrPSDmonMono = rot_apply(mcrotaPSDmonMono, mctc1);
  mcDEBUG_COMPONENT("PSDmonMono", mcposaPSDmonMono, mcrotaPSDmonMono)
  mccomp_posa[20] = mcposaPSDmonMono;
  mccomp_posr[20] = mcposrPSDmonMono;
  mcNCounter[20]  = mcPCounter[20] = mcP2Counter[20] = 0;
  mcAbsorbProp[20]= 0;
    /* Component MONOTOF. */
  /* Setting parameters for component MONOTOF. */
  SIG_MESSAGE("MONOTOF (Init:SetPar)");
#line 201 "PSI_Focus.instr"
  mccMONOTOF_xwidth = 0;
#line 201 "PSI_Focus.instr"
  mccMONOTOF_yheight = 0;
#line 201 "PSI_Focus.instr"
  mccMONOTOF_zdepth = 0;
#line 212 "PSI_Focus.instr"
  mccMONOTOF_xmin = -0.05;
#line 213 "PSI_Focus.instr"
  mccMONOTOF_xmax = 0.05;
#line 213 "PSI_Focus.instr"
  mccMONOTOF_ymin = -0.1;
#line 213 "PSI_Focus.instr"
  mccMONOTOF_ymax = 0.1;
#line 202 "PSI_Focus.instr"
  mccMONOTOF_zmin = 0;
#line 202 "PSI_Focus.instr"
  mccMONOTOF_zmax = 0;
#line 212 "PSI_Focus.instr"
  mccMONOTOF_bins = 30;
#line 203 "PSI_Focus.instr"
  mccMONOTOF_min = -1e40;
#line 203 "PSI_Focus.instr"
  mccMONOTOF_max = 1e40;
#line 203 "PSI_Focus.instr"
  mccMONOTOF_restore_neutron = 0;
#line 203 "PSI_Focus.instr"
  mccMONOTOF_radius = 0;
#line 213 "PSI_Focus.instr"
  if("auto, time") strncpy(mccMONOTOF_options, "auto, time" ? "auto, time" : "", 16384); else mccMONOTOF_options[0]='\0';
#line 212 "PSI_Focus.instr"
  if("MONO_TOF.dat") strncpy(mccMONOTOF_filename, "MONO_TOF.dat" ? "MONO_TOF.dat" : "", 16384); else mccMONOTOF_filename[0]='\0';
#line 204 "PSI_Focus.instr"
  if("NULL") strncpy(mccMONOTOF_geometry, "NULL" ? "NULL" : "", 16384); else mccMONOTOF_geometry[0]='\0';
#line 205 "PSI_Focus.instr"
  if("NULL") strncpy(mccMONOTOF_username1, "NULL" ? "NULL" : "", 16384); else mccMONOTOF_username1[0]='\0';
#line 205 "PSI_Focus.instr"
  if("NULL") strncpy(mccMONOTOF_username2, "NULL" ? "NULL" : "", 16384); else mccMONOTOF_username2[0]='\0';
#line 205 "PSI_Focus.instr"
  if("NULL") strncpy(mccMONOTOF_username3, "NULL" ? "NULL" : "", 16384); else mccMONOTOF_username3[0]='\0';
#line 206 "PSI_Focus.instr"
  mccMONOTOF_nowritefile = 0;
#line 14259 "PSI_Focus.c"

  SIG_MESSAGE("MONOTOF (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 14266 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa1, mcrotaMONOTOF);
  rot_transpose(mcrotaPSDmonMono, mctr1);
  rot_mul(mcrotaMONOTOF, mctr1, mcrotrMONOTOF);
  mctc1 = coords_set(
#line 214 "PSI_Focus.instr"
    0,
#line 214 "PSI_Focus.instr"
    0,
#line 214 "PSI_Focus.instr"
    72.41);
#line 14277 "PSI_Focus.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMONOTOF = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaPSDmonMono, mcposaMONOTOF);
  mcposrMONOTOF = rot_apply(mcrotaMONOTOF, mctc1);
  mcDEBUG_COMPONENT("MONOTOF", mcposaMONOTOF, mcrotaMONOTOF)
  mccomp_posa[21] = mcposaMONOTOF;
  mccomp_posr[21] = mcposrMONOTOF;
  mcNCounter[21]  = mcPCounter[21] = mcP2Counter[21] = 0;
  mcAbsorbProp[21]= 0;
    /* Component DivMonMono. */
  /* Setting parameters for component DivMonMono. */
  SIG_MESSAGE("DivMonMono (Init:SetPar)");
#line 217 "PSI_Focus.instr"
  if("divmono.dat") strncpy(mccDivMonMono_filename, "divmono.dat" ? "divmono.dat" : "", 16384); else mccDivMonMono_filename[0]='\0';
#line 55 "PSI_Focus.instr"
  mccDivMonMono_xmin = -0.05;
#line 55 "PSI_Focus.instr"
  mccDivMonMono_xmax = 0.05;
#line 55 "PSI_Focus.instr"
  mccDivMonMono_ymin = -0.05;
#line 55 "PSI_Focus.instr"
  mccDivMonMono_ymax = 0.05;
#line 218 "PSI_Focus.instr"
  mccDivMonMono_xwidth = 0.1;
#line 218 "PSI_Focus.instr"
  mccDivMonMono_yheight = 0.2;
#line 219 "PSI_Focus.instr"
  mccDivMonMono_maxdiv_h = 3;
#line 219 "PSI_Focus.instr"
  mccDivMonMono_maxdiv_v = 3;
#line 56 "PSI_Focus.instr"
  mccDivMonMono_restore_neutron = 0;
#line 56 "PSI_Focus.instr"
  mccDivMonMono_nx = 0;
#line 56 "PSI_Focus.instr"
  mccDivMonMono_ny = 0;
#line 56 "PSI_Focus.instr"
  mccDivMonMono_nz = 1;
#line 56 "PSI_Focus.instr"
  mccDivMonMono_nowritefile = 0;
#line 14319 "PSI_Focus.c"

  SIG_MESSAGE("DivMonMono (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 220 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 220 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 220 "PSI_Focus.instr"
    (0)*DEG2RAD);
#line 14329 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa1, mcrotaDivMonMono);
  rot_transpose(mcrotaMONOTOF, mctr1);
  rot_mul(mcrotaDivMonMono, mctr1, mcrotrDivMonMono);
  mctc1 = coords_set(
#line 220 "PSI_Focus.instr"
    0,
#line 220 "PSI_Focus.instr"
    0,
#line 220 "PSI_Focus.instr"
    72.42);
#line 14340 "PSI_Focus.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaDivMonMono = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaMONOTOF, mcposaDivMonMono);
  mcposrDivMonMono = rot_apply(mcrotaDivMonMono, mctc1);
  mcDEBUG_COMPONENT("DivMonMono", mcposaDivMonMono, mcrotaDivMonMono)
  mccomp_posa[22] = mcposaDivMonMono;
  mccomp_posr[22] = mcposrDivMonMono;
  mcNCounter[22]  = mcPCounter[22] = mcP2Counter[22] = 0;
  mcAbsorbProp[22]= 0;
    /* Component focus_mono. */
  /* Setting parameters for component focus_mono. */
  SIG_MESSAGE("focus_mono (Init:SetPar)");

  SIG_MESSAGE("focus_mono (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 224 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 224 "PSI_Focus.instr"
    (PHM)*DEG2RAD,
#line 224 "PSI_Focus.instr"
    (0)*DEG2RAD);
#line 14363 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa1, mcrotafocus_mono);
  rot_transpose(mcrotaDivMonMono, mctr1);
  rot_mul(mcrotafocus_mono, mctr1, mcrotrfocus_mono);
  mctc1 = coords_set(
#line 224 "PSI_Focus.instr"
    0,
#line 224 "PSI_Focus.instr"
    0,
#line 224 "PSI_Focus.instr"
    72.617);
#line 14374 "PSI_Focus.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposafocus_mono = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaDivMonMono, mcposafocus_mono);
  mcposrfocus_mono = rot_apply(mcrotafocus_mono, mctc1);
  mcDEBUG_COMPONENT("focus_mono", mcposafocus_mono, mcrotafocus_mono)
  mccomp_posa[23] = mcposafocus_mono;
  mccomp_posr[23] = mcposrfocus_mono;
  mcNCounter[23]  = mcPCounter[23] = mcP2Counter[23] = 0;
  mcAbsorbProp[23]= 0;
    /* Component mono. */
  /* Setting parameters for component mono. */
  SIG_MESSAGE("mono (Init:SetPar)");
#line 78 "PSI_Focus.instr"
  if(0) strncpy(mccmono_reflect, 0 ? 0 : "", 16384); else mccmono_reflect[0]='\0';
#line 229 "PSI_Focus.instr"
  mccmono_zwidth = 0.019;
#line 229 "PSI_Focus.instr"
  mccmono_yheight = 0.025;
#line 229 "PSI_Focus.instr"
  mccmono_gap = 0.001;
#line 230 "PSI_Focus.instr"
  mccmono_NH = 9;
#line 230 "PSI_Focus.instr"
  mccmono_NV = 7;
#line 231 "PSI_Focus.instr"
  mccmono_mosaich = 48;
#line 231 "PSI_Focus.instr"
  mccmono_mosaicv = 48;
#line 232 "PSI_Focus.instr"
  mccmono_r0 = 0.99;
#line 233 "PSI_Focus.instr"
  mccmono_Q = 1.873;
#line 234 "PSI_Focus.instr"
  mccmono_RV = RV_2;
#line 234 "PSI_Focus.instr"
  mccmono_RH = RH_2;
#line 78 "PSI_Focus.instr"
  mccmono_DM = 0;
#line 78 "PSI_Focus.instr"
  mccmono_mosaic = 0;
#line 78 "PSI_Focus.instr"
  mccmono_width = 0;
#line 78 "PSI_Focus.instr"
  mccmono_height = 0;
#line 78 "PSI_Focus.instr"
  mccmono_verbose = 0;
#line 14422 "PSI_Focus.c"

  SIG_MESSAGE("mono (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 14429 "PSI_Focus.c"
  rot_mul(mctr1, mcrotafocus_mono, mcrotamono);
  rot_transpose(mcrotafocus_mono, mctr1);
  rot_mul(mcrotamono, mctr1, mcrotrmono);
  mctc1 = coords_set(
#line 235 "PSI_Focus.instr"
    0,
#line 235 "PSI_Focus.instr"
    0,
#line 235 "PSI_Focus.instr"
    0);
#line 14440 "PSI_Focus.c"
  rot_transpose(mcrotafocus_mono, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposamono = coords_add(mcposafocus_mono, mctc2);
  mctc1 = coords_sub(mcposafocus_mono, mcposamono);
  mcposrmono = rot_apply(mcrotamono, mctc1);
  mcDEBUG_COMPONENT("mono", mcposamono, mcrotamono)
  mccomp_posa[24] = mcposamono;
  mccomp_posr[24] = mcposrmono;
  mcNCounter[24]  = mcPCounter[24] = mcP2Counter[24] = 0;
  mcAbsorbProp[24]= 0;
    /* Component a2. */
  /* Setting parameters for component a2. */
  SIG_MESSAGE("a2 (Init:SetPar)");

  SIG_MESSAGE("a2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 238 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 238 "PSI_Focus.instr"
    (TTM)*DEG2RAD,
#line 238 "PSI_Focus.instr"
    (0)*DEG2RAD);
#line 14463 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa1, mcrotaa2);
  rot_transpose(mcrotamono, mctr1);
  rot_mul(mcrotaa2, mctr1, mcrotra2);
  mctc1 = coords_set(
#line 238 "PSI_Focus.instr"
    0,
#line 238 "PSI_Focus.instr"
    0,
#line 238 "PSI_Focus.instr"
    0);
#line 14474 "PSI_Focus.c"
  rot_transpose(mcrotafocus_mono, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaa2 = coords_add(mcposafocus_mono, mctc2);
  mctc1 = coords_sub(mcposamono, mcposaa2);
  mcposra2 = rot_apply(mcrotaa2, mctc1);
  mcDEBUG_COMPONENT("a2", mcposaa2, mcrotaa2)
  mccomp_posa[25] = mcposaa2;
  mccomp_posr[25] = mcposra2;
  mcNCounter[25]  = mcPCounter[25] = mcP2Counter[25] = 0;
  mcAbsorbProp[25]= 0;
    /* Component FERMITOF_before. */
  /* Setting parameters for component FERMITOF_before. */
  SIG_MESSAGE("FERMITOF_before (Init:SetPar)");
#line 201 "PSI_Focus.instr"
  mccFERMITOF_before_xwidth = 0;
#line 201 "PSI_Focus.instr"
  mccFERMITOF_before_yheight = 0;
#line 201 "PSI_Focus.instr"
  mccFERMITOF_before_zdepth = 0;
#line 241 "PSI_Focus.instr"
  mccFERMITOF_before_xmin = -0.05;
#line 242 "PSI_Focus.instr"
  mccFERMITOF_before_xmax = 0.05;
#line 242 "PSI_Focus.instr"
  mccFERMITOF_before_ymin = -0.1;
#line 242 "PSI_Focus.instr"
  mccFERMITOF_before_ymax = 0.1;
#line 202 "PSI_Focus.instr"
  mccFERMITOF_before_zmin = 0;
#line 202 "PSI_Focus.instr"
  mccFERMITOF_before_zmax = 0;
#line 241 "PSI_Focus.instr"
  mccFERMITOF_before_bins = 30;
#line 203 "PSI_Focus.instr"
  mccFERMITOF_before_min = -1e40;
#line 203 "PSI_Focus.instr"
  mccFERMITOF_before_max = 1e40;
#line 203 "PSI_Focus.instr"
  mccFERMITOF_before_restore_neutron = 0;
#line 203 "PSI_Focus.instr"
  mccFERMITOF_before_radius = 0;
#line 242 "PSI_Focus.instr"
  if("auto, time") strncpy(mccFERMITOF_before_options, "auto, time" ? "auto, time" : "", 16384); else mccFERMITOF_before_options[0]='\0';
#line 241 "PSI_Focus.instr"
  if("FERMI_TOF_before.dat") strncpy(mccFERMITOF_before_filename, "FERMI_TOF_before.dat" ? "FERMI_TOF_before.dat" : "", 16384); else mccFERMITOF_before_filename[0]='\0';
#line 204 "PSI_Focus.instr"
  if("NULL") strncpy(mccFERMITOF_before_geometry, "NULL" ? "NULL" : "", 16384); else mccFERMITOF_before_geometry[0]='\0';
#line 205 "PSI_Focus.instr"
  if("NULL") strncpy(mccFERMITOF_before_username1, "NULL" ? "NULL" : "", 16384); else mccFERMITOF_before_username1[0]='\0';
#line 205 "PSI_Focus.instr"
  if("NULL") strncpy(mccFERMITOF_before_username2, "NULL" ? "NULL" : "", 16384); else mccFERMITOF_before_username2[0]='\0';
#line 205 "PSI_Focus.instr"
  if("NULL") strncpy(mccFERMITOF_before_username3, "NULL" ? "NULL" : "", 16384); else mccFERMITOF_before_username3[0]='\0';
#line 206 "PSI_Focus.instr"
  mccFERMITOF_before_nowritefile = 0;
#line 14530 "PSI_Focus.c"

  SIG_MESSAGE("FERMITOF_before (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 14537 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa2, mcrotaFERMITOF_before);
  rot_transpose(mcrotaa2, mctr1);
  rot_mul(mcrotaFERMITOF_before, mctr1, mcrotrFERMITOF_before);
  mctc1 = coords_set(
#line 243 "PSI_Focus.instr"
    0,
#line 243 "PSI_Focus.instr"
    0,
#line 243 "PSI_Focus.instr"
    0.9);
#line 14548 "PSI_Focus.c"
  rot_transpose(mcrotaa2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaFERMITOF_before = coords_add(mcposaa2, mctc2);
  mctc1 = coords_sub(mcposaa2, mcposaFERMITOF_before);
  mcposrFERMITOF_before = rot_apply(mcrotaFERMITOF_before, mctc1);
  mcDEBUG_COMPONENT("FERMITOF_before", mcposaFERMITOF_before, mcrotaFERMITOF_before)
  mccomp_posa[26] = mcposaFERMITOF_before;
  mccomp_posr[26] = mcposrFERMITOF_before;
  mcNCounter[26]  = mcPCounter[26] = mcP2Counter[26] = 0;
  mcAbsorbProp[26]= 0;
    /* Component lambdaFermi. */
  /* Setting parameters for component lambdaFermi. */
  SIG_MESSAGE("lambdaFermi (Init:SetPar)");
#line 248 "PSI_Focus.instr"
  if("lambdafermi.dat") strncpy(mcclambdaFermi_filename, "lambdafermi.dat" ? "lambdafermi.dat" : "", 16384); else mcclambdaFermi_filename[0]='\0';
#line 50 "PSI_Focus.instr"
  mcclambdaFermi_xmin = -0.05;
#line 50 "PSI_Focus.instr"
  mcclambdaFermi_xmax = 0.05;
#line 50 "PSI_Focus.instr"
  mcclambdaFermi_ymin = -0.05;
#line 50 "PSI_Focus.instr"
  mcclambdaFermi_ymax = 0.05;
#line 246 "PSI_Focus.instr"
  mcclambdaFermi_xwidth = 0.1;
#line 246 "PSI_Focus.instr"
  mcclambdaFermi_yheight = 0.2;
#line 247 "PSI_Focus.instr"
  mcclambdaFermi_Lmin = LMIN;
#line 247 "PSI_Focus.instr"
  mcclambdaFermi_Lmax = LMAX;
#line 51 "PSI_Focus.instr"
  mcclambdaFermi_restore_neutron = 0;
#line 51 "PSI_Focus.instr"
  mcclambdaFermi_nowritefile = 0;
#line 14584 "PSI_Focus.c"

  SIG_MESSAGE("lambdaFermi (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 14591 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa2, mcrotalambdaFermi);
  rot_transpose(mcrotaFERMITOF_before, mctr1);
  rot_mul(mcrotalambdaFermi, mctr1, mcrotrlambdaFermi);
  mctc1 = coords_set(
#line 249 "PSI_Focus.instr"
    0,
#line 249 "PSI_Focus.instr"
    0,
#line 249 "PSI_Focus.instr"
    0.901);
#line 14602 "PSI_Focus.c"
  rot_transpose(mcrotaa2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposalambdaFermi = coords_add(mcposaa2, mctc2);
  mctc1 = coords_sub(mcposaFERMITOF_before, mcposalambdaFermi);
  mcposrlambdaFermi = rot_apply(mcrotalambdaFermi, mctc1);
  mcDEBUG_COMPONENT("lambdaFermi", mcposalambdaFermi, mcrotalambdaFermi)
  mccomp_posa[27] = mcposalambdaFermi;
  mccomp_posr[27] = mcposrlambdaFermi;
  mcNCounter[27]  = mcPCounter[27] = mcP2Counter[27] = 0;
  mcAbsorbProp[27]= 0;
    /* Component EMON_Fermi. */
  /* Setting parameters for component EMON_Fermi. */
  SIG_MESSAGE("EMON_Fermi (Init:SetPar)");
#line 254 "PSI_Focus.instr"
  if("emon_fermi.dat") strncpy(mccEMON_Fermi_filename, "emon_fermi.dat" ? "emon_fermi.dat" : "", 16384); else mccEMON_Fermi_filename[0]='\0';
#line 53 "PSI_Focus.instr"
  mccEMON_Fermi_xmin = -0.05;
#line 53 "PSI_Focus.instr"
  mccEMON_Fermi_xmax = 0.05;
#line 53 "PSI_Focus.instr"
  mccEMON_Fermi_ymin = -0.05;
#line 53 "PSI_Focus.instr"
  mccEMON_Fermi_ymax = 0.05;
#line 252 "PSI_Focus.instr"
  mccEMON_Fermi_xwidth = 0.06;
#line 252 "PSI_Focus.instr"
  mccEMON_Fermi_yheight = 0.1;
#line 253 "PSI_Focus.instr"
  mccEMON_Fermi_Emin = EMIN;
#line 253 "PSI_Focus.instr"
  mccEMON_Fermi_Emax = EMAX;
#line 54 "PSI_Focus.instr"
  mccEMON_Fermi_restore_neutron = 0;
#line 54 "PSI_Focus.instr"
  mccEMON_Fermi_nowritefile = 0;
#line 14638 "PSI_Focus.c"

  SIG_MESSAGE("EMON_Fermi (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 255 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 255 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 255 "PSI_Focus.instr"
    (0)*DEG2RAD);
#line 14648 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa2, mcrotaEMON_Fermi);
  rot_transpose(mcrotalambdaFermi, mctr1);
  rot_mul(mcrotaEMON_Fermi, mctr1, mcrotrEMON_Fermi);
  mctc1 = coords_set(
#line 255 "PSI_Focus.instr"
    0,
#line 255 "PSI_Focus.instr"
    0,
#line 255 "PSI_Focus.instr"
    0.9397);
#line 14659 "PSI_Focus.c"
  rot_transpose(mcrotaa2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaEMON_Fermi = coords_add(mcposaa2, mctc2);
  mctc1 = coords_sub(mcposalambdaFermi, mcposaEMON_Fermi);
  mcposrEMON_Fermi = rot_apply(mcrotaEMON_Fermi, mctc1);
  mcDEBUG_COMPONENT("EMON_Fermi", mcposaEMON_Fermi, mcrotaEMON_Fermi)
  mccomp_posa[28] = mcposaEMON_Fermi;
  mccomp_posr[28] = mcposrEMON_Fermi;
  mcNCounter[28]  = mcPCounter[28] = mcP2Counter[28] = 0;
  mcAbsorbProp[28]= 0;
    /* Component DivMonfermi1. */
  /* Setting parameters for component DivMonfermi1. */
  SIG_MESSAGE("DivMonfermi1 (Init:SetPar)");
#line 258 "PSI_Focus.instr"
  if("divfermi1.dat") strncpy(mccDivMonfermi1_filename, "divfermi1.dat" ? "divfermi1.dat" : "", 16384); else mccDivMonfermi1_filename[0]='\0';
#line 55 "PSI_Focus.instr"
  mccDivMonfermi1_xmin = -0.05;
#line 55 "PSI_Focus.instr"
  mccDivMonfermi1_xmax = 0.05;
#line 55 "PSI_Focus.instr"
  mccDivMonfermi1_ymin = -0.05;
#line 55 "PSI_Focus.instr"
  mccDivMonfermi1_ymax = 0.05;
#line 259 "PSI_Focus.instr"
  mccDivMonfermi1_xwidth = 0.06;
#line 259 "PSI_Focus.instr"
  mccDivMonfermi1_yheight = 0.1;
#line 260 "PSI_Focus.instr"
  mccDivMonfermi1_maxdiv_h = 2;
#line 260 "PSI_Focus.instr"
  mccDivMonfermi1_maxdiv_v = 2;
#line 56 "PSI_Focus.instr"
  mccDivMonfermi1_restore_neutron = 0;
#line 56 "PSI_Focus.instr"
  mccDivMonfermi1_nx = 0;
#line 56 "PSI_Focus.instr"
  mccDivMonfermi1_ny = 0;
#line 56 "PSI_Focus.instr"
  mccDivMonfermi1_nz = 1;
#line 56 "PSI_Focus.instr"
  mccDivMonfermi1_nowritefile = 0;
#line 14701 "PSI_Focus.c"

  SIG_MESSAGE("DivMonfermi1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 14708 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa2, mcrotaDivMonfermi1);
  rot_transpose(mcrotaEMON_Fermi, mctr1);
  rot_mul(mcrotaDivMonfermi1, mctr1, mcrotrDivMonfermi1);
  mctc1 = coords_set(
#line 261 "PSI_Focus.instr"
    0,
#line 261 "PSI_Focus.instr"
    0,
#line 261 "PSI_Focus.instr"
    0.9398);
#line 14719 "PSI_Focus.c"
  rot_transpose(mcrotaa2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaDivMonfermi1 = coords_add(mcposaa2, mctc2);
  mctc1 = coords_sub(mcposaEMON_Fermi, mcposaDivMonfermi1);
  mcposrDivMonfermi1 = rot_apply(mcrotaDivMonfermi1, mctc1);
  mcDEBUG_COMPONENT("DivMonfermi1", mcposaDivMonfermi1, mcrotaDivMonfermi1)
  mccomp_posa[29] = mcposaDivMonfermi1;
  mccomp_posr[29] = mcposrDivMonfermi1;
  mcNCounter[29]  = mcPCounter[29] = mcP2Counter[29] = 0;
  mcAbsorbProp[29]= 0;
    /* Component PSD_Fermi1. */
  /* Setting parameters for component PSD_Fermi1. */
  SIG_MESSAGE("PSD_Fermi1 (Init:SetPar)");
#line 265 "PSI_Focus.instr"
  if("psdfermi1.dat") strncpy(mccPSD_Fermi1_filename, "psdfermi1.dat" ? "psdfermi1.dat" : "", 16384); else mccPSD_Fermi1_filename[0]='\0';
#line 50 "PSI_Focus.instr"
  mccPSD_Fermi1_xmin = -0.05;
#line 50 "PSI_Focus.instr"
  mccPSD_Fermi1_xmax = 0.05;
#line 50 "PSI_Focus.instr"
  mccPSD_Fermi1_ymin = -0.05;
#line 50 "PSI_Focus.instr"
  mccPSD_Fermi1_ymax = 0.05;
#line 264 "PSI_Focus.instr"
  mccPSD_Fermi1_xwidth = 0.06;
#line 264 "PSI_Focus.instr"
  mccPSD_Fermi1_yheight = 0.11;
#line 50 "PSI_Focus.instr"
  mccPSD_Fermi1_restore_neutron = 0;
#line 50 "PSI_Focus.instr"
  mccPSD_Fermi1_nowritefile = 0;
#line 14751 "PSI_Focus.c"

  SIG_MESSAGE("PSD_Fermi1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 266 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 266 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 266 "PSI_Focus.instr"
    (0)*DEG2RAD);
#line 14761 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa2, mcrotaPSD_Fermi1);
  rot_transpose(mcrotaDivMonfermi1, mctr1);
  rot_mul(mcrotaPSD_Fermi1, mctr1, mcrotrPSD_Fermi1);
  mctc1 = coords_set(
#line 266 "PSI_Focus.instr"
    0,
#line 266 "PSI_Focus.instr"
    0,
#line 266 "PSI_Focus.instr"
    0.9399);
#line 14772 "PSI_Focus.c"
  rot_transpose(mcrotaa2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPSD_Fermi1 = coords_add(mcposaa2, mctc2);
  mctc1 = coords_sub(mcposaDivMonfermi1, mcposaPSD_Fermi1);
  mcposrPSD_Fermi1 = rot_apply(mcrotaPSD_Fermi1, mctc1);
  mcDEBUG_COMPONENT("PSD_Fermi1", mcposaPSD_Fermi1, mcrotaPSD_Fermi1)
  mccomp_posa[30] = mcposaPSD_Fermi1;
  mccomp_posr[30] = mcposrPSD_Fermi1;
  mcNCounter[30]  = mcPCounter[30] = mcP2Counter[30] = 0;
  mcAbsorbProp[30]= 0;
    /* Component FoChopper. */
  /* Setting parameters for component FoChopper. */
  SIG_MESSAGE("FoChopper (Init:SetPar)");
#line 88 "PSI_Focus.instr"
  mccFoChopper_phase = 0;
#line 269 "PSI_Focus.instr"
  mccFoChopper_radius = 0.06;
#line 269 "PSI_Focus.instr"
  mccFoChopper_nu = - FERMI_SPEED;
#line 270 "PSI_Focus.instr"
  mccFoChopper_w = 0.0005;
#line 270 "PSI_Focus.instr"
  mccFoChopper_nslit = 120;
#line 271 "PSI_Focus.instr"
  mccFoChopper_R0 = 0;
#line 271 "PSI_Focus.instr"
  mccFoChopper_Qc = 0;
#line 270 "PSI_Focus.instr"
  mccFoChopper_alpha = 0;
#line 271 "PSI_Focus.instr"
  mccFoChopper_m = 0;
#line 271 "PSI_Focus.instr"
  mccFoChopper_W = 0.0001;
#line 271 "PSI_Focus.instr"
  mccFoChopper_length = 0.012;
#line 90 "PSI_Focus.instr"
  mccFoChopper_eff = 0.95;
#line 91 "PSI_Focus.instr"
  mccFoChopper_zero_time = 0;
#line 91 "PSI_Focus.instr"
  mccFoChopper_xwidth = 0;
#line 91 "PSI_Focus.instr"
  mccFoChopper_verbose = 0;
#line 269 "PSI_Focus.instr"
  mccFoChopper_yheight = 0.11;
#line 92 "PSI_Focus.instr"
  mccFoChopper_curvature = 0;
#line 92 "PSI_Focus.instr"
  mccFoChopper_delay = 0;
#line 14822 "PSI_Focus.c"

  SIG_MESSAGE("FoChopper (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 14829 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa2, mcrotaFoChopper);
  rot_transpose(mcrotaPSD_Fermi1, mctr1);
  rot_mul(mcrotaFoChopper, mctr1, mcrotrFoChopper);
  mctc1 = coords_set(
#line 272 "PSI_Focus.instr"
    0,
#line 272 "PSI_Focus.instr"
    0,
#line 272 "PSI_Focus.instr"
    1.002);
#line 14840 "PSI_Focus.c"
  rot_transpose(mcrotaa2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaFoChopper = coords_add(mcposaa2, mctc2);
  mctc1 = coords_sub(mcposaPSD_Fermi1, mcposaFoChopper);
  mcposrFoChopper = rot_apply(mcrotaFoChopper, mctc1);
  mcDEBUG_COMPONENT("FoChopper", mcposaFoChopper, mcrotaFoChopper)
  mccomp_posa[31] = mcposaFoChopper;
  mccomp_posr[31] = mcposrFoChopper;
  mcNCounter[31]  = mcPCounter[31] = mcP2Counter[31] = 0;
  mcAbsorbProp[31]= 0;
    /* Component PSD_Fermi2. */
  /* Setting parameters for component PSD_Fermi2. */
  SIG_MESSAGE("PSD_Fermi2 (Init:SetPar)");
#line 276 "PSI_Focus.instr"
  if("psdfermi2.dat") strncpy(mccPSD_Fermi2_filename, "psdfermi2.dat" ? "psdfermi2.dat" : "", 16384); else mccPSD_Fermi2_filename[0]='\0';
#line 50 "PSI_Focus.instr"
  mccPSD_Fermi2_xmin = -0.05;
#line 50 "PSI_Focus.instr"
  mccPSD_Fermi2_xmax = 0.05;
#line 50 "PSI_Focus.instr"
  mccPSD_Fermi2_ymin = -0.05;
#line 50 "PSI_Focus.instr"
  mccPSD_Fermi2_ymax = 0.05;
#line 275 "PSI_Focus.instr"
  mccPSD_Fermi2_xwidth = 0.06;
#line 275 "PSI_Focus.instr"
  mccPSD_Fermi2_yheight = 0.11;
#line 50 "PSI_Focus.instr"
  mccPSD_Fermi2_restore_neutron = 0;
#line 50 "PSI_Focus.instr"
  mccPSD_Fermi2_nowritefile = 0;
#line 14872 "PSI_Focus.c"

  SIG_MESSAGE("PSD_Fermi2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 277 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 277 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 277 "PSI_Focus.instr"
    (0)*DEG2RAD);
#line 14882 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa2, mcrotaPSD_Fermi2);
  rot_transpose(mcrotaFoChopper, mctr1);
  rot_mul(mcrotaPSD_Fermi2, mctr1, mcrotrPSD_Fermi2);
  mctc1 = coords_set(
#line 277 "PSI_Focus.instr"
    0,
#line 277 "PSI_Focus.instr"
    0,
#line 277 "PSI_Focus.instr"
    1.063);
#line 14893 "PSI_Focus.c"
  rot_transpose(mcrotaa2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPSD_Fermi2 = coords_add(mcposaa2, mctc2);
  mctc1 = coords_sub(mcposaFoChopper, mcposaPSD_Fermi2);
  mcposrPSD_Fermi2 = rot_apply(mcrotaPSD_Fermi2, mctc1);
  mcDEBUG_COMPONENT("PSD_Fermi2", mcposaPSD_Fermi2, mcrotaPSD_Fermi2)
  mccomp_posa[32] = mcposaPSD_Fermi2;
  mccomp_posr[32] = mcposrPSD_Fermi2;
  mcNCounter[32]  = mcPCounter[32] = mcP2Counter[32] = 0;
  mcAbsorbProp[32]= 0;
    /* Component DivMonfermi2. */
  /* Setting parameters for component DivMonfermi2. */
  SIG_MESSAGE("DivMonfermi2 (Init:SetPar)");
#line 280 "PSI_Focus.instr"
  if("divfermi2.dat") strncpy(mccDivMonfermi2_filename, "divfermi2.dat" ? "divfermi2.dat" : "", 16384); else mccDivMonfermi2_filename[0]='\0';
#line 55 "PSI_Focus.instr"
  mccDivMonfermi2_xmin = -0.05;
#line 55 "PSI_Focus.instr"
  mccDivMonfermi2_xmax = 0.05;
#line 55 "PSI_Focus.instr"
  mccDivMonfermi2_ymin = -0.05;
#line 55 "PSI_Focus.instr"
  mccDivMonfermi2_ymax = 0.05;
#line 281 "PSI_Focus.instr"
  mccDivMonfermi2_xwidth = 0.06;
#line 281 "PSI_Focus.instr"
  mccDivMonfermi2_yheight = 0.1;
#line 282 "PSI_Focus.instr"
  mccDivMonfermi2_maxdiv_h = 2;
#line 282 "PSI_Focus.instr"
  mccDivMonfermi2_maxdiv_v = 2;
#line 56 "PSI_Focus.instr"
  mccDivMonfermi2_restore_neutron = 0;
#line 56 "PSI_Focus.instr"
  mccDivMonfermi2_nx = 0;
#line 56 "PSI_Focus.instr"
  mccDivMonfermi2_ny = 0;
#line 56 "PSI_Focus.instr"
  mccDivMonfermi2_nz = 1;
#line 56 "PSI_Focus.instr"
  mccDivMonfermi2_nowritefile = 0;
#line 14935 "PSI_Focus.c"

  SIG_MESSAGE("DivMonfermi2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 14942 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa2, mcrotaDivMonfermi2);
  rot_transpose(mcrotaPSD_Fermi2, mctr1);
  rot_mul(mcrotaDivMonfermi2, mctr1, mcrotrDivMonfermi2);
  mctc1 = coords_set(
#line 283 "PSI_Focus.instr"
    0,
#line 283 "PSI_Focus.instr"
    0,
#line 283 "PSI_Focus.instr"
    1.064);
#line 14953 "PSI_Focus.c"
  rot_transpose(mcrotaa2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaDivMonfermi2 = coords_add(mcposaa2, mctc2);
  mctc1 = coords_sub(mcposaPSD_Fermi2, mcposaDivMonfermi2);
  mcposrDivMonfermi2 = rot_apply(mcrotaDivMonfermi2, mctc1);
  mcDEBUG_COMPONENT("DivMonfermi2", mcposaDivMonfermi2, mcrotaDivMonfermi2)
  mccomp_posa[33] = mcposaDivMonfermi2;
  mccomp_posr[33] = mcposrDivMonfermi2;
  mcNCounter[33]  = mcPCounter[33] = mcP2Counter[33] = 0;
  mcAbsorbProp[33]= 0;
    /* Component FERMITOF1. */
  /* Setting parameters for component FERMITOF1. */
  SIG_MESSAGE("FERMITOF1 (Init:SetPar)");
#line 201 "PSI_Focus.instr"
  mccFERMITOF1_xwidth = 0;
#line 201 "PSI_Focus.instr"
  mccFERMITOF1_yheight = 0;
#line 201 "PSI_Focus.instr"
  mccFERMITOF1_zdepth = 0;
#line 286 "PSI_Focus.instr"
  mccFERMITOF1_xmin = -0.05;
#line 287 "PSI_Focus.instr"
  mccFERMITOF1_xmax = 0.05;
#line 287 "PSI_Focus.instr"
  mccFERMITOF1_ymin = -0.1;
#line 287 "PSI_Focus.instr"
  mccFERMITOF1_ymax = 0.1;
#line 202 "PSI_Focus.instr"
  mccFERMITOF1_zmin = 0;
#line 202 "PSI_Focus.instr"
  mccFERMITOF1_zmax = 0;
#line 286 "PSI_Focus.instr"
  mccFERMITOF1_bins = 30;
#line 203 "PSI_Focus.instr"
  mccFERMITOF1_min = -1e40;
#line 203 "PSI_Focus.instr"
  mccFERMITOF1_max = 1e40;
#line 203 "PSI_Focus.instr"
  mccFERMITOF1_restore_neutron = 0;
#line 203 "PSI_Focus.instr"
  mccFERMITOF1_radius = 0;
#line 287 "PSI_Focus.instr"
  if("auto, time") strncpy(mccFERMITOF1_options, "auto, time" ? "auto, time" : "", 16384); else mccFERMITOF1_options[0]='\0';
#line 286 "PSI_Focus.instr"
  if("FERMI_TOF1.dat.dat") strncpy(mccFERMITOF1_filename, "FERMI_TOF1.dat.dat" ? "FERMI_TOF1.dat.dat" : "", 16384); else mccFERMITOF1_filename[0]='\0';
#line 204 "PSI_Focus.instr"
  if("NULL") strncpy(mccFERMITOF1_geometry, "NULL" ? "NULL" : "", 16384); else mccFERMITOF1_geometry[0]='\0';
#line 205 "PSI_Focus.instr"
  if("NULL") strncpy(mccFERMITOF1_username1, "NULL" ? "NULL" : "", 16384); else mccFERMITOF1_username1[0]='\0';
#line 205 "PSI_Focus.instr"
  if("NULL") strncpy(mccFERMITOF1_username2, "NULL" ? "NULL" : "", 16384); else mccFERMITOF1_username2[0]='\0';
#line 205 "PSI_Focus.instr"
  if("NULL") strncpy(mccFERMITOF1_username3, "NULL" ? "NULL" : "", 16384); else mccFERMITOF1_username3[0]='\0';
#line 206 "PSI_Focus.instr"
  mccFERMITOF1_nowritefile = 0;
#line 15009 "PSI_Focus.c"

  SIG_MESSAGE("FERMITOF1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15016 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa2, mcrotaFERMITOF1);
  rot_transpose(mcrotaDivMonfermi2, mctr1);
  rot_mul(mcrotaFERMITOF1, mctr1, mcrotrFERMITOF1);
  mctc1 = coords_set(
#line 288 "PSI_Focus.instr"
    0,
#line 288 "PSI_Focus.instr"
    0,
#line 288 "PSI_Focus.instr"
    1.154);
#line 15027 "PSI_Focus.c"
  rot_transpose(mcrotaa2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaFERMITOF1 = coords_add(mcposaa2, mctc2);
  mctc1 = coords_sub(mcposaDivMonfermi2, mcposaFERMITOF1);
  mcposrFERMITOF1 = rot_apply(mcrotaFERMITOF1, mctc1);
  mcDEBUG_COMPONENT("FERMITOF1", mcposaFERMITOF1, mcrotaFERMITOF1)
  mccomp_posa[34] = mcposaFERMITOF1;
  mccomp_posr[34] = mcposrFERMITOF1;
  mcNCounter[34]  = mcPCounter[34] = mcP2Counter[34] = 0;
  mcAbsorbProp[34]= 0;
    /* Component SAMPLE_SLIT. */
  /* Setting parameters for component SAMPLE_SLIT. */
  SIG_MESSAGE("SAMPLE_SLIT (Init:SetPar)");
#line 46 "PSI_Focus.instr"
  mccSAMPLE_SLIT_xmin = 0;
#line 46 "PSI_Focus.instr"
  mccSAMPLE_SLIT_xmax = 0;
#line 46 "PSI_Focus.instr"
  mccSAMPLE_SLIT_ymin = 0;
#line 46 "PSI_Focus.instr"
  mccSAMPLE_SLIT_ymax = 0;
#line 46 "PSI_Focus.instr"
  mccSAMPLE_SLIT_radius = 0;
#line 291 "PSI_Focus.instr"
  mccSAMPLE_SLIT_xwidth = 0.02;
#line 291 "PSI_Focus.instr"
  mccSAMPLE_SLIT_yheight = 0.06;
#line 15055 "PSI_Focus.c"

  SIG_MESSAGE("SAMPLE_SLIT (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15062 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa2, mcrotaSAMPLE_SLIT);
  rot_transpose(mcrotaFERMITOF1, mctr1);
  rot_mul(mcrotaSAMPLE_SLIT, mctr1, mcrotrSAMPLE_SLIT);
  mctc1 = coords_set(
#line 292 "PSI_Focus.instr"
    0,
#line 292 "PSI_Focus.instr"
    0,
#line 292 "PSI_Focus.instr"
    1.155);
#line 15073 "PSI_Focus.c"
  rot_transpose(mcrotaa2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaSAMPLE_SLIT = coords_add(mcposaa2, mctc2);
  mctc1 = coords_sub(mcposaFERMITOF1, mcposaSAMPLE_SLIT);
  mcposrSAMPLE_SLIT = rot_apply(mcrotaSAMPLE_SLIT, mctc1);
  mcDEBUG_COMPONENT("SAMPLE_SLIT", mcposaSAMPLE_SLIT, mcrotaSAMPLE_SLIT)
  mccomp_posa[35] = mcposaSAMPLE_SLIT;
  mccomp_posr[35] = mcposrSAMPLE_SLIT;
  mcNCounter[35]  = mcPCounter[35] = mcP2Counter[35] = 0;
  mcAbsorbProp[35]= 0;
    /* Component FERMITOF2. */
  /* Setting parameters for component FERMITOF2. */
  SIG_MESSAGE("FERMITOF2 (Init:SetPar)");
#line 201 "PSI_Focus.instr"
  mccFERMITOF2_xwidth = 0;
#line 201 "PSI_Focus.instr"
  mccFERMITOF2_yheight = 0;
#line 201 "PSI_Focus.instr"
  mccFERMITOF2_zdepth = 0;
#line 295 "PSI_Focus.instr"
  mccFERMITOF2_xmin = -0.05;
#line 296 "PSI_Focus.instr"
  mccFERMITOF2_xmax = 0.05;
#line 296 "PSI_Focus.instr"
  mccFERMITOF2_ymin = -0.1;
#line 296 "PSI_Focus.instr"
  mccFERMITOF2_ymax = 0.1;
#line 202 "PSI_Focus.instr"
  mccFERMITOF2_zmin = 0;
#line 202 "PSI_Focus.instr"
  mccFERMITOF2_zmax = 0;
#line 295 "PSI_Focus.instr"
  mccFERMITOF2_bins = 30;
#line 203 "PSI_Focus.instr"
  mccFERMITOF2_min = -1e40;
#line 203 "PSI_Focus.instr"
  mccFERMITOF2_max = 1e40;
#line 203 "PSI_Focus.instr"
  mccFERMITOF2_restore_neutron = 0;
#line 203 "PSI_Focus.instr"
  mccFERMITOF2_radius = 0;
#line 296 "PSI_Focus.instr"
  if("auto, time") strncpy(mccFERMITOF2_options, "auto, time" ? "auto, time" : "", 16384); else mccFERMITOF2_options[0]='\0';
#line 295 "PSI_Focus.instr"
  if("FERMI_TOF2.dat.dat") strncpy(mccFERMITOF2_filename, "FERMI_TOF2.dat.dat" ? "FERMI_TOF2.dat.dat" : "", 16384); else mccFERMITOF2_filename[0]='\0';
#line 204 "PSI_Focus.instr"
  if("NULL") strncpy(mccFERMITOF2_geometry, "NULL" ? "NULL" : "", 16384); else mccFERMITOF2_geometry[0]='\0';
#line 205 "PSI_Focus.instr"
  if("NULL") strncpy(mccFERMITOF2_username1, "NULL" ? "NULL" : "", 16384); else mccFERMITOF2_username1[0]='\0';
#line 205 "PSI_Focus.instr"
  if("NULL") strncpy(mccFERMITOF2_username2, "NULL" ? "NULL" : "", 16384); else mccFERMITOF2_username2[0]='\0';
#line 205 "PSI_Focus.instr"
  if("NULL") strncpy(mccFERMITOF2_username3, "NULL" ? "NULL" : "", 16384); else mccFERMITOF2_username3[0]='\0';
#line 206 "PSI_Focus.instr"
  mccFERMITOF2_nowritefile = 0;
#line 15129 "PSI_Focus.c"

  SIG_MESSAGE("FERMITOF2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15136 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa2, mcrotaFERMITOF2);
  rot_transpose(mcrotaSAMPLE_SLIT, mctr1);
  rot_mul(mcrotaFERMITOF2, mctr1, mcrotrFERMITOF2);
  mctc1 = coords_set(
#line 297 "PSI_Focus.instr"
    0,
#line 297 "PSI_Focus.instr"
    0,
#line 297 "PSI_Focus.instr"
    1.2);
#line 15147 "PSI_Focus.c"
  rot_transpose(mcrotaa2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaFERMITOF2 = coords_add(mcposaa2, mctc2);
  mctc1 = coords_sub(mcposaSAMPLE_SLIT, mcposaFERMITOF2);
  mcposrFERMITOF2 = rot_apply(mcrotaFERMITOF2, mctc1);
  mcDEBUG_COMPONENT("FERMITOF2", mcposaFERMITOF2, mcrotaFERMITOF2)
  mccomp_posa[36] = mcposaFERMITOF2;
  mccomp_posr[36] = mcposrFERMITOF2;
  mcNCounter[36]  = mcPCounter[36] = mcP2Counter[36] = 0;
  mcAbsorbProp[36]= 0;
    /* Component PSD_SAMPLE. */
  /* Setting parameters for component PSD_SAMPLE. */
  SIG_MESSAGE("PSD_SAMPLE (Init:SetPar)");
#line 301 "PSI_Focus.instr"
  if("psdsample.dat") strncpy(mccPSD_SAMPLE_filename, "psdsample.dat" ? "psdsample.dat" : "", 16384); else mccPSD_SAMPLE_filename[0]='\0';
#line 50 "PSI_Focus.instr"
  mccPSD_SAMPLE_xmin = -0.05;
#line 50 "PSI_Focus.instr"
  mccPSD_SAMPLE_xmax = 0.05;
#line 50 "PSI_Focus.instr"
  mccPSD_SAMPLE_ymin = -0.05;
#line 50 "PSI_Focus.instr"
  mccPSD_SAMPLE_ymax = 0.05;
#line 300 "PSI_Focus.instr"
  mccPSD_SAMPLE_xwidth = 0.02;
#line 300 "PSI_Focus.instr"
  mccPSD_SAMPLE_yheight = 0.1;
#line 50 "PSI_Focus.instr"
  mccPSD_SAMPLE_restore_neutron = 0;
#line 50 "PSI_Focus.instr"
  mccPSD_SAMPLE_nowritefile = 0;
#line 15179 "PSI_Focus.c"

  SIG_MESSAGE("PSD_SAMPLE (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 302 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 302 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 302 "PSI_Focus.instr"
    (0)*DEG2RAD);
#line 15189 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa2, mcrotaPSD_SAMPLE);
  rot_transpose(mcrotaFERMITOF2, mctr1);
  rot_mul(mcrotaPSD_SAMPLE, mctr1, mcrotrPSD_SAMPLE);
  mctc1 = coords_set(
#line 302 "PSI_Focus.instr"
    0,
#line 302 "PSI_Focus.instr"
    0,
#line 302 "PSI_Focus.instr"
    1.45);
#line 15200 "PSI_Focus.c"
  rot_transpose(mcrotaa2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPSD_SAMPLE = coords_add(mcposaa2, mctc2);
  mctc1 = coords_sub(mcposaFERMITOF2, mcposaPSD_SAMPLE);
  mcposrPSD_SAMPLE = rot_apply(mcrotaPSD_SAMPLE, mctc1);
  mcDEBUG_COMPONENT("PSD_SAMPLE", mcposaPSD_SAMPLE, mcrotaPSD_SAMPLE)
  mccomp_posa[37] = mcposaPSD_SAMPLE;
  mccomp_posr[37] = mcposrPSD_SAMPLE;
  mcNCounter[37]  = mcPCounter[37] = mcP2Counter[37] = 0;
  mcAbsorbProp[37]= 0;
    /* Component DivMon_Sample. */
  /* Setting parameters for component DivMon_Sample. */
  SIG_MESSAGE("DivMon_Sample (Init:SetPar)");
#line 305 "PSI_Focus.instr"
  if("div2.dat") strncpy(mccDivMon_Sample_filename, "div2.dat" ? "div2.dat" : "", 16384); else mccDivMon_Sample_filename[0]='\0';
#line 55 "PSI_Focus.instr"
  mccDivMon_Sample_xmin = -0.05;
#line 55 "PSI_Focus.instr"
  mccDivMon_Sample_xmax = 0.05;
#line 55 "PSI_Focus.instr"
  mccDivMon_Sample_ymin = -0.05;
#line 55 "PSI_Focus.instr"
  mccDivMon_Sample_ymax = 0.05;
#line 306 "PSI_Focus.instr"
  mccDivMon_Sample_xwidth = 0.01;
#line 306 "PSI_Focus.instr"
  mccDivMon_Sample_yheight = 0.06;
#line 307 "PSI_Focus.instr"
  mccDivMon_Sample_maxdiv_h = 3;
#line 307 "PSI_Focus.instr"
  mccDivMon_Sample_maxdiv_v = 3;
#line 56 "PSI_Focus.instr"
  mccDivMon_Sample_restore_neutron = 0;
#line 56 "PSI_Focus.instr"
  mccDivMon_Sample_nx = 0;
#line 56 "PSI_Focus.instr"
  mccDivMon_Sample_ny = 0;
#line 56 "PSI_Focus.instr"
  mccDivMon_Sample_nz = 1;
#line 56 "PSI_Focus.instr"
  mccDivMon_Sample_nowritefile = 0;
#line 15242 "PSI_Focus.c"

  SIG_MESSAGE("DivMon_Sample (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 308 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 308 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 308 "PSI_Focus.instr"
    (0)*DEG2RAD);
#line 15252 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa2, mcrotaDivMon_Sample);
  rot_transpose(mcrotaPSD_SAMPLE, mctr1);
  rot_mul(mcrotaDivMon_Sample, mctr1, mcrotrDivMon_Sample);
  mctc1 = coords_set(
#line 308 "PSI_Focus.instr"
    0,
#line 308 "PSI_Focus.instr"
    0,
#line 308 "PSI_Focus.instr"
    1.451);
#line 15263 "PSI_Focus.c"
  rot_transpose(mcrotaa2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaDivMon_Sample = coords_add(mcposaa2, mctc2);
  mctc1 = coords_sub(mcposaPSD_SAMPLE, mcposaDivMon_Sample);
  mcposrDivMon_Sample = rot_apply(mcrotaDivMon_Sample, mctc1);
  mcDEBUG_COMPONENT("DivMon_Sample", mcposaDivMon_Sample, mcrotaDivMon_Sample)
  mccomp_posa[38] = mcposaDivMon_Sample;
  mccomp_posr[38] = mcposrDivMon_Sample;
  mcNCounter[38]  = mcPCounter[38] = mcP2Counter[38] = 0;
  mcAbsorbProp[38]= 0;
    /* Component EMON_SAMPLE. */
  /* Setting parameters for component EMON_SAMPLE. */
  SIG_MESSAGE("EMON_SAMPLE (Init:SetPar)");
#line 313 "PSI_Focus.instr"
  if("emon_sample.dat") strncpy(mccEMON_SAMPLE_filename, "emon_sample.dat" ? "emon_sample.dat" : "", 16384); else mccEMON_SAMPLE_filename[0]='\0';
#line 53 "PSI_Focus.instr"
  mccEMON_SAMPLE_xmin = -0.05;
#line 53 "PSI_Focus.instr"
  mccEMON_SAMPLE_xmax = 0.05;
#line 53 "PSI_Focus.instr"
  mccEMON_SAMPLE_ymin = -0.05;
#line 53 "PSI_Focus.instr"
  mccEMON_SAMPLE_ymax = 0.05;
#line 311 "PSI_Focus.instr"
  mccEMON_SAMPLE_xwidth = 0.01;
#line 311 "PSI_Focus.instr"
  mccEMON_SAMPLE_yheight = 0.06;
#line 312 "PSI_Focus.instr"
  mccEMON_SAMPLE_Emin = EMIN;
#line 312 "PSI_Focus.instr"
  mccEMON_SAMPLE_Emax = EMAX;
#line 54 "PSI_Focus.instr"
  mccEMON_SAMPLE_restore_neutron = 0;
#line 54 "PSI_Focus.instr"
  mccEMON_SAMPLE_nowritefile = 0;
#line 15299 "PSI_Focus.c"

  SIG_MESSAGE("EMON_SAMPLE (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 314 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 314 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 314 "PSI_Focus.instr"
    (0)*DEG2RAD);
#line 15309 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa2, mcrotaEMON_SAMPLE);
  rot_transpose(mcrotaDivMon_Sample, mctr1);
  rot_mul(mcrotaEMON_SAMPLE, mctr1, mcrotrEMON_SAMPLE);
  mctc1 = coords_set(
#line 314 "PSI_Focus.instr"
    0,
#line 314 "PSI_Focus.instr"
    0,
#line 314 "PSI_Focus.instr"
    1.452);
#line 15320 "PSI_Focus.c"
  rot_transpose(mcrotaa2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaEMON_SAMPLE = coords_add(mcposaa2, mctc2);
  mctc1 = coords_sub(mcposaDivMon_Sample, mcposaEMON_SAMPLE);
  mcposrEMON_SAMPLE = rot_apply(mcrotaEMON_SAMPLE, mctc1);
  mcDEBUG_COMPONENT("EMON_SAMPLE", mcposaEMON_SAMPLE, mcrotaEMON_SAMPLE)
  mccomp_posa[39] = mcposaEMON_SAMPLE;
  mccomp_posr[39] = mcposrEMON_SAMPLE;
  mcNCounter[39]  = mcPCounter[39] = mcP2Counter[39] = 0;
  mcAbsorbProp[39]= 0;
    /* Component a3. */
  /* Setting parameters for component a3. */
  SIG_MESSAGE("a3 (Init:SetPar)");

  SIG_MESSAGE("a3 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15340 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa2, mcrotaa3);
  rot_transpose(mcrotaEMON_SAMPLE, mctr1);
  rot_mul(mcrotaa3, mctr1, mcrotra3);
  mctc1 = coords_set(
#line 319 "PSI_Focus.instr"
    0,
#line 319 "PSI_Focus.instr"
    0,
#line 319 "PSI_Focus.instr"
    1.5);
#line 15351 "PSI_Focus.c"
  rot_transpose(mcrotaa2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaa3 = coords_add(mcposaa2, mctc2);
  mctc1 = coords_sub(mcposaEMON_SAMPLE, mcposaa3);
  mcposra3 = rot_apply(mcrotaa3, mctc1);
  mcDEBUG_COMPONENT("a3", mcposaa3, mcrotaa3)
  mccomp_posa[40] = mcposaa3;
  mccomp_posr[40] = mcposra3;
  mcNCounter[40]  = mcPCounter[40] = mcP2Counter[40] = 0;
  mcAbsorbProp[40]= 0;
    /* Component Sample. */
  /* Setting parameters for component Sample. */
  SIG_MESSAGE("Sample (Init:SetPar)");
#line 322 "PSI_Focus.instr"
  mccSample_radius = 0.008;
#line 91 "PSI_Focus.instr"
  mccSample_thickness = 0;
#line 91 "PSI_Focus.instr"
  mccSample_zdepth = 0;
#line 91 "PSI_Focus.instr"
  mccSample_Vc = 13.827;
#line 91 "PSI_Focus.instr"
  mccSample_sigma_abs = 5.08;
#line 91 "PSI_Focus.instr"
  mccSample_sigma_inc = 5.08;
#line 92 "PSI_Focus.instr"
  mccSample_radius_i = 0;
#line 92 "PSI_Focus.instr"
  mccSample_radius_o = 0;
#line 92 "PSI_Focus.instr"
  mccSample_h = 0;
#line 92 "PSI_Focus.instr"
  mccSample_focus_r = 0;
#line 322 "PSI_Focus.instr"
  mccSample_pack = 1;
#line 92 "PSI_Focus.instr"
  mccSample_frac = 1;
#line 92 "PSI_Focus.instr"
  mccSample_f_QE = 0;
#line 92 "PSI_Focus.instr"
  mccSample_gamma = 0;
#line 93 "PSI_Focus.instr"
  mccSample_target_x = 0;
#line 93 "PSI_Focus.instr"
  mccSample_target_y = 0;
#line 93 "PSI_Focus.instr"
  mccSample_target_z = 0;
#line 322 "PSI_Focus.instr"
  mccSample_focus_xw = 0.336;
#line 322 "PSI_Focus.instr"
  mccSample_focus_yh = 0.4;
#line 94 "PSI_Focus.instr"
  mccSample_focus_aw = 0;
#line 94 "PSI_Focus.instr"
  mccSample_focus_ah = 0;
#line 94 "PSI_Focus.instr"
  mccSample_xwidth = 0;
#line 322 "PSI_Focus.instr"
  mccSample_yheight = 0.055;
#line 94 "PSI_Focus.instr"
  mccSample_zthick = 0;
#line 94 "PSI_Focus.instr"
  mccSample_rad_sphere = 0;
#line 94 "PSI_Focus.instr"
  mccSample_sig_a = 0;
#line 94 "PSI_Focus.instr"
  mccSample_sig_i = 0;
#line 94 "PSI_Focus.instr"
  mccSample_V0 = 0;
#line 323 "PSI_Focus.instr"
  mccSample_target_index = + 3;
#line 94 "PSI_Focus.instr"
  mccSample_multiples = 1;
#line 15425 "PSI_Focus.c"

  SIG_MESSAGE("Sample (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15432 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa3, mcrotaSample);
  rot_transpose(mcrotaa3, mctr1);
  rot_mul(mcrotaSample, mctr1, mcrotrSample);
  mctc1 = coords_set(
#line 324 "PSI_Focus.instr"
    0,
#line 324 "PSI_Focus.instr"
    0,
#line 324 "PSI_Focus.instr"
    0);
#line 15443 "PSI_Focus.c"
  rot_transpose(mcrotaa3, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaSample = coords_add(mcposaa3, mctc2);
  mctc1 = coords_sub(mcposaa3, mcposaSample);
  mcposrSample = rot_apply(mcrotaSample, mctc1);
  mcDEBUG_COMPONENT("Sample", mcposaSample, mcrotaSample)
  mccomp_posa[41] = mcposaSample;
  mccomp_posr[41] = mcposrSample;
  mcNCounter[41]  = mcPCounter[41] = mcP2Counter[41] = 0;
  mcAbsorbProp[41]= 0;
    /* Component TOF_Det. */
  /* Setting parameters for component TOF_Det. */
  SIG_MESSAGE("TOF_Det (Init:SetPar)");
#line 327 "PSI_Focus.instr"
  mccTOF_Det_xwidth = 1.5;
#line 327 "PSI_Focus.instr"
  mccTOF_Det_yheight = 0.2;
#line 201 "PSI_Focus.instr"
  mccTOF_Det_zdepth = 0;
#line 202 "PSI_Focus.instr"
  mccTOF_Det_xmin = 0;
#line 202 "PSI_Focus.instr"
  mccTOF_Det_xmax = 0;
#line 202 "PSI_Focus.instr"
  mccTOF_Det_ymin = 0;
#line 202 "PSI_Focus.instr"
  mccTOF_Det_ymax = 0;
#line 202 "PSI_Focus.instr"
  mccTOF_Det_zmin = 0;
#line 202 "PSI_Focus.instr"
  mccTOF_Det_zmax = 0;
#line 329 "PSI_Focus.instr"
  mccTOF_Det_bins = 100;
#line 203 "PSI_Focus.instr"
  mccTOF_Det_min = -1e40;
#line 203 "PSI_Focus.instr"
  mccTOF_Det_max = 1e40;
#line 203 "PSI_Focus.instr"
  mccTOF_Det_restore_neutron = 0;
#line 203 "PSI_Focus.instr"
  mccTOF_Det_radius = 0;
#line 328 "PSI_Focus.instr"
  if("auto time, abs angle limits=[10 180], banana, parallel") strncpy(mccTOF_Det_options, "auto time, abs angle limits=[10 180], banana, parallel" ? "auto time, abs angle limits=[10 180], banana, parallel" : "", 16384); else mccTOF_Det_options[0]='\0';
#line 204 "PSI_Focus.instr"
  if("NULL") strncpy(mccTOF_Det_filename, "NULL" ? "NULL" : "", 16384); else mccTOF_Det_filename[0]='\0';
#line 204 "PSI_Focus.instr"
  if("NULL") strncpy(mccTOF_Det_geometry, "NULL" ? "NULL" : "", 16384); else mccTOF_Det_geometry[0]='\0';
#line 205 "PSI_Focus.instr"
  if("NULL") strncpy(mccTOF_Det_username1, "NULL" ? "NULL" : "", 16384); else mccTOF_Det_username1[0]='\0';
#line 205 "PSI_Focus.instr"
  if("NULL") strncpy(mccTOF_Det_username2, "NULL" ? "NULL" : "", 16384); else mccTOF_Det_username2[0]='\0';
#line 205 "PSI_Focus.instr"
  if("NULL") strncpy(mccTOF_Det_username3, "NULL" ? "NULL" : "", 16384); else mccTOF_Det_username3[0]='\0';
#line 206 "PSI_Focus.instr"
  mccTOF_Det_nowritefile = 0;
#line 15499 "PSI_Focus.c"

  SIG_MESSAGE("TOF_Det (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15506 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaSample, mcrotaTOF_Det);
  rot_transpose(mcrotaSample, mctr1);
  rot_mul(mcrotaTOF_Det, mctr1, mcrotrTOF_Det);
  mctc1 = coords_set(
#line 330 "PSI_Focus.instr"
    0,
#line 330 "PSI_Focus.instr"
    0,
#line 330 "PSI_Focus.instr"
    0);
#line 15517 "PSI_Focus.c"
  rot_transpose(mcrotaSample, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaTOF_Det = coords_add(mcposaSample, mctc2);
  mctc1 = coords_sub(mcposaSample, mcposaTOF_Det);
  mcposrTOF_Det = rot_apply(mcrotaTOF_Det, mctc1);
  mcDEBUG_COMPONENT("TOF_Det", mcposaTOF_Det, mcrotaTOF_Det)
  mccomp_posa[42] = mcposaTOF_Det;
  mccomp_posr[42] = mcposrTOF_Det;
  mcNCounter[42]  = mcPCounter[42] = mcP2Counter[42] = 0;
  mcAbsorbProp[42]= 0;
    /* Component FoDet. */
  /* Setting parameters for component FoDet. */
  SIG_MESSAGE("FoDet (Init:SetPar)");
#line 333 "PSI_Focus.instr"
  mccFoDet_xwidth = 0.336;
#line 333 "PSI_Focus.instr"
  mccFoDet_yheight = 0.4;
#line 201 "PSI_Focus.instr"
  mccFoDet_zdepth = 0;
#line 202 "PSI_Focus.instr"
  mccFoDet_xmin = 0;
#line 202 "PSI_Focus.instr"
  mccFoDet_xmax = 0;
#line 202 "PSI_Focus.instr"
  mccFoDet_ymin = 0;
#line 202 "PSI_Focus.instr"
  mccFoDet_ymax = 0;
#line 202 "PSI_Focus.instr"
  mccFoDet_zmin = 0;
#line 202 "PSI_Focus.instr"
  mccFoDet_zmax = 0;
#line 203 "PSI_Focus.instr"
  mccFoDet_bins = 0;
#line 203 "PSI_Focus.instr"
  mccFoDet_min = -1e40;
#line 203 "PSI_Focus.instr"
  mccFoDet_max = 1e40;
#line 203 "PSI_Focus.instr"
  mccFoDet_restore_neutron = 0;
#line 203 "PSI_Focus.instr"
  mccFoDet_radius = 0;
#line 334 "PSI_Focus.instr"
  if("t auto file=TofDet.dat") strncpy(mccFoDet_options, "t auto file=TofDet.dat" ? "t auto file=TofDet.dat" : "", 16384); else mccFoDet_options[0]='\0';
#line 204 "PSI_Focus.instr"
  if("NULL") strncpy(mccFoDet_filename, "NULL" ? "NULL" : "", 16384); else mccFoDet_filename[0]='\0';
#line 204 "PSI_Focus.instr"
  if("NULL") strncpy(mccFoDet_geometry, "NULL" ? "NULL" : "", 16384); else mccFoDet_geometry[0]='\0';
#line 205 "PSI_Focus.instr"
  if("NULL") strncpy(mccFoDet_username1, "NULL" ? "NULL" : "", 16384); else mccFoDet_username1[0]='\0';
#line 205 "PSI_Focus.instr"
  if("NULL") strncpy(mccFoDet_username2, "NULL" ? "NULL" : "", 16384); else mccFoDet_username2[0]='\0';
#line 205 "PSI_Focus.instr"
  if("NULL") strncpy(mccFoDet_username3, "NULL" ? "NULL" : "", 16384); else mccFoDet_username3[0]='\0';
#line 206 "PSI_Focus.instr"
  mccFoDet_nowritefile = 0;
#line 15573 "PSI_Focus.c"

  SIG_MESSAGE("FoDet (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 335 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 335 "PSI_Focus.instr"
    (mcipDET)*DEG2RAD,
#line 335 "PSI_Focus.instr"
    (0)*DEG2RAD);
#line 15583 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa3, mcrotaFoDet);
  rot_transpose(mcrotaTOF_Det, mctr1);
  rot_mul(mcrotaFoDet, mctr1, mcrotrFoDet);
  mctc1 = coords_set(
#line 335 "PSI_Focus.instr"
    0,
#line 335 "PSI_Focus.instr"
    0,
#line 335 "PSI_Focus.instr"
    2.5);
#line 15594 "PSI_Focus.c"
  rot_transpose(mcrotaa3, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaFoDet = coords_add(mcposaa3, mctc2);
  mctc1 = coords_sub(mcposaTOF_Det, mcposaFoDet);
  mcposrFoDet = rot_apply(mcrotaFoDet, mctc1);
  mcDEBUG_COMPONENT("FoDet", mcposaFoDet, mcrotaFoDet)
  mccomp_posa[43] = mcposaFoDet;
  mccomp_posr[43] = mcposrFoDet;
  mcNCounter[43]  = mcPCounter[43] = mcP2Counter[43] = 0;
  mcAbsorbProp[43]= 0;
    /* Component EMON_DET. */
  /* Setting parameters for component EMON_DET. */
  SIG_MESSAGE("EMON_DET (Init:SetPar)");
#line 339 "PSI_Focus.instr"
  if("emon_det.dat") strncpy(mccEMON_DET_filename, "emon_det.dat" ? "emon_det.dat" : "", 16384); else mccEMON_DET_filename[0]='\0';
#line 53 "PSI_Focus.instr"
  mccEMON_DET_xmin = -0.05;
#line 53 "PSI_Focus.instr"
  mccEMON_DET_xmax = 0.05;
#line 53 "PSI_Focus.instr"
  mccEMON_DET_ymin = -0.05;
#line 53 "PSI_Focus.instr"
  mccEMON_DET_ymax = 0.05;
#line 338 "PSI_Focus.instr"
  mccEMON_DET_xwidth = 0.4;
#line 338 "PSI_Focus.instr"
  mccEMON_DET_yheight = 0.4;
#line 339 "PSI_Focus.instr"
  mccEMON_DET_Emin = EMIN;
#line 339 "PSI_Focus.instr"
  mccEMON_DET_Emax = EMAX;
#line 54 "PSI_Focus.instr"
  mccEMON_DET_restore_neutron = 0;
#line 54 "PSI_Focus.instr"
  mccEMON_DET_nowritefile = 0;
#line 15630 "PSI_Focus.c"

  SIG_MESSAGE("EMON_DET (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 340 "PSI_Focus.instr"
    (0)*DEG2RAD,
#line 340 "PSI_Focus.instr"
    (mcipDET)*DEG2RAD,
#line 340 "PSI_Focus.instr"
    (0)*DEG2RAD);
#line 15640 "PSI_Focus.c"
  rot_mul(mctr1, mcrotaa3, mcrotaEMON_DET);
  rot_transpose(mcrotaFoDet, mctr1);
  rot_mul(mcrotaEMON_DET, mctr1, mcrotrEMON_DET);
  mctc1 = coords_set(
#line 340 "PSI_Focus.instr"
    0,
#line 340 "PSI_Focus.instr"
    0,
#line 340 "PSI_Focus.instr"
    2.501);
#line 15651 "PSI_Focus.c"
  rot_transpose(mcrotaa3, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaEMON_DET = coords_add(mcposaa3, mctc2);
  mctc1 = coords_sub(mcposaFoDet, mcposaEMON_DET);
  mcposrEMON_DET = rot_apply(mcrotaEMON_DET, mctc1);
  mcDEBUG_COMPONENT("EMON_DET", mcposaEMON_DET, mcrotaEMON_DET)
  mccomp_posa[44] = mcposaEMON_DET;
  mccomp_posr[44] = mcposrEMON_DET;
  mcNCounter[44]  = mcPCounter[44] = mcP2Counter[44] = 0;
  mcAbsorbProp[44]= 0;
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
#line 15688 "PSI_Focus.c"
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

  /* Initializations for component csource. */
  SIG_MESSAGE("csource (Init)");
#define mccompcurname  csource
#define mccompcurtype  Source_gen
#define mccompcurindex 2
#define p_in mcccsource_p_in
#define lambda1 mcccsource_lambda1
#define lambda2 mcccsource_lambda2
#define lambda3 mcccsource_lambda3
#define pTable mcccsource_pTable
#define pTable_x mcccsource_pTable_x
#define pTable_y mcccsource_pTable_y
#define pTable_xmin mcccsource_pTable_xmin
#define pTable_xmax mcccsource_pTable_xmax
#define pTable_xsum mcccsource_pTable_xsum
#define pTable_ymin mcccsource_pTable_ymin
#define pTable_ymax mcccsource_pTable_ymax
#define pTable_ysum mcccsource_pTable_ysum
#define pTable_dxmin mcccsource_pTable_dxmin
#define pTable_dxmax mcccsource_pTable_dxmax
#define pTable_dymin mcccsource_pTable_dymin
#define pTable_dymax mcccsource_pTable_dymax
#define flux_file mcccsource_flux_file
#define xdiv_file mcccsource_xdiv_file
#define ydiv_file mcccsource_ydiv_file
#define radius mcccsource_radius
#define dist mcccsource_dist
#define focus_xw mcccsource_focus_xw
#define focus_yh mcccsource_focus_yh
#define focus_aw mcccsource_focus_aw
#define focus_ah mcccsource_focus_ah
#define E0 mcccsource_E0
#define dE mcccsource_dE
#define lambda0 mcccsource_lambda0
#define dlambda mcccsource_dlambda
#define I1 mcccsource_I1
#define yheight mcccsource_yheight
#define xwidth mcccsource_xwidth
#define verbose mcccsource_verbose
#define T1 mcccsource_T1
#define flux_file_perAA mcccsource_flux_file_perAA
#define flux_file_log mcccsource_flux_file_log
#define Lmin mcccsource_Lmin
#define Lmax mcccsource_Lmax
#define Emin mcccsource_Emin
#define Emax mcccsource_Emax
#define T2 mcccsource_T2
#define I2 mcccsource_I2
#define T3 mcccsource_T3
#define I3 mcccsource_I3
#define zdepth mcccsource_zdepth
#define target_index mcccsource_target_index
#line 206 "/usr/share/mcstas/2.5/sources/Source_gen.comp"
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
#line 16025 "PSI_Focus.c"
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

  /* Initializations for component guide1. */
  SIG_MESSAGE("guide1 (Init)");
#define mccompcurname  guide1
#define mccompcurtype  Guide
#define mccompcurindex 3
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
#line 74 "/usr/share/mcstas/2.5/optics/Guide.comp"
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
#line 16113 "PSI_Focus.c"
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

  /* Initializations for component guide2. */
  SIG_MESSAGE("guide2 (Init)");
#define mccompcurname  guide2
#define mccompcurtype  Bender
#define mccompcurindex 4
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
#line 112 "/usr/share/mcstas/2.5/optics/Bender.comp"
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
#line 16179 "PSI_Focus.c"
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

  /* Initializations for component bunker. */
  SIG_MESSAGE("bunker (Init)");
#define mccompcurname  bunker
#define mccompcurtype  Guide
#define mccompcurindex 5
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
#line 74 "/usr/share/mcstas/2.5/optics/Guide.comp"
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
#line 16244 "PSI_Focus.c"
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
#define mccompcurindex 6
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
#line 74 "/usr/share/mcstas/2.5/optics/Guide.comp"
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
#line 16297 "PSI_Focus.c"
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

  /* Initializations for component lambdaGuideExit. */
  SIG_MESSAGE("lambdaGuideExit (Init)");
#define mccompcurname  lambdaGuideExit
#define mccompcurtype  L_monitor
#define mccompcurindex 7
#define nL mcclambdaGuideExit_nL
#define L_N mcclambdaGuideExit_L_N
#define L_p mcclambdaGuideExit_L_p
#define L_p2 mcclambdaGuideExit_L_p2
#define filename mcclambdaGuideExit_filename
#define xmin mcclambdaGuideExit_xmin
#define xmax mcclambdaGuideExit_xmax
#define ymin mcclambdaGuideExit_ymin
#define ymax mcclambdaGuideExit_ymax
#define xwidth mcclambdaGuideExit_xwidth
#define yheight mcclambdaGuideExit_yheight
#define Lmin mcclambdaGuideExit_Lmin
#define Lmax mcclambdaGuideExit_Lmax
#define restore_neutron mcclambdaGuideExit_restore_neutron
#define nowritefile mcclambdaGuideExit_nowritefile
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
#line 16355 "PSI_Focus.c"
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

  /* Initializations for component DivMonGuideExit. */
  SIG_MESSAGE("DivMonGuideExit (Init)");
#define mccompcurname  DivMonGuideExit
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 8
#define nh mccDivMonGuideExit_nh
#define nv mccDivMonGuideExit_nv
#define Div_N mccDivMonGuideExit_Div_N
#define Div_p mccDivMonGuideExit_Div_p
#define Div_p2 mccDivMonGuideExit_Div_p2
#define filename mccDivMonGuideExit_filename
#define xmin mccDivMonGuideExit_xmin
#define xmax mccDivMonGuideExit_xmax
#define ymin mccDivMonGuideExit_ymin
#define ymax mccDivMonGuideExit_ymax
#define xwidth mccDivMonGuideExit_xwidth
#define yheight mccDivMonGuideExit_yheight
#define maxdiv_h mccDivMonGuideExit_maxdiv_h
#define maxdiv_v mccDivMonGuideExit_maxdiv_v
#define restore_neutron mccDivMonGuideExit_restore_neutron
#define nx mccDivMonGuideExit_nx
#define ny mccDivMonGuideExit_ny
#define nz mccDivMonGuideExit_nz
#define nowritefile mccDivMonGuideExit_nowritefile
#line 66 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
{
int i,j;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("Divergence_monitor: %s: Null detection area !\n"
                   "ERROR               (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    for (i=0; i<nh; i++)
     for (j=0; j<nv; j++)
     {
      Div_N[i][j] = 0;
      Div_p[i][j] = 0;
      Div_p2[i][j] = 0;
     }
    NORM(nx,ny,nz);
}
#line 16422 "PSI_Focus.c"
#undef nowritefile
#undef nz
#undef ny
#undef nx
#undef restore_neutron
#undef maxdiv_v
#undef maxdiv_h
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
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component PSDGuideExit. */
  SIG_MESSAGE("PSDGuideExit (Init)");
#define mccompcurname  PSDGuideExit
#define mccompcurtype  PSD_monitor
#define mccompcurindex 9
#define nx mccPSDGuideExit_nx
#define ny mccPSDGuideExit_ny
#define PSD_N mccPSDGuideExit_PSD_N
#define PSD_p mccPSDGuideExit_PSD_p
#define PSD_p2 mccPSDGuideExit_PSD_p2
#define filename mccPSDGuideExit_filename
#define xmin mccPSDGuideExit_xmin
#define xmax mccPSDGuideExit_xmax
#define ymin mccPSDGuideExit_ymin
#define ymax mccPSDGuideExit_ymax
#define xwidth mccPSDGuideExit_xwidth
#define yheight mccPSDGuideExit_yheight
#define restore_neutron mccPSDGuideExit_restore_neutron
#define nowritefile mccPSDGuideExit_nowritefile
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
#line 16487 "PSI_Focus.c"
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

  /* Initializations for component FOCUSguide. */
  SIG_MESSAGE("FOCUSguide (Init)");
#define mccompcurname  FOCUSguide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 10
#define w1c mccFOCUSguide_w1c
#define w2c mccFOCUSguide_w2c
#define ww mccFOCUSguide_ww
#define hh mccFOCUSguide_hh
#define whalf mccFOCUSguide_whalf
#define hhalf mccFOCUSguide_hhalf
#define lwhalf mccFOCUSguide_lwhalf
#define lhhalf mccFOCUSguide_lhhalf
#define w1 mccFOCUSguide_w1
#define h1 mccFOCUSguide_h1
#define w2 mccFOCUSguide_w2
#define h2 mccFOCUSguide_h2
#define l mccFOCUSguide_l
#define R0 mccFOCUSguide_R0
#define Qc mccFOCUSguide_Qc
#define alpha mccFOCUSguide_alpha
#define m mccFOCUSguide_m
#define nslit mccFOCUSguide_nslit
#define d mccFOCUSguide_d
#define Qcx mccFOCUSguide_Qcx
#define Qcy mccFOCUSguide_Qcy
#define alphax mccFOCUSguide_alphax
#define alphay mccFOCUSguide_alphay
#define W mccFOCUSguide_W
#define mx mccFOCUSguide_mx
#define my mccFOCUSguide_my
#define nu mccFOCUSguide_nu
#define phase mccFOCUSguide_phase
#line 89 "/usr/share/mcstas/2.5/optics/Guide_channeled.comp"
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
#line 16580 "PSI_Focus.c"
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

  /* Initializations for component FirstChopper. */
  SIG_MESSAGE("FirstChopper (Init)");
#define mccompcurname  FirstChopper
#define mccompcurtype  DiskChopper
#define mccompcurindex 11
#define Tg mccFirstChopper_Tg
#define To mccFirstChopper_To
#define delta_y mccFirstChopper_delta_y
#define height mccFirstChopper_height
#define omega mccFirstChopper_omega
#define theta_0 mccFirstChopper_theta_0
#define radius mccFirstChopper_radius
#define yheight mccFirstChopper_yheight
#define nu mccFirstChopper_nu
#define nslit mccFirstChopper_nslit
#define jitter mccFirstChopper_jitter
#define delay mccFirstChopper_delay
#define isfirst mccFirstChopper_isfirst
#define n_pulse mccFirstChopper_n_pulse
#define abs_out mccFirstChopper_abs_out
#define phase mccFirstChopper_phase
#define xwidth mccFirstChopper_xwidth
#define verbose mccFirstChopper_verbose
#line 70 "/usr/share/mcstas/2.5/optics/DiskChopper.comp"
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
#line 16698 "PSI_Focus.c"
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

  /* Initializations for component DISCTOF. */
  SIG_MESSAGE("DISCTOF (Init)");
#define mccompcurname  DISCTOF
#define mccompcurtype  Monitor_nD
#define mccompcurindex 12
#define user1 mccDISCTOF_user1
#define user2 mccDISCTOF_user2
#define user3 mccDISCTOF_user3
#define DEFS mccDISCTOF_DEFS
#define Vars mccDISCTOF_Vars
#define detector mccDISCTOF_detector
#define offdata mccDISCTOF_offdata
#define xwidth mccDISCTOF_xwidth
#define yheight mccDISCTOF_yheight
#define zdepth mccDISCTOF_zdepth
#define xmin mccDISCTOF_xmin
#define xmax mccDISCTOF_xmax
#define ymin mccDISCTOF_ymin
#define ymax mccDISCTOF_ymax
#define zmin mccDISCTOF_zmin
#define zmax mccDISCTOF_zmax
#define bins mccDISCTOF_bins
#define min mccDISCTOF_min
#define max mccDISCTOF_max
#define restore_neutron mccDISCTOF_restore_neutron
#define radius mccDISCTOF_radius
#define options mccDISCTOF_options
#define filename mccDISCTOF_filename
#define geometry mccDISCTOF_geometry
#define username1 mccDISCTOF_username1
#define username2 mccDISCTOF_username2
#define username3 mccDISCTOF_username3
#define nowritefile mccDISCTOF_nowritefile
#line 229 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
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
#line 16833 "PSI_Focus.c"
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

  /* Initializations for component PSDmon1Chopper. */
  SIG_MESSAGE("PSDmon1Chopper (Init)");
#define mccompcurname  PSDmon1Chopper
#define mccompcurtype  PSD_monitor
#define mccompcurindex 13
#define nx mccPSDmon1Chopper_nx
#define ny mccPSDmon1Chopper_ny
#define PSD_N mccPSDmon1Chopper_PSD_N
#define PSD_p mccPSDmon1Chopper_PSD_p
#define PSD_p2 mccPSDmon1Chopper_PSD_p2
#define filename mccPSDmon1Chopper_filename
#define xmin mccPSDmon1Chopper_xmin
#define xmax mccPSDmon1Chopper_xmax
#define ymin mccPSDmon1Chopper_ymin
#define ymax mccPSDmon1Chopper_ymax
#define xwidth mccPSDmon1Chopper_xwidth
#define yheight mccPSDmon1Chopper_yheight
#define restore_neutron mccPSDmon1Chopper_restore_neutron
#define nowritefile mccPSDmon1Chopper_nowritefile
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
#line 16907 "PSI_Focus.c"
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

  /* Initializations for component VacuumTube1entry. */
  SIG_MESSAGE("VacuumTube1entry (Init)");
#define mccompcurname  VacuumTube1entry
#define mccompcurtype  Slit
#define mccompcurindex 14
#define xmin mccVacuumTube1entry_xmin
#define xmax mccVacuumTube1entry_xmax
#define ymin mccVacuumTube1entry_ymin
#define ymax mccVacuumTube1entry_ymax
#define radius mccVacuumTube1entry_radius
#define xwidth mccVacuumTube1entry_xwidth
#define yheight mccVacuumTube1entry_yheight
#line 50 "/usr/share/mcstas/2.5/optics/Slit.comp"
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
#line 16958 "PSI_Focus.c"
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

  /* Initializations for component VacuumTube1exit. */
  SIG_MESSAGE("VacuumTube1exit (Init)");
#define mccompcurname  VacuumTube1exit
#define mccompcurtype  Slit
#define mccompcurindex 15
#define xmin mccVacuumTube1exit_xmin
#define xmax mccVacuumTube1exit_xmax
#define ymin mccVacuumTube1exit_ymin
#define ymax mccVacuumTube1exit_ymax
#define radius mccVacuumTube1exit_radius
#define xwidth mccVacuumTube1exit_xwidth
#define yheight mccVacuumTube1exit_yheight
#line 50 "/usr/share/mcstas/2.5/optics/Slit.comp"
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
#line 17002 "PSI_Focus.c"
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

  /* Initializations for component VacuumTube2entry. */
  SIG_MESSAGE("VacuumTube2entry (Init)");
#define mccompcurname  VacuumTube2entry
#define mccompcurtype  Slit
#define mccompcurindex 16
#define xmin mccVacuumTube2entry_xmin
#define xmax mccVacuumTube2entry_xmax
#define ymin mccVacuumTube2entry_ymin
#define ymax mccVacuumTube2entry_ymax
#define radius mccVacuumTube2entry_radius
#define xwidth mccVacuumTube2entry_xwidth
#define yheight mccVacuumTube2entry_yheight
#line 50 "/usr/share/mcstas/2.5/optics/Slit.comp"
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
#line 17046 "PSI_Focus.c"
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

  /* Initializations for component VacuumTube2exit. */
  SIG_MESSAGE("VacuumTube2exit (Init)");
#define mccompcurname  VacuumTube2exit
#define mccompcurtype  Slit
#define mccompcurindex 17
#define xmin mccVacuumTube2exit_xmin
#define xmax mccVacuumTube2exit_xmax
#define ymin mccVacuumTube2exit_ymin
#define ymax mccVacuumTube2exit_ymax
#define radius mccVacuumTube2exit_radius
#define xwidth mccVacuumTube2exit_xwidth
#define yheight mccVacuumTube2exit_yheight
#line 50 "/usr/share/mcstas/2.5/optics/Slit.comp"
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
#line 17090 "PSI_Focus.c"
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

  /* Initializations for component VacuumTube3entry. */
  SIG_MESSAGE("VacuumTube3entry (Init)");
#define mccompcurname  VacuumTube3entry
#define mccompcurtype  Slit
#define mccompcurindex 18
#define xmin mccVacuumTube3entry_xmin
#define xmax mccVacuumTube3entry_xmax
#define ymin mccVacuumTube3entry_ymin
#define ymax mccVacuumTube3entry_ymax
#define radius mccVacuumTube3entry_radius
#define xwidth mccVacuumTube3entry_xwidth
#define yheight mccVacuumTube3entry_yheight
#line 50 "/usr/share/mcstas/2.5/optics/Slit.comp"
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
#line 17134 "PSI_Focus.c"
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

  /* Initializations for component VacuumTube3exit. */
  SIG_MESSAGE("VacuumTube3exit (Init)");
#define mccompcurname  VacuumTube3exit
#define mccompcurtype  Slit
#define mccompcurindex 19
#define xmin mccVacuumTube3exit_xmin
#define xmax mccVacuumTube3exit_xmax
#define ymin mccVacuumTube3exit_ymin
#define ymax mccVacuumTube3exit_ymax
#define radius mccVacuumTube3exit_radius
#define xwidth mccVacuumTube3exit_xwidth
#define yheight mccVacuumTube3exit_yheight
#line 50 "/usr/share/mcstas/2.5/optics/Slit.comp"
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
#line 17178 "PSI_Focus.c"
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

  /* Initializations for component PSDmonMono. */
  SIG_MESSAGE("PSDmonMono (Init)");
#define mccompcurname  PSDmonMono
#define mccompcurtype  PSD_monitor
#define mccompcurindex 20
#define nx mccPSDmonMono_nx
#define ny mccPSDmonMono_ny
#define PSD_N mccPSDmonMono_PSD_N
#define PSD_p mccPSDmonMono_PSD_p
#define PSD_p2 mccPSDmonMono_PSD_p2
#define filename mccPSDmonMono_filename
#define xmin mccPSDmonMono_xmin
#define xmax mccPSDmonMono_xmax
#define ymin mccPSDmonMono_ymin
#define ymax mccPSDmonMono_ymax
#define xwidth mccPSDmonMono_xwidth
#define yheight mccPSDmonMono_yheight
#define restore_neutron mccPSDmonMono_restore_neutron
#define nowritefile mccPSDmonMono_nowritefile
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
#line 17231 "PSI_Focus.c"
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

  /* Initializations for component MONOTOF. */
  SIG_MESSAGE("MONOTOF (Init)");
#define mccompcurname  MONOTOF
#define mccompcurtype  Monitor_nD
#define mccompcurindex 21
#define user1 mccMONOTOF_user1
#define user2 mccMONOTOF_user2
#define user3 mccMONOTOF_user3
#define DEFS mccMONOTOF_DEFS
#define Vars mccMONOTOF_Vars
#define detector mccMONOTOF_detector
#define offdata mccMONOTOF_offdata
#define xwidth mccMONOTOF_xwidth
#define yheight mccMONOTOF_yheight
#define zdepth mccMONOTOF_zdepth
#define xmin mccMONOTOF_xmin
#define xmax mccMONOTOF_xmax
#define ymin mccMONOTOF_ymin
#define ymax mccMONOTOF_ymax
#define zmin mccMONOTOF_zmin
#define zmax mccMONOTOF_zmax
#define bins mccMONOTOF_bins
#define min mccMONOTOF_min
#define max mccMONOTOF_max
#define restore_neutron mccMONOTOF_restore_neutron
#define radius mccMONOTOF_radius
#define options mccMONOTOF_options
#define filename mccMONOTOF_filename
#define geometry mccMONOTOF_geometry
#define username1 mccMONOTOF_username1
#define username2 mccMONOTOF_username2
#define username3 mccMONOTOF_username3
#define nowritefile mccMONOTOF_nowritefile
#line 229 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
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
#line 17362 "PSI_Focus.c"
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

  /* Initializations for component DivMonMono. */
  SIG_MESSAGE("DivMonMono (Init)");
#define mccompcurname  DivMonMono
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 22
#define nh mccDivMonMono_nh
#define nv mccDivMonMono_nv
#define Div_N mccDivMonMono_Div_N
#define Div_p mccDivMonMono_Div_p
#define Div_p2 mccDivMonMono_Div_p2
#define filename mccDivMonMono_filename
#define xmin mccDivMonMono_xmin
#define xmax mccDivMonMono_xmax
#define ymin mccDivMonMono_ymin
#define ymax mccDivMonMono_ymax
#define xwidth mccDivMonMono_xwidth
#define yheight mccDivMonMono_yheight
#define maxdiv_h mccDivMonMono_maxdiv_h
#define maxdiv_v mccDivMonMono_maxdiv_v
#define restore_neutron mccDivMonMono_restore_neutron
#define nx mccDivMonMono_nx
#define ny mccDivMonMono_ny
#define nz mccDivMonMono_nz
#define nowritefile mccDivMonMono_nowritefile
#line 66 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
{
int i,j;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("Divergence_monitor: %s: Null detection area !\n"
                   "ERROR               (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    for (i=0; i<nh; i++)
     for (j=0; j<nv; j++)
     {
      Div_N[i][j] = 0;
      Div_p[i][j] = 0;
      Div_p2[i][j] = 0;
     }
    NORM(nx,ny,nz);
}
#line 17442 "PSI_Focus.c"
#undef nowritefile
#undef nz
#undef ny
#undef nx
#undef restore_neutron
#undef maxdiv_v
#undef maxdiv_h
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
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component focus_mono. */
  SIG_MESSAGE("focus_mono (Init)");

  /* Initializations for component mono. */
  SIG_MESSAGE("mono (Init)");
#define mccompcurname  mono
#define mccompcurtype  Monochromator_2foc
#define mccompcurindex 24
#define mos_y mccmono_mos_y
#define mos_z mccmono_mos_z
#define mono_Q mccmono_mono_Q
#define SlabWidth mccmono_SlabWidth
#define SlabHeight mccmono_SlabHeight
#define rTable mccmono_rTable
#define reflect mccmono_reflect
#define zwidth mccmono_zwidth
#define yheight mccmono_yheight
#define gap mccmono_gap
#define NH mccmono_NH
#define NV mccmono_NV
#define mosaich mccmono_mosaich
#define mosaicv mccmono_mosaicv
#define r0 mccmono_r0
#define Q mccmono_Q
#define RV mccmono_RV
#define RH mccmono_RH
#define DM mccmono_DM
#define mosaic mccmono_mosaic
#define width mccmono_width
#define height mccmono_height
#define verbose mccmono_verbose
#line 99 "/usr/share/mcstas/2.5/contrib/Monochromator_2foc.comp"
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
#line 17545 "PSI_Focus.c"
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

  /* Initializations for component a2. */
  SIG_MESSAGE("a2 (Init)");

  /* Initializations for component FERMITOF_before. */
  SIG_MESSAGE("FERMITOF_before (Init)");
#define mccompcurname  FERMITOF_before
#define mccompcurtype  Monitor_nD
#define mccompcurindex 26
#define user1 mccFERMITOF_before_user1
#define user2 mccFERMITOF_before_user2
#define user3 mccFERMITOF_before_user3
#define DEFS mccFERMITOF_before_DEFS
#define Vars mccFERMITOF_before_Vars
#define detector mccFERMITOF_before_detector
#define offdata mccFERMITOF_before_offdata
#define xwidth mccFERMITOF_before_xwidth
#define yheight mccFERMITOF_before_yheight
#define zdepth mccFERMITOF_before_zdepth
#define xmin mccFERMITOF_before_xmin
#define xmax mccFERMITOF_before_xmax
#define ymin mccFERMITOF_before_ymin
#define ymax mccFERMITOF_before_ymax
#define zmin mccFERMITOF_before_zmin
#define zmax mccFERMITOF_before_zmax
#define bins mccFERMITOF_before_bins
#define min mccFERMITOF_before_min
#define max mccFERMITOF_before_max
#define restore_neutron mccFERMITOF_before_restore_neutron
#define radius mccFERMITOF_before_radius
#define options mccFERMITOF_before_options
#define filename mccFERMITOF_before_filename
#define geometry mccFERMITOF_before_geometry
#define username1 mccFERMITOF_before_username1
#define username2 mccFERMITOF_before_username2
#define username3 mccFERMITOF_before_username3
#define nowritefile mccFERMITOF_before_nowritefile
#line 229 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
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
#line 17688 "PSI_Focus.c"
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

  /* Initializations for component lambdaFermi. */
  SIG_MESSAGE("lambdaFermi (Init)");
#define mccompcurname  lambdaFermi
#define mccompcurtype  L_monitor
#define mccompcurindex 27
#define nL mcclambdaFermi_nL
#define L_N mcclambdaFermi_L_N
#define L_p mcclambdaFermi_L_p
#define L_p2 mcclambdaFermi_L_p2
#define filename mcclambdaFermi_filename
#define xmin mcclambdaFermi_xmin
#define xmax mcclambdaFermi_xmax
#define ymin mcclambdaFermi_ymin
#define ymax mcclambdaFermi_ymax
#define xwidth mcclambdaFermi_xwidth
#define yheight mcclambdaFermi_yheight
#define Lmin mcclambdaFermi_Lmin
#define Lmax mcclambdaFermi_Lmax
#define restore_neutron mcclambdaFermi_restore_neutron
#define nowritefile mcclambdaFermi_nowritefile
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
#line 17762 "PSI_Focus.c"
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

  /* Initializations for component EMON_Fermi. */
  SIG_MESSAGE("EMON_Fermi (Init)");
#define mccompcurname  EMON_Fermi
#define mccompcurtype  E_monitor
#define mccompcurindex 28
#define nE mccEMON_Fermi_nE
#define E_N mccEMON_Fermi_E_N
#define E_p mccEMON_Fermi_E_p
#define E_p2 mccEMON_Fermi_E_p2
#define S_p mccEMON_Fermi_S_p
#define S_pE mccEMON_Fermi_S_pE
#define S_pE2 mccEMON_Fermi_S_pE2
#define filename mccEMON_Fermi_filename
#define xmin mccEMON_Fermi_xmin
#define xmax mccEMON_Fermi_xmax
#define ymin mccEMON_Fermi_ymin
#define ymax mccEMON_Fermi_ymax
#define xwidth mccEMON_Fermi_xwidth
#define yheight mccEMON_Fermi_yheight
#define Emin mccEMON_Fermi_Emin
#define Emax mccEMON_Fermi_Emax
#define restore_neutron mccEMON_Fermi_restore_neutron
#define nowritefile mccEMON_Fermi_nowritefile
#line 66 "/usr/share/mcstas/2.5/monitors/E_monitor.comp"
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
#line 17827 "PSI_Focus.c"
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

  /* Initializations for component DivMonfermi1. */
  SIG_MESSAGE("DivMonfermi1 (Init)");
#define mccompcurname  DivMonfermi1
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 29
#define nh mccDivMonfermi1_nh
#define nv mccDivMonfermi1_nv
#define Div_N mccDivMonfermi1_Div_N
#define Div_p mccDivMonfermi1_Div_p
#define Div_p2 mccDivMonfermi1_Div_p2
#define filename mccDivMonfermi1_filename
#define xmin mccDivMonfermi1_xmin
#define xmax mccDivMonfermi1_xmax
#define ymin mccDivMonfermi1_ymin
#define ymax mccDivMonfermi1_ymax
#define xwidth mccDivMonfermi1_xwidth
#define yheight mccDivMonfermi1_yheight
#define maxdiv_h mccDivMonfermi1_maxdiv_h
#define maxdiv_v mccDivMonfermi1_maxdiv_v
#define restore_neutron mccDivMonfermi1_restore_neutron
#define nx mccDivMonfermi1_nx
#define ny mccDivMonfermi1_ny
#define nz mccDivMonfermi1_nz
#define nowritefile mccDivMonfermi1_nowritefile
#line 66 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
{
int i,j;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("Divergence_monitor: %s: Null detection area !\n"
                   "ERROR               (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    for (i=0; i<nh; i++)
     for (j=0; j<nv; j++)
     {
      Div_N[i][j] = 0;
      Div_p[i][j] = 0;
      Div_p2[i][j] = 0;
     }
    NORM(nx,ny,nz);
}
#line 17897 "PSI_Focus.c"
#undef nowritefile
#undef nz
#undef ny
#undef nx
#undef restore_neutron
#undef maxdiv_v
#undef maxdiv_h
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
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component PSD_Fermi1. */
  SIG_MESSAGE("PSD_Fermi1 (Init)");
#define mccompcurname  PSD_Fermi1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 30
#define nx mccPSD_Fermi1_nx
#define ny mccPSD_Fermi1_ny
#define PSD_N mccPSD_Fermi1_PSD_N
#define PSD_p mccPSD_Fermi1_PSD_p
#define PSD_p2 mccPSD_Fermi1_PSD_p2
#define filename mccPSD_Fermi1_filename
#define xmin mccPSD_Fermi1_xmin
#define xmax mccPSD_Fermi1_xmax
#define ymin mccPSD_Fermi1_ymin
#define ymax mccPSD_Fermi1_ymax
#define xwidth mccPSD_Fermi1_xwidth
#define yheight mccPSD_Fermi1_yheight
#define restore_neutron mccPSD_Fermi1_restore_neutron
#define nowritefile mccPSD_Fermi1_nowritefile
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
#line 17962 "PSI_Focus.c"
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

  /* Initializations for component FoChopper. */
  SIG_MESSAGE("FoChopper (Init)");
#define mccompcurname  FoChopper
#define mccompcurtype  FermiChopper
#define mccompcurindex 31
#define FCVars mccFoChopper_FCVars
#define phase mccFoChopper_phase
#define radius mccFoChopper_radius
#define nu mccFoChopper_nu
#define w mccFoChopper_w
#define nslit mccFoChopper_nslit
#define R0 mccFoChopper_R0
#define Qc mccFoChopper_Qc
#define alpha mccFoChopper_alpha
#define m mccFoChopper_m
#define W mccFoChopper_W
#define length mccFoChopper_length
#define eff mccFoChopper_eff
#define zero_time mccFoChopper_zero_time
#define xwidth mccFoChopper_xwidth
#define verbose mccFoChopper_verbose
#define yheight mccFoChopper_yheight
#define curvature mccFoChopper_curvature
#define delay mccFoChopper_delay
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
#line 18083 "PSI_Focus.c"
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

  /* Initializations for component PSD_Fermi2. */
  SIG_MESSAGE("PSD_Fermi2 (Init)");
#define mccompcurname  PSD_Fermi2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 32
#define nx mccPSD_Fermi2_nx
#define ny mccPSD_Fermi2_ny
#define PSD_N mccPSD_Fermi2_PSD_N
#define PSD_p mccPSD_Fermi2_PSD_p
#define PSD_p2 mccPSD_Fermi2_PSD_p2
#define filename mccPSD_Fermi2_filename
#define xmin mccPSD_Fermi2_xmin
#define xmax mccPSD_Fermi2_xmax
#define ymin mccPSD_Fermi2_ymin
#define ymax mccPSD_Fermi2_ymax
#define xwidth mccPSD_Fermi2_xwidth
#define yheight mccPSD_Fermi2_yheight
#define restore_neutron mccPSD_Fermi2_restore_neutron
#define nowritefile mccPSD_Fermi2_nowritefile
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
#line 18148 "PSI_Focus.c"
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

  /* Initializations for component DivMonfermi2. */
  SIG_MESSAGE("DivMonfermi2 (Init)");
#define mccompcurname  DivMonfermi2
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 33
#define nh mccDivMonfermi2_nh
#define nv mccDivMonfermi2_nv
#define Div_N mccDivMonfermi2_Div_N
#define Div_p mccDivMonfermi2_Div_p
#define Div_p2 mccDivMonfermi2_Div_p2
#define filename mccDivMonfermi2_filename
#define xmin mccDivMonfermi2_xmin
#define xmax mccDivMonfermi2_xmax
#define ymin mccDivMonfermi2_ymin
#define ymax mccDivMonfermi2_ymax
#define xwidth mccDivMonfermi2_xwidth
#define yheight mccDivMonfermi2_yheight
#define maxdiv_h mccDivMonfermi2_maxdiv_h
#define maxdiv_v mccDivMonfermi2_maxdiv_v
#define restore_neutron mccDivMonfermi2_restore_neutron
#define nx mccDivMonfermi2_nx
#define ny mccDivMonfermi2_ny
#define nz mccDivMonfermi2_nz
#define nowritefile mccDivMonfermi2_nowritefile
#line 66 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
{
int i,j;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("Divergence_monitor: %s: Null detection area !\n"
                   "ERROR               (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    for (i=0; i<nh; i++)
     for (j=0; j<nv; j++)
     {
      Div_N[i][j] = 0;
      Div_p[i][j] = 0;
      Div_p2[i][j] = 0;
     }
    NORM(nx,ny,nz);
}
#line 18214 "PSI_Focus.c"
#undef nowritefile
#undef nz
#undef ny
#undef nx
#undef restore_neutron
#undef maxdiv_v
#undef maxdiv_h
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
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component FERMITOF1. */
  SIG_MESSAGE("FERMITOF1 (Init)");
#define mccompcurname  FERMITOF1
#define mccompcurtype  Monitor_nD
#define mccompcurindex 34
#define user1 mccFERMITOF1_user1
#define user2 mccFERMITOF1_user2
#define user3 mccFERMITOF1_user3
#define DEFS mccFERMITOF1_DEFS
#define Vars mccFERMITOF1_Vars
#define detector mccFERMITOF1_detector
#define offdata mccFERMITOF1_offdata
#define xwidth mccFERMITOF1_xwidth
#define yheight mccFERMITOF1_yheight
#define zdepth mccFERMITOF1_zdepth
#define xmin mccFERMITOF1_xmin
#define xmax mccFERMITOF1_xmax
#define ymin mccFERMITOF1_ymin
#define ymax mccFERMITOF1_ymax
#define zmin mccFERMITOF1_zmin
#define zmax mccFERMITOF1_zmax
#define bins mccFERMITOF1_bins
#define min mccFERMITOF1_min
#define max mccFERMITOF1_max
#define restore_neutron mccFERMITOF1_restore_neutron
#define radius mccFERMITOF1_radius
#define options mccFERMITOF1_options
#define filename mccFERMITOF1_filename
#define geometry mccFERMITOF1_geometry
#define username1 mccFERMITOF1_username1
#define username2 mccFERMITOF1_username2
#define username3 mccFERMITOF1_username3
#define nowritefile mccFERMITOF1_nowritefile
#line 229 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
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
#line 18350 "PSI_Focus.c"
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

  /* Initializations for component SAMPLE_SLIT. */
  SIG_MESSAGE("SAMPLE_SLIT (Init)");
#define mccompcurname  SAMPLE_SLIT
#define mccompcurtype  Slit
#define mccompcurindex 35
#define xmin mccSAMPLE_SLIT_xmin
#define xmax mccSAMPLE_SLIT_xmax
#define ymin mccSAMPLE_SLIT_ymin
#define ymax mccSAMPLE_SLIT_ymax
#define radius mccSAMPLE_SLIT_radius
#define xwidth mccSAMPLE_SLIT_xwidth
#define yheight mccSAMPLE_SLIT_yheight
#line 50 "/usr/share/mcstas/2.5/optics/Slit.comp"
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
#line 18415 "PSI_Focus.c"
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

  /* Initializations for component FERMITOF2. */
  SIG_MESSAGE("FERMITOF2 (Init)");
#define mccompcurname  FERMITOF2
#define mccompcurtype  Monitor_nD
#define mccompcurindex 36
#define user1 mccFERMITOF2_user1
#define user2 mccFERMITOF2_user2
#define user3 mccFERMITOF2_user3
#define DEFS mccFERMITOF2_DEFS
#define Vars mccFERMITOF2_Vars
#define detector mccFERMITOF2_detector
#define offdata mccFERMITOF2_offdata
#define xwidth mccFERMITOF2_xwidth
#define yheight mccFERMITOF2_yheight
#define zdepth mccFERMITOF2_zdepth
#define xmin mccFERMITOF2_xmin
#define xmax mccFERMITOF2_xmax
#define ymin mccFERMITOF2_ymin
#define ymax mccFERMITOF2_ymax
#define zmin mccFERMITOF2_zmin
#define zmax mccFERMITOF2_zmax
#define bins mccFERMITOF2_bins
#define min mccFERMITOF2_min
#define max mccFERMITOF2_max
#define restore_neutron mccFERMITOF2_restore_neutron
#define radius mccFERMITOF2_radius
#define options mccFERMITOF2_options
#define filename mccFERMITOF2_filename
#define geometry mccFERMITOF2_geometry
#define username1 mccFERMITOF2_username1
#define username2 mccFERMITOF2_username2
#define username3 mccFERMITOF2_username3
#define nowritefile mccFERMITOF2_nowritefile
#line 229 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
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
#line 18539 "PSI_Focus.c"
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

  /* Initializations for component PSD_SAMPLE. */
  SIG_MESSAGE("PSD_SAMPLE (Init)");
#define mccompcurname  PSD_SAMPLE
#define mccompcurtype  PSD_monitor
#define mccompcurindex 37
#define nx mccPSD_SAMPLE_nx
#define ny mccPSD_SAMPLE_ny
#define PSD_N mccPSD_SAMPLE_PSD_N
#define PSD_p mccPSD_SAMPLE_PSD_p
#define PSD_p2 mccPSD_SAMPLE_PSD_p2
#define filename mccPSD_SAMPLE_filename
#define xmin mccPSD_SAMPLE_xmin
#define xmax mccPSD_SAMPLE_xmax
#define ymin mccPSD_SAMPLE_ymin
#define ymax mccPSD_SAMPLE_ymax
#define xwidth mccPSD_SAMPLE_xwidth
#define yheight mccPSD_SAMPLE_yheight
#define restore_neutron mccPSD_SAMPLE_restore_neutron
#define nowritefile mccPSD_SAMPLE_nowritefile
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
#line 18613 "PSI_Focus.c"
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

  /* Initializations for component DivMon_Sample. */
  SIG_MESSAGE("DivMon_Sample (Init)");
#define mccompcurname  DivMon_Sample
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 38
#define nh mccDivMon_Sample_nh
#define nv mccDivMon_Sample_nv
#define Div_N mccDivMon_Sample_Div_N
#define Div_p mccDivMon_Sample_Div_p
#define Div_p2 mccDivMon_Sample_Div_p2
#define filename mccDivMon_Sample_filename
#define xmin mccDivMon_Sample_xmin
#define xmax mccDivMon_Sample_xmax
#define ymin mccDivMon_Sample_ymin
#define ymax mccDivMon_Sample_ymax
#define xwidth mccDivMon_Sample_xwidth
#define yheight mccDivMon_Sample_yheight
#define maxdiv_h mccDivMon_Sample_maxdiv_h
#define maxdiv_v mccDivMon_Sample_maxdiv_v
#define restore_neutron mccDivMon_Sample_restore_neutron
#define nx mccDivMon_Sample_nx
#define ny mccDivMon_Sample_ny
#define nz mccDivMon_Sample_nz
#define nowritefile mccDivMon_Sample_nowritefile
#line 66 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
{
int i,j;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("Divergence_monitor: %s: Null detection area !\n"
                   "ERROR               (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    for (i=0; i<nh; i++)
     for (j=0; j<nv; j++)
     {
      Div_N[i][j] = 0;
      Div_p[i][j] = 0;
      Div_p2[i][j] = 0;
     }
    NORM(nx,ny,nz);
}
#line 18679 "PSI_Focus.c"
#undef nowritefile
#undef nz
#undef ny
#undef nx
#undef restore_neutron
#undef maxdiv_v
#undef maxdiv_h
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
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component EMON_SAMPLE. */
  SIG_MESSAGE("EMON_SAMPLE (Init)");
#define mccompcurname  EMON_SAMPLE
#define mccompcurtype  E_monitor
#define mccompcurindex 39
#define nE mccEMON_SAMPLE_nE
#define E_N mccEMON_SAMPLE_E_N
#define E_p mccEMON_SAMPLE_E_p
#define E_p2 mccEMON_SAMPLE_E_p2
#define S_p mccEMON_SAMPLE_S_p
#define S_pE mccEMON_SAMPLE_S_pE
#define S_pE2 mccEMON_SAMPLE_S_pE2
#define filename mccEMON_SAMPLE_filename
#define xmin mccEMON_SAMPLE_xmin
#define xmax mccEMON_SAMPLE_xmax
#define ymin mccEMON_SAMPLE_ymin
#define ymax mccEMON_SAMPLE_ymax
#define xwidth mccEMON_SAMPLE_xwidth
#define yheight mccEMON_SAMPLE_yheight
#define Emin mccEMON_SAMPLE_Emin
#define Emax mccEMON_SAMPLE_Emax
#define restore_neutron mccEMON_SAMPLE_restore_neutron
#define nowritefile mccEMON_SAMPLE_nowritefile
#line 66 "/usr/share/mcstas/2.5/monitors/E_monitor.comp"
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
#line 18748 "PSI_Focus.c"
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

  /* Initializations for component a3. */
  SIG_MESSAGE("a3 (Init)");

  /* Initializations for component Sample. */
  SIG_MESSAGE("Sample (Init)");
#define mccompcurname  Sample
#define mccompcurtype  V_sample
#define mccompcurindex 41
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
#line 18868 "PSI_Focus.c"
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

  /* Initializations for component TOF_Det. */
  SIG_MESSAGE("TOF_Det (Init)");
#define mccompcurname  TOF_Det
#define mccompcurtype  Monitor_nD
#define mccompcurindex 42
#define user1 mccTOF_Det_user1
#define user2 mccTOF_Det_user2
#define user3 mccTOF_Det_user3
#define DEFS mccTOF_Det_DEFS
#define Vars mccTOF_Det_Vars
#define detector mccTOF_Det_detector
#define offdata mccTOF_Det_offdata
#define xwidth mccTOF_Det_xwidth
#define yheight mccTOF_Det_yheight
#define zdepth mccTOF_Det_zdepth
#define xmin mccTOF_Det_xmin
#define xmax mccTOF_Det_xmax
#define ymin mccTOF_Det_ymin
#define ymax mccTOF_Det_ymax
#define zmin mccTOF_Det_zmin
#define zmax mccTOF_Det_zmax
#define bins mccTOF_Det_bins
#define min mccTOF_Det_min
#define max mccTOF_Det_max
#define restore_neutron mccTOF_Det_restore_neutron
#define radius mccTOF_Det_radius
#define options mccTOF_Det_options
#define filename mccTOF_Det_filename
#define geometry mccTOF_Det_geometry
#define username1 mccTOF_Det_username1
#define username2 mccTOF_Det_username2
#define username3 mccTOF_Det_username3
#define nowritefile mccTOF_Det_nowritefile
#line 229 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
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
#line 19016 "PSI_Focus.c"
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

  /* Initializations for component FoDet. */
  SIG_MESSAGE("FoDet (Init)");
#define mccompcurname  FoDet
#define mccompcurtype  Monitor_nD
#define mccompcurindex 43
#define user1 mccFoDet_user1
#define user2 mccFoDet_user2
#define user3 mccFoDet_user3
#define DEFS mccFoDet_DEFS
#define Vars mccFoDet_Vars
#define detector mccFoDet_detector
#define offdata mccFoDet_offdata
#define xwidth mccFoDet_xwidth
#define yheight mccFoDet_yheight
#define zdepth mccFoDet_zdepth
#define xmin mccFoDet_xmin
#define xmax mccFoDet_xmax
#define ymin mccFoDet_ymin
#define ymax mccFoDet_ymax
#define zmin mccFoDet_zmin
#define zmax mccFoDet_zmax
#define bins mccFoDet_bins
#define min mccFoDet_min
#define max mccFoDet_max
#define restore_neutron mccFoDet_restore_neutron
#define radius mccFoDet_radius
#define options mccFoDet_options
#define filename mccFoDet_filename
#define geometry mccFoDet_geometry
#define username1 mccFoDet_username1
#define username2 mccFoDet_username2
#define username3 mccFoDet_username3
#define nowritefile mccFoDet_nowritefile
#line 229 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
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
#line 19161 "PSI_Focus.c"
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

  /* Initializations for component EMON_DET. */
  SIG_MESSAGE("EMON_DET (Init)");
#define mccompcurname  EMON_DET
#define mccompcurtype  E_monitor
#define mccompcurindex 44
#define nE mccEMON_DET_nE
#define E_N mccEMON_DET_E_N
#define E_p mccEMON_DET_E_p
#define E_p2 mccEMON_DET_E_p2
#define S_p mccEMON_DET_S_p
#define S_pE mccEMON_DET_S_pE
#define S_pE2 mccEMON_DET_S_pE2
#define filename mccEMON_DET_filename
#define xmin mccEMON_DET_xmin
#define xmax mccEMON_DET_xmax
#define ymin mccEMON_DET_ymin
#define ymax mccEMON_DET_ymax
#define xwidth mccEMON_DET_xwidth
#define yheight mccEMON_DET_yheight
#define Emin mccEMON_DET_Emin
#define Emax mccEMON_DET_Emax
#define restore_neutron mccEMON_DET_restore_neutron
#define nowritefile mccEMON_DET_nowritefile
#line 66 "/usr/share/mcstas/2.5/monitors/E_monitor.comp"
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
#line 19239 "PSI_Focus.c"
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
  /* SPLIT counter for component mono */
  int mcSplit_mono=0;
  /* SPLIT counter for component Sample */
  int mcSplit_Sample=0;
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
#line 19419 "PSI_Focus.c"
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

  /* TRACE Component csource [2] */
  mccoordschange(mcposrcsource, mcrotrcsource,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component csource (without coords transformations) */
  mcJumpTrace_csource:
  SIG_MESSAGE("csource (Trace)");
  mcDEBUG_COMP("csource")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompcsource
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
#define mccompcurname  csource
#define mccompcurtype  Source_gen
#define mccompcurindex 2
#define p_in mcccsource_p_in
#define lambda1 mcccsource_lambda1
#define lambda2 mcccsource_lambda2
#define lambda3 mcccsource_lambda3
#define pTable mcccsource_pTable
#define pTable_x mcccsource_pTable_x
#define pTable_y mcccsource_pTable_y
#define pTable_xmin mcccsource_pTable_xmin
#define pTable_xmax mcccsource_pTable_xmax
#define pTable_xsum mcccsource_pTable_xsum
#define pTable_ymin mcccsource_pTable_ymin
#define pTable_ymax mcccsource_pTable_ymax
#define pTable_ysum mcccsource_pTable_ysum
#define pTable_dxmin mcccsource_pTable_dxmin
#define pTable_dxmax mcccsource_pTable_dxmax
#define pTable_dymin mcccsource_pTable_dymin
#define pTable_dymax mcccsource_pTable_dymax
{   /* Declarations of csource=Source_gen() SETTING parameters. */
char* flux_file = mcccsource_flux_file;
char* xdiv_file = mcccsource_xdiv_file;
char* ydiv_file = mcccsource_ydiv_file;
MCNUM radius = mcccsource_radius;
MCNUM dist = mcccsource_dist;
MCNUM focus_xw = mcccsource_focus_xw;
MCNUM focus_yh = mcccsource_focus_yh;
MCNUM focus_aw = mcccsource_focus_aw;
MCNUM focus_ah = mcccsource_focus_ah;
MCNUM E0 = mcccsource_E0;
MCNUM dE = mcccsource_dE;
MCNUM lambda0 = mcccsource_lambda0;
MCNUM dlambda = mcccsource_dlambda;
MCNUM I1 = mcccsource_I1;
MCNUM yheight = mcccsource_yheight;
MCNUM xwidth = mcccsource_xwidth;
MCNUM verbose = mcccsource_verbose;
MCNUM T1 = mcccsource_T1;
MCNUM flux_file_perAA = mcccsource_flux_file_perAA;
MCNUM flux_file_log = mcccsource_flux_file_log;
MCNUM Lmin = mcccsource_Lmin;
MCNUM Lmax = mcccsource_Lmax;
MCNUM Emin = mcccsource_Emin;
MCNUM Emax = mcccsource_Emax;
MCNUM T2 = mcccsource_T2;
MCNUM I2 = mcccsource_I2;
MCNUM T3 = mcccsource_T3;
MCNUM I3 = mcccsource_I3;
MCNUM zdepth = mcccsource_zdepth;
int target_index = mcccsource_target_index;
#line 479 "/usr/share/mcstas/2.5/sources/Source_gen.comp"
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
#line 19667 "PSI_Focus.c"
}   /* End of csource=Source_gen() SETTING parameter declarations. */
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
  mcabsorbCompcsource:
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

  /* TRACE Component guide1 [3] */
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
#define mccompcurname  guide1
#define mccompcurtype  Guide
#define mccompcurindex 3
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
#line 94 "/usr/share/mcstas/2.5/optics/Guide.comp"
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
#line 19909 "PSI_Focus.c"
}   /* End of guide1=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompguide1:
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

  /* TRACE Component guide2 [4] */
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
#define mccompcurname  guide2
#define mccompcurtype  Bender
#define mccompcurindex 4
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
#line 133 "/usr/share/mcstas/2.5/optics/Bender.comp"
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
#line 20175 "PSI_Focus.c"
}   /* End of guide2=Bender() SETTING parameter declarations. */
#undef mWin
#undef bk
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompguide2:
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

  /* TRACE Component bunker [5] */
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
#define mccompcurname  bunker
#define mccompcurtype  Guide
#define mccompcurindex 5
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
#line 94 "/usr/share/mcstas/2.5/optics/Guide.comp"
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
#line 20402 "PSI_Focus.c"
}   /* End of bunker=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompbunker:
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

  /* TRACE Component guide3 [6] */
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
#define mccompcurname  guide3
#define mccompcurtype  Guide
#define mccompcurindex 6
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
#line 94 "/usr/share/mcstas/2.5/optics/Guide.comp"
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
#line 20628 "PSI_Focus.c"
}   /* End of guide3=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompguide3:
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

  /* TRACE Component lambdaGuideExit [7] */
  mccoordschange(mcposrlambdaGuideExit, mcrotrlambdaGuideExit,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component lambdaGuideExit (without coords transformations) */
  mcJumpTrace_lambdaGuideExit:
  SIG_MESSAGE("lambdaGuideExit (Trace)");
  mcDEBUG_COMP("lambdaGuideExit")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComplambdaGuideExit
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
#define mccompcurname  lambdaGuideExit
#define mccompcurtype  L_monitor
#define mccompcurindex 7
#define nL mcclambdaGuideExit_nL
#define L_N mcclambdaGuideExit_L_N
#define L_p mcclambdaGuideExit_L_p
#define L_p2 mcclambdaGuideExit_L_p2
{   /* Declarations of lambdaGuideExit=L_monitor() SETTING parameters. */
char* filename = mcclambdaGuideExit_filename;
MCNUM xmin = mcclambdaGuideExit_xmin;
MCNUM xmax = mcclambdaGuideExit_xmax;
MCNUM ymin = mcclambdaGuideExit_ymin;
MCNUM ymax = mcclambdaGuideExit_ymax;
MCNUM xwidth = mcclambdaGuideExit_xwidth;
MCNUM yheight = mcclambdaGuideExit_yheight;
MCNUM Lmin = mcclambdaGuideExit_Lmin;
MCNUM Lmax = mcclambdaGuideExit_Lmax;
MCNUM restore_neutron = mcclambdaGuideExit_restore_neutron;
int nowritefile = mcclambdaGuideExit_nowritefile;
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
#line 20772 "PSI_Focus.c"
}   /* End of lambdaGuideExit=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComplambdaGuideExit:
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

  /* TRACE Component DivMonGuideExit [8] */
  mccoordschange(mcposrDivMonGuideExit, mcrotrDivMonGuideExit,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component DivMonGuideExit (without coords transformations) */
  mcJumpTrace_DivMonGuideExit:
  SIG_MESSAGE("DivMonGuideExit (Trace)");
  mcDEBUG_COMP("DivMonGuideExit")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompDivMonGuideExit
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
#define mccompcurname  DivMonGuideExit
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 8
#define nh mccDivMonGuideExit_nh
#define nv mccDivMonGuideExit_nv
#define Div_N mccDivMonGuideExit_Div_N
#define Div_p mccDivMonGuideExit_Div_p
#define Div_p2 mccDivMonGuideExit_Div_p2
{   /* Declarations of DivMonGuideExit=Divergence_monitor() SETTING parameters. */
char* filename = mccDivMonGuideExit_filename;
MCNUM xmin = mccDivMonGuideExit_xmin;
MCNUM xmax = mccDivMonGuideExit_xmax;
MCNUM ymin = mccDivMonGuideExit_ymin;
MCNUM ymax = mccDivMonGuideExit_ymax;
MCNUM xwidth = mccDivMonGuideExit_xwidth;
MCNUM yheight = mccDivMonGuideExit_yheight;
MCNUM maxdiv_h = mccDivMonGuideExit_maxdiv_h;
MCNUM maxdiv_v = mccDivMonGuideExit_maxdiv_v;
MCNUM restore_neutron = mccDivMonGuideExit_restore_neutron;
MCNUM nx = mccDivMonGuideExit_nx;
MCNUM ny = mccDivMonGuideExit_ny;
MCNUM nz = mccDivMonGuideExit_nz;
int nowritefile = mccDivMonGuideExit_nowritefile;
#line 89 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
{
    int i,j;
    double h_div, v_div;
    double v, vn;

    PROP_Z0;
    if (x>xmin && x<xmax && y>ymin && y<ymax)
    {
      /* Find length of projection onto the [nx ny nz] axis */
      vn = scalar_prod(vx, vy, vz, nx, ny, nz);
      h_div = RAD2DEG*atan2(vx,vn);
      v_div = RAD2DEG*atan2(vy,vn);
      if (h_div < maxdiv_h && h_div > -maxdiv_h &&
          v_div < maxdiv_v && v_div > -maxdiv_v)
      {
        i = floor((h_div + maxdiv_h)*nh/(2.0*maxdiv_h));
        j = floor((v_div + maxdiv_v)*nv/(2.0*maxdiv_v));
        Div_N[i][j]++;
        Div_p[i][j] += p;
        Div_p2[i][j] += p*p;
        SCATTER;
      }
    }
    if (restore_neutron) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
    }
}
#line 20929 "PSI_Focus.c"
}   /* End of DivMonGuideExit=Divergence_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompDivMonGuideExit:
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

  /* TRACE Component PSDGuideExit [9] */
  mccoordschange(mcposrPSDGuideExit, mcrotrPSDGuideExit,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PSDGuideExit (without coords transformations) */
  mcJumpTrace_PSDGuideExit:
  SIG_MESSAGE("PSDGuideExit (Trace)");
  mcDEBUG_COMP("PSDGuideExit")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompPSDGuideExit
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
#define mccompcurname  PSDGuideExit
#define mccompcurtype  PSD_monitor
#define mccompcurindex 9
#define nx mccPSDGuideExit_nx
#define ny mccPSDGuideExit_ny
#define PSD_N mccPSDGuideExit_PSD_N
#define PSD_p mccPSDGuideExit_PSD_p
#define PSD_p2 mccPSDGuideExit_PSD_p2
{   /* Declarations of PSDGuideExit=PSD_monitor() SETTING parameters. */
char* filename = mccPSDGuideExit_filename;
MCNUM xmin = mccPSDGuideExit_xmin;
MCNUM xmax = mccPSDGuideExit_xmax;
MCNUM ymin = mccPSDGuideExit_ymin;
MCNUM ymax = mccPSDGuideExit_ymax;
MCNUM xwidth = mccPSDGuideExit_xwidth;
MCNUM yheight = mccPSDGuideExit_yheight;
MCNUM restore_neutron = mccPSDGuideExit_restore_neutron;
int nowritefile = mccPSDGuideExit_nowritefile;
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
#line 21072 "PSI_Focus.c"
}   /* End of PSDGuideExit=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPSDGuideExit:
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

  /* TRACE Component FOCUSguide [10] */
  mccoordschange(mcposrFOCUSguide, mcrotrFOCUSguide,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component FOCUSguide (without coords transformations) */
  mcJumpTrace_FOCUSguide:
  SIG_MESSAGE("FOCUSguide (Trace)");
  mcDEBUG_COMP("FOCUSguide")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompFOCUSguide
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
#define mccompcurname  FOCUSguide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 10
#define w1c mccFOCUSguide_w1c
#define w2c mccFOCUSguide_w2c
#define ww mccFOCUSguide_ww
#define hh mccFOCUSguide_hh
#define whalf mccFOCUSguide_whalf
#define hhalf mccFOCUSguide_hhalf
#define lwhalf mccFOCUSguide_lwhalf
#define lhhalf mccFOCUSguide_lhhalf
{   /* Declarations of FOCUSguide=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccFOCUSguide_w1;
MCNUM h1 = mccFOCUSguide_h1;
MCNUM w2 = mccFOCUSguide_w2;
MCNUM h2 = mccFOCUSguide_h2;
MCNUM l = mccFOCUSguide_l;
MCNUM R0 = mccFOCUSguide_R0;
MCNUM Qc = mccFOCUSguide_Qc;
MCNUM alpha = mccFOCUSguide_alpha;
MCNUM m = mccFOCUSguide_m;
MCNUM nslit = mccFOCUSguide_nslit;
MCNUM d = mccFOCUSguide_d;
MCNUM Qcx = mccFOCUSguide_Qcx;
MCNUM Qcy = mccFOCUSguide_Qcy;
MCNUM alphax = mccFOCUSguide_alphax;
MCNUM alphay = mccFOCUSguide_alphay;
MCNUM W = mccFOCUSguide_W;
MCNUM mx = mccFOCUSguide_mx;
MCNUM my = mccFOCUSguide_my;
MCNUM nu = mccFOCUSguide_nu;
MCNUM phase = mccFOCUSguide_phase;
#line 131 "/usr/share/mcstas/2.5/optics/Guide_channeled.comp"
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
#line 21376 "PSI_Focus.c"
}   /* End of FOCUSguide=Guide_channeled() SETTING parameter declarations. */
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
  mcabsorbCompFOCUSguide:
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

  /* TRACE Component FirstChopper [11] */
  mccoordschange(mcposrFirstChopper, mcrotrFirstChopper,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component FirstChopper (without coords transformations) */
  mcJumpTrace_FirstChopper:
  SIG_MESSAGE("FirstChopper (Trace)");
  mcDEBUG_COMP("FirstChopper")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompFirstChopper
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
#define mccompcurname  FirstChopper
#define mccompcurtype  DiskChopper
#define mccompcurindex 11
#define Tg mccFirstChopper_Tg
#define To mccFirstChopper_To
#define delta_y mccFirstChopper_delta_y
#define height mccFirstChopper_height
#define omega mccFirstChopper_omega
{   /* Declarations of FirstChopper=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccFirstChopper_theta_0;
MCNUM radius = mccFirstChopper_radius;
MCNUM yheight = mccFirstChopper_yheight;
MCNUM nu = mccFirstChopper_nu;
MCNUM nslit = mccFirstChopper_nslit;
MCNUM jitter = mccFirstChopper_jitter;
MCNUM delay = mccFirstChopper_delay;
MCNUM isfirst = mccFirstChopper_isfirst;
MCNUM n_pulse = mccFirstChopper_n_pulse;
MCNUM abs_out = mccFirstChopper_abs_out;
MCNUM phase = mccFirstChopper_phase;
MCNUM xwidth = mccFirstChopper_xwidth;
MCNUM verbose = mccFirstChopper_verbose;
#line 133 "/usr/share/mcstas/2.5/optics/DiskChopper.comp"
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
#line 21542 "PSI_Focus.c"
}   /* End of FirstChopper=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompFirstChopper:
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

  /* TRACE Component DISCTOF [12] */
  mccoordschange(mcposrDISCTOF, mcrotrDISCTOF,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component DISCTOF (without coords transformations) */
  mcJumpTrace_DISCTOF:
  SIG_MESSAGE("DISCTOF (Trace)");
  mcDEBUG_COMP("DISCTOF")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompDISCTOF
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
#define mccompcurname  DISCTOF
#define mccompcurtype  Monitor_nD
#define mccompcurindex 12
#define user1 mccDISCTOF_user1
#define user2 mccDISCTOF_user2
#define user3 mccDISCTOF_user3
#define DEFS mccDISCTOF_DEFS
#define Vars mccDISCTOF_Vars
#define detector mccDISCTOF_detector
#define offdata mccDISCTOF_offdata
{   /* Declarations of DISCTOF=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccDISCTOF_xwidth;
MCNUM yheight = mccDISCTOF_yheight;
MCNUM zdepth = mccDISCTOF_zdepth;
MCNUM xmin = mccDISCTOF_xmin;
MCNUM xmax = mccDISCTOF_xmax;
MCNUM ymin = mccDISCTOF_ymin;
MCNUM ymax = mccDISCTOF_ymax;
MCNUM zmin = mccDISCTOF_zmin;
MCNUM zmax = mccDISCTOF_zmax;
MCNUM bins = mccDISCTOF_bins;
MCNUM min = mccDISCTOF_min;
MCNUM max = mccDISCTOF_max;
MCNUM restore_neutron = mccDISCTOF_restore_neutron;
MCNUM radius = mccDISCTOF_radius;
char* options = mccDISCTOF_options;
char* filename = mccDISCTOF_filename;
char* geometry = mccDISCTOF_geometry;
char* username1 = mccDISCTOF_username1;
char* username2 = mccDISCTOF_username2;
char* username3 = mccDISCTOF_username3;
int nowritefile = mccDISCTOF_nowritefile;
#line 309 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
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
#line 21850 "PSI_Focus.c"
}   /* End of DISCTOF=Monitor_nD() SETTING parameter declarations. */
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
  mcabsorbCompDISCTOF:
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

  /* TRACE Component PSDmon1Chopper [13] */
  mccoordschange(mcposrPSDmon1Chopper, mcrotrPSDmon1Chopper,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PSDmon1Chopper (without coords transformations) */
  mcJumpTrace_PSDmon1Chopper:
  SIG_MESSAGE("PSDmon1Chopper (Trace)");
  mcDEBUG_COMP("PSDmon1Chopper")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompPSDmon1Chopper
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
#define mccompcurname  PSDmon1Chopper
#define mccompcurtype  PSD_monitor
#define mccompcurindex 13
#define nx mccPSDmon1Chopper_nx
#define ny mccPSDmon1Chopper_ny
#define PSD_N mccPSDmon1Chopper_PSD_N
#define PSD_p mccPSDmon1Chopper_PSD_p
#define PSD_p2 mccPSDmon1Chopper_PSD_p2
{   /* Declarations of PSDmon1Chopper=PSD_monitor() SETTING parameters. */
char* filename = mccPSDmon1Chopper_filename;
MCNUM xmin = mccPSDmon1Chopper_xmin;
MCNUM xmax = mccPSDmon1Chopper_xmax;
MCNUM ymin = mccPSDmon1Chopper_ymin;
MCNUM ymax = mccPSDmon1Chopper_ymax;
MCNUM xwidth = mccPSDmon1Chopper_xwidth;
MCNUM yheight = mccPSDmon1Chopper_yheight;
MCNUM restore_neutron = mccPSDmon1Chopper_restore_neutron;
int nowritefile = mccPSDmon1Chopper_nowritefile;
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
#line 21995 "PSI_Focus.c"
}   /* End of PSDmon1Chopper=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPSDmon1Chopper:
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

  /* TRACE Component VacuumTube1entry [14] */
  mccoordschange(mcposrVacuumTube1entry, mcrotrVacuumTube1entry,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component VacuumTube1entry (without coords transformations) */
  mcJumpTrace_VacuumTube1entry:
  SIG_MESSAGE("VacuumTube1entry (Trace)");
  mcDEBUG_COMP("VacuumTube1entry")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompVacuumTube1entry
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
#define mccompcurname  VacuumTube1entry
#define mccompcurtype  Slit
#define mccompcurindex 14
{   /* Declarations of VacuumTube1entry=Slit() SETTING parameters. */
MCNUM xmin = mccVacuumTube1entry_xmin;
MCNUM xmax = mccVacuumTube1entry_xmax;
MCNUM ymin = mccVacuumTube1entry_ymin;
MCNUM ymax = mccVacuumTube1entry_ymax;
MCNUM radius = mccVacuumTube1entry_radius;
MCNUM xwidth = mccVacuumTube1entry_xwidth;
MCNUM yheight = mccVacuumTube1entry_yheight;
#line 71 "/usr/share/mcstas/2.5/optics/Slit.comp"
{
    PROP_Z0;
    if (((radius == 0) && (x<xmin || x>xmax || y<ymin || y>ymax))
    || ((radius != 0) && (x*x + y*y > radius*radius)))
      ABSORB;
    else
        SCATTER;
}
#line 22122 "PSI_Focus.c"
}   /* End of VacuumTube1entry=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompVacuumTube1entry:
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

  /* TRACE Component VacuumTube1exit [15] */
  mccoordschange(mcposrVacuumTube1exit, mcrotrVacuumTube1exit,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component VacuumTube1exit (without coords transformations) */
  mcJumpTrace_VacuumTube1exit:
  SIG_MESSAGE("VacuumTube1exit (Trace)");
  mcDEBUG_COMP("VacuumTube1exit")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompVacuumTube1exit
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
#define mccompcurname  VacuumTube1exit
#define mccompcurtype  Slit
#define mccompcurindex 15
{   /* Declarations of VacuumTube1exit=Slit() SETTING parameters. */
MCNUM xmin = mccVacuumTube1exit_xmin;
MCNUM xmax = mccVacuumTube1exit_xmax;
MCNUM ymin = mccVacuumTube1exit_ymin;
MCNUM ymax = mccVacuumTube1exit_ymax;
MCNUM radius = mccVacuumTube1exit_radius;
MCNUM xwidth = mccVacuumTube1exit_xwidth;
MCNUM yheight = mccVacuumTube1exit_yheight;
#line 71 "/usr/share/mcstas/2.5/optics/Slit.comp"
{
    PROP_Z0;
    if (((radius == 0) && (x<xmin || x>xmax || y<ymin || y>ymax))
    || ((radius != 0) && (x*x + y*y > radius*radius)))
      ABSORB;
    else
        SCATTER;
}
#line 22244 "PSI_Focus.c"
}   /* End of VacuumTube1exit=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompVacuumTube1exit:
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

  /* TRACE Component VacuumTube2entry [16] */
  mccoordschange(mcposrVacuumTube2entry, mcrotrVacuumTube2entry,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component VacuumTube2entry (without coords transformations) */
  mcJumpTrace_VacuumTube2entry:
  SIG_MESSAGE("VacuumTube2entry (Trace)");
  mcDEBUG_COMP("VacuumTube2entry")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompVacuumTube2entry
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
#define mccompcurname  VacuumTube2entry
#define mccompcurtype  Slit
#define mccompcurindex 16
{   /* Declarations of VacuumTube2entry=Slit() SETTING parameters. */
MCNUM xmin = mccVacuumTube2entry_xmin;
MCNUM xmax = mccVacuumTube2entry_xmax;
MCNUM ymin = mccVacuumTube2entry_ymin;
MCNUM ymax = mccVacuumTube2entry_ymax;
MCNUM radius = mccVacuumTube2entry_radius;
MCNUM xwidth = mccVacuumTube2entry_xwidth;
MCNUM yheight = mccVacuumTube2entry_yheight;
#line 71 "/usr/share/mcstas/2.5/optics/Slit.comp"
{
    PROP_Z0;
    if (((radius == 0) && (x<xmin || x>xmax || y<ymin || y>ymax))
    || ((radius != 0) && (x*x + y*y > radius*radius)))
      ABSORB;
    else
        SCATTER;
}
#line 22366 "PSI_Focus.c"
}   /* End of VacuumTube2entry=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompVacuumTube2entry:
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

  /* TRACE Component VacuumTube2exit [17] */
  mccoordschange(mcposrVacuumTube2exit, mcrotrVacuumTube2exit,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component VacuumTube2exit (without coords transformations) */
  mcJumpTrace_VacuumTube2exit:
  SIG_MESSAGE("VacuumTube2exit (Trace)");
  mcDEBUG_COMP("VacuumTube2exit")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompVacuumTube2exit
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
#define mccompcurname  VacuumTube2exit
#define mccompcurtype  Slit
#define mccompcurindex 17
{   /* Declarations of VacuumTube2exit=Slit() SETTING parameters. */
MCNUM xmin = mccVacuumTube2exit_xmin;
MCNUM xmax = mccVacuumTube2exit_xmax;
MCNUM ymin = mccVacuumTube2exit_ymin;
MCNUM ymax = mccVacuumTube2exit_ymax;
MCNUM radius = mccVacuumTube2exit_radius;
MCNUM xwidth = mccVacuumTube2exit_xwidth;
MCNUM yheight = mccVacuumTube2exit_yheight;
#line 71 "/usr/share/mcstas/2.5/optics/Slit.comp"
{
    PROP_Z0;
    if (((radius == 0) && (x<xmin || x>xmax || y<ymin || y>ymax))
    || ((radius != 0) && (x*x + y*y > radius*radius)))
      ABSORB;
    else
        SCATTER;
}
#line 22488 "PSI_Focus.c"
}   /* End of VacuumTube2exit=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompVacuumTube2exit:
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

  /* TRACE Component VacuumTube3entry [18] */
  mccoordschange(mcposrVacuumTube3entry, mcrotrVacuumTube3entry,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component VacuumTube3entry (without coords transformations) */
  mcJumpTrace_VacuumTube3entry:
  SIG_MESSAGE("VacuumTube3entry (Trace)");
  mcDEBUG_COMP("VacuumTube3entry")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompVacuumTube3entry
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
#define mccompcurname  VacuumTube3entry
#define mccompcurtype  Slit
#define mccompcurindex 18
{   /* Declarations of VacuumTube3entry=Slit() SETTING parameters. */
MCNUM xmin = mccVacuumTube3entry_xmin;
MCNUM xmax = mccVacuumTube3entry_xmax;
MCNUM ymin = mccVacuumTube3entry_ymin;
MCNUM ymax = mccVacuumTube3entry_ymax;
MCNUM radius = mccVacuumTube3entry_radius;
MCNUM xwidth = mccVacuumTube3entry_xwidth;
MCNUM yheight = mccVacuumTube3entry_yheight;
#line 71 "/usr/share/mcstas/2.5/optics/Slit.comp"
{
    PROP_Z0;
    if (((radius == 0) && (x<xmin || x>xmax || y<ymin || y>ymax))
    || ((radius != 0) && (x*x + y*y > radius*radius)))
      ABSORB;
    else
        SCATTER;
}
#line 22610 "PSI_Focus.c"
}   /* End of VacuumTube3entry=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompVacuumTube3entry:
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

  /* TRACE Component VacuumTube3exit [19] */
  mccoordschange(mcposrVacuumTube3exit, mcrotrVacuumTube3exit,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component VacuumTube3exit (without coords transformations) */
  mcJumpTrace_VacuumTube3exit:
  SIG_MESSAGE("VacuumTube3exit (Trace)");
  mcDEBUG_COMP("VacuumTube3exit")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompVacuumTube3exit
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
#define mccompcurname  VacuumTube3exit
#define mccompcurtype  Slit
#define mccompcurindex 19
{   /* Declarations of VacuumTube3exit=Slit() SETTING parameters. */
MCNUM xmin = mccVacuumTube3exit_xmin;
MCNUM xmax = mccVacuumTube3exit_xmax;
MCNUM ymin = mccVacuumTube3exit_ymin;
MCNUM ymax = mccVacuumTube3exit_ymax;
MCNUM radius = mccVacuumTube3exit_radius;
MCNUM xwidth = mccVacuumTube3exit_xwidth;
MCNUM yheight = mccVacuumTube3exit_yheight;
#line 71 "/usr/share/mcstas/2.5/optics/Slit.comp"
{
    PROP_Z0;
    if (((radius == 0) && (x<xmin || x>xmax || y<ymin || y>ymax))
    || ((radius != 0) && (x*x + y*y > radius*radius)))
      ABSORB;
    else
        SCATTER;
}
#line 22732 "PSI_Focus.c"
}   /* End of VacuumTube3exit=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompVacuumTube3exit:
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

  /* TRACE Component PSDmonMono [20] */
  mccoordschange(mcposrPSDmonMono, mcrotrPSDmonMono,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PSDmonMono (without coords transformations) */
  mcJumpTrace_PSDmonMono:
  SIG_MESSAGE("PSDmonMono (Trace)");
  mcDEBUG_COMP("PSDmonMono")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompPSDmonMono
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
#define mccompcurname  PSDmonMono
#define mccompcurtype  PSD_monitor
#define mccompcurindex 20
#define nx mccPSDmonMono_nx
#define ny mccPSDmonMono_ny
#define PSD_N mccPSDmonMono_PSD_N
#define PSD_p mccPSDmonMono_PSD_p
#define PSD_p2 mccPSDmonMono_PSD_p2
{   /* Declarations of PSDmonMono=PSD_monitor() SETTING parameters. */
char* filename = mccPSDmonMono_filename;
MCNUM xmin = mccPSDmonMono_xmin;
MCNUM xmax = mccPSDmonMono_xmax;
MCNUM ymin = mccPSDmonMono_ymin;
MCNUM ymax = mccPSDmonMono_ymax;
MCNUM xwidth = mccPSDmonMono_xwidth;
MCNUM yheight = mccPSDmonMono_yheight;
MCNUM restore_neutron = mccPSDmonMono_restore_neutron;
int nowritefile = mccPSDmonMono_nowritefile;
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
#line 22870 "PSI_Focus.c"
}   /* End of PSDmonMono=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPSDmonMono:
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

  /* TRACE Component MONOTOF [21] */
  mccoordschange(mcposrMONOTOF, mcrotrMONOTOF,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component MONOTOF (without coords transformations) */
  mcJumpTrace_MONOTOF:
  SIG_MESSAGE("MONOTOF (Trace)");
  mcDEBUG_COMP("MONOTOF")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompMONOTOF
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
#define mccompcurname  MONOTOF
#define mccompcurtype  Monitor_nD
#define mccompcurindex 21
#define user1 mccMONOTOF_user1
#define user2 mccMONOTOF_user2
#define user3 mccMONOTOF_user3
#define DEFS mccMONOTOF_DEFS
#define Vars mccMONOTOF_Vars
#define detector mccMONOTOF_detector
#define offdata mccMONOTOF_offdata
{   /* Declarations of MONOTOF=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccMONOTOF_xwidth;
MCNUM yheight = mccMONOTOF_yheight;
MCNUM zdepth = mccMONOTOF_zdepth;
MCNUM xmin = mccMONOTOF_xmin;
MCNUM xmax = mccMONOTOF_xmax;
MCNUM ymin = mccMONOTOF_ymin;
MCNUM ymax = mccMONOTOF_ymax;
MCNUM zmin = mccMONOTOF_zmin;
MCNUM zmax = mccMONOTOF_zmax;
MCNUM bins = mccMONOTOF_bins;
MCNUM min = mccMONOTOF_min;
MCNUM max = mccMONOTOF_max;
MCNUM restore_neutron = mccMONOTOF_restore_neutron;
MCNUM radius = mccMONOTOF_radius;
char* options = mccMONOTOF_options;
char* filename = mccMONOTOF_filename;
char* geometry = mccMONOTOF_geometry;
char* username1 = mccMONOTOF_username1;
char* username2 = mccMONOTOF_username2;
char* username3 = mccMONOTOF_username3;
int nowritefile = mccMONOTOF_nowritefile;
#line 309 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
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
#line 23178 "PSI_Focus.c"
}   /* End of MONOTOF=Monitor_nD() SETTING parameter declarations. */
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
  mcabsorbCompMONOTOF:
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

  /* TRACE Component DivMonMono [22] */
  mccoordschange(mcposrDivMonMono, mcrotrDivMonMono,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component DivMonMono (without coords transformations) */
  mcJumpTrace_DivMonMono:
  SIG_MESSAGE("DivMonMono (Trace)");
  mcDEBUG_COMP("DivMonMono")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompDivMonMono
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
#define mccompcurname  DivMonMono
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 22
#define nh mccDivMonMono_nh
#define nv mccDivMonMono_nv
#define Div_N mccDivMonMono_Div_N
#define Div_p mccDivMonMono_Div_p
#define Div_p2 mccDivMonMono_Div_p2
{   /* Declarations of DivMonMono=Divergence_monitor() SETTING parameters. */
char* filename = mccDivMonMono_filename;
MCNUM xmin = mccDivMonMono_xmin;
MCNUM xmax = mccDivMonMono_xmax;
MCNUM ymin = mccDivMonMono_ymin;
MCNUM ymax = mccDivMonMono_ymax;
MCNUM xwidth = mccDivMonMono_xwidth;
MCNUM yheight = mccDivMonMono_yheight;
MCNUM maxdiv_h = mccDivMonMono_maxdiv_h;
MCNUM maxdiv_v = mccDivMonMono_maxdiv_v;
MCNUM restore_neutron = mccDivMonMono_restore_neutron;
MCNUM nx = mccDivMonMono_nx;
MCNUM ny = mccDivMonMono_ny;
MCNUM nz = mccDivMonMono_nz;
int nowritefile = mccDivMonMono_nowritefile;
#line 89 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
{
    int i,j;
    double h_div, v_div;
    double v, vn;

    PROP_Z0;
    if (x>xmin && x<xmax && y>ymin && y<ymax)
    {
      /* Find length of projection onto the [nx ny nz] axis */
      vn = scalar_prod(vx, vy, vz, nx, ny, nz);
      h_div = RAD2DEG*atan2(vx,vn);
      v_div = RAD2DEG*atan2(vy,vn);
      if (h_div < maxdiv_h && h_div > -maxdiv_h &&
          v_div < maxdiv_v && v_div > -maxdiv_v)
      {
        i = floor((h_div + maxdiv_h)*nh/(2.0*maxdiv_h));
        j = floor((v_div + maxdiv_v)*nv/(2.0*maxdiv_v));
        Div_N[i][j]++;
        Div_p[i][j] += p;
        Div_p2[i][j] += p*p;
        SCATTER;
      }
    }
    if (restore_neutron) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
    }
}
#line 23338 "PSI_Focus.c"
}   /* End of DivMonMono=Divergence_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompDivMonMono:
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

  /* TRACE Component focus_mono [23] */
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
#define mccompcurname  focus_mono
#define mccompcurtype  Arm
#define mccompcurindex 23
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompfocus_mono:
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

  /* TRACE Component mono [24] */
  mccoordschange(mcposrmono, mcrotrmono,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component mono (without coords transformations) */
  mcJumpTrace_mono:
  SIG_MESSAGE("mono (Trace)");
  mcDEBUG_COMP("mono")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompmono
  if (!mcSplit_mono) {                   /* STORE only the first time */
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
  mcSplit_mono++; /* SPLIT number */
  mcScattered=0;
  mcRestore=0;
  mcNCounter[24]++;
  mcPCounter[24] += p;
  mcP2Counter[24] += p*p;
#define mccompcurname  mono
#define mccompcurtype  Monochromator_2foc
#define mccompcurindex 24
#define mos_y mccmono_mos_y
#define mos_z mccmono_mos_z
#define mono_Q mccmono_mono_Q
#define SlabWidth mccmono_SlabWidth
#define SlabHeight mccmono_SlabHeight
#define rTable mccmono_rTable
{   /* Declarations of mono=Monochromator_2foc() SETTING parameters. */
char* reflect = mccmono_reflect;
MCNUM zwidth = mccmono_zwidth;
MCNUM yheight = mccmono_yheight;
MCNUM gap = mccmono_gap;
MCNUM NH = mccmono_NH;
MCNUM NV = mccmono_NV;
MCNUM mosaich = mccmono_mosaich;
MCNUM mosaicv = mccmono_mosaicv;
MCNUM r0 = mccmono_r0;
MCNUM Q = mccmono_Q;
MCNUM RV = mccmono_RV;
MCNUM RH = mccmono_RH;
MCNUM DM = mccmono_DM;
MCNUM mosaic = mccmono_mosaic;
MCNUM width = mccmono_width;
MCNUM height = mccmono_height;
MCNUM verbose = mccmono_verbose;
#line 148 "/usr/share/mcstas/2.5/contrib/Monochromator_2foc.comp"
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
#line 23713 "PSI_Focus.c"
}   /* End of mono=Monochromator_2foc() SETTING parameter declarations. */
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
  mcabsorbCompmono:
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

  /* TRACE Component a2 [25] */
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
#define mccompcurname  a2
#define mccompcurtype  Arm
#define mccompcurindex 25
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompa2:
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

  /* TRACE Component FERMITOF_before [26] */
  mccoordschange(mcposrFERMITOF_before, mcrotrFERMITOF_before,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component FERMITOF_before (without coords transformations) */
  mcJumpTrace_FERMITOF_before:
  SIG_MESSAGE("FERMITOF_before (Trace)");
  mcDEBUG_COMP("FERMITOF_before")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompFERMITOF_before
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
#define mccompcurname  FERMITOF_before
#define mccompcurtype  Monitor_nD
#define mccompcurindex 26
#define user1 mccFERMITOF_before_user1
#define user2 mccFERMITOF_before_user2
#define user3 mccFERMITOF_before_user3
#define DEFS mccFERMITOF_before_DEFS
#define Vars mccFERMITOF_before_Vars
#define detector mccFERMITOF_before_detector
#define offdata mccFERMITOF_before_offdata
{   /* Declarations of FERMITOF_before=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccFERMITOF_before_xwidth;
MCNUM yheight = mccFERMITOF_before_yheight;
MCNUM zdepth = mccFERMITOF_before_zdepth;
MCNUM xmin = mccFERMITOF_before_xmin;
MCNUM xmax = mccFERMITOF_before_xmax;
MCNUM ymin = mccFERMITOF_before_ymin;
MCNUM ymax = mccFERMITOF_before_ymax;
MCNUM zmin = mccFERMITOF_before_zmin;
MCNUM zmax = mccFERMITOF_before_zmax;
MCNUM bins = mccFERMITOF_before_bins;
MCNUM min = mccFERMITOF_before_min;
MCNUM max = mccFERMITOF_before_max;
MCNUM restore_neutron = mccFERMITOF_before_restore_neutron;
MCNUM radius = mccFERMITOF_before_radius;
char* options = mccFERMITOF_before_options;
char* filename = mccFERMITOF_before_filename;
char* geometry = mccFERMITOF_before_geometry;
char* username1 = mccFERMITOF_before_username1;
char* username2 = mccFERMITOF_before_username2;
char* username3 = mccFERMITOF_before_username3;
int nowritefile = mccFERMITOF_before_nowritefile;
#line 309 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
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
#line 24125 "PSI_Focus.c"
}   /* End of FERMITOF_before=Monitor_nD() SETTING parameter declarations. */
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
  mcabsorbCompFERMITOF_before:
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

  /* TRACE Component lambdaFermi [27] */
  mccoordschange(mcposrlambdaFermi, mcrotrlambdaFermi,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component lambdaFermi (without coords transformations) */
  mcJumpTrace_lambdaFermi:
  SIG_MESSAGE("lambdaFermi (Trace)");
  mcDEBUG_COMP("lambdaFermi")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComplambdaFermi
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
#define mccompcurname  lambdaFermi
#define mccompcurtype  L_monitor
#define mccompcurindex 27
#define nL mcclambdaFermi_nL
#define L_N mcclambdaFermi_L_N
#define L_p mcclambdaFermi_L_p
#define L_p2 mcclambdaFermi_L_p2
{   /* Declarations of lambdaFermi=L_monitor() SETTING parameters. */
char* filename = mcclambdaFermi_filename;
MCNUM xmin = mcclambdaFermi_xmin;
MCNUM xmax = mcclambdaFermi_xmax;
MCNUM ymin = mcclambdaFermi_ymin;
MCNUM ymax = mcclambdaFermi_ymax;
MCNUM xwidth = mcclambdaFermi_xwidth;
MCNUM yheight = mcclambdaFermi_yheight;
MCNUM Lmin = mcclambdaFermi_Lmin;
MCNUM Lmax = mcclambdaFermi_Lmax;
MCNUM restore_neutron = mcclambdaFermi_restore_neutron;
int nowritefile = mcclambdaFermi_nowritefile;
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
#line 24275 "PSI_Focus.c"
}   /* End of lambdaFermi=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComplambdaFermi:
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

  /* TRACE Component EMON_Fermi [28] */
  mccoordschange(mcposrEMON_Fermi, mcrotrEMON_Fermi,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component EMON_Fermi (without coords transformations) */
  mcJumpTrace_EMON_Fermi:
  SIG_MESSAGE("EMON_Fermi (Trace)");
  mcDEBUG_COMP("EMON_Fermi")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompEMON_Fermi
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
#define mccompcurname  EMON_Fermi
#define mccompcurtype  E_monitor
#define mccompcurindex 28
#define nE mccEMON_Fermi_nE
#define E_N mccEMON_Fermi_E_N
#define E_p mccEMON_Fermi_E_p
#define E_p2 mccEMON_Fermi_E_p2
#define S_p mccEMON_Fermi_S_p
#define S_pE mccEMON_Fermi_S_pE
#define S_pE2 mccEMON_Fermi_S_pE2
{   /* Declarations of EMON_Fermi=E_monitor() SETTING parameters. */
char* filename = mccEMON_Fermi_filename;
MCNUM xmin = mccEMON_Fermi_xmin;
MCNUM xmax = mccEMON_Fermi_xmax;
MCNUM ymin = mccEMON_Fermi_ymin;
MCNUM ymax = mccEMON_Fermi_ymax;
MCNUM xwidth = mccEMON_Fermi_xwidth;
MCNUM yheight = mccEMON_Fermi_yheight;
MCNUM Emin = mccEMON_Fermi_Emin;
MCNUM Emax = mccEMON_Fermi_Emax;
MCNUM restore_neutron = mccEMON_Fermi_restore_neutron;
int nowritefile = mccEMON_Fermi_nowritefile;
#line 89 "/usr/share/mcstas/2.5/monitors/E_monitor.comp"
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
#line 24430 "PSI_Focus.c"
}   /* End of EMON_Fermi=E_monitor() SETTING parameter declarations. */
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
  mcabsorbCompEMON_Fermi:
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

  /* TRACE Component DivMonfermi1 [29] */
  mccoordschange(mcposrDivMonfermi1, mcrotrDivMonfermi1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component DivMonfermi1 (without coords transformations) */
  mcJumpTrace_DivMonfermi1:
  SIG_MESSAGE("DivMonfermi1 (Trace)");
  mcDEBUG_COMP("DivMonfermi1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompDivMonfermi1
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
#define mccompcurname  DivMonfermi1
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 29
#define nh mccDivMonfermi1_nh
#define nv mccDivMonfermi1_nv
#define Div_N mccDivMonfermi1_Div_N
#define Div_p mccDivMonfermi1_Div_p
#define Div_p2 mccDivMonfermi1_Div_p2
{   /* Declarations of DivMonfermi1=Divergence_monitor() SETTING parameters. */
char* filename = mccDivMonfermi1_filename;
MCNUM xmin = mccDivMonfermi1_xmin;
MCNUM xmax = mccDivMonfermi1_xmax;
MCNUM ymin = mccDivMonfermi1_ymin;
MCNUM ymax = mccDivMonfermi1_ymax;
MCNUM xwidth = mccDivMonfermi1_xwidth;
MCNUM yheight = mccDivMonfermi1_yheight;
MCNUM maxdiv_h = mccDivMonfermi1_maxdiv_h;
MCNUM maxdiv_v = mccDivMonfermi1_maxdiv_v;
MCNUM restore_neutron = mccDivMonfermi1_restore_neutron;
MCNUM nx = mccDivMonfermi1_nx;
MCNUM ny = mccDivMonfermi1_ny;
MCNUM nz = mccDivMonfermi1_nz;
int nowritefile = mccDivMonfermi1_nowritefile;
#line 89 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
{
    int i,j;
    double h_div, v_div;
    double v, vn;

    PROP_Z0;
    if (x>xmin && x<xmax && y>ymin && y<ymax)
    {
      /* Find length of projection onto the [nx ny nz] axis */
      vn = scalar_prod(vx, vy, vz, nx, ny, nz);
      h_div = RAD2DEG*atan2(vx,vn);
      v_div = RAD2DEG*atan2(vy,vn);
      if (h_div < maxdiv_h && h_div > -maxdiv_h &&
          v_div < maxdiv_v && v_div > -maxdiv_v)
      {
        i = floor((h_div + maxdiv_h)*nh/(2.0*maxdiv_h));
        j = floor((v_div + maxdiv_v)*nv/(2.0*maxdiv_v));
        Div_N[i][j]++;
        Div_p[i][j] += p;
        Div_p2[i][j] += p*p;
        SCATTER;
      }
    }
    if (restore_neutron) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
    }
}
#line 24590 "PSI_Focus.c"
}   /* End of DivMonfermi1=Divergence_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompDivMonfermi1:
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

  /* TRACE Component PSD_Fermi1 [30] */
  mccoordschange(mcposrPSD_Fermi1, mcrotrPSD_Fermi1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PSD_Fermi1 (without coords transformations) */
  mcJumpTrace_PSD_Fermi1:
  SIG_MESSAGE("PSD_Fermi1 (Trace)");
  mcDEBUG_COMP("PSD_Fermi1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompPSD_Fermi1
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
#define mccompcurname  PSD_Fermi1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 30
#define nx mccPSD_Fermi1_nx
#define ny mccPSD_Fermi1_ny
#define PSD_N mccPSD_Fermi1_PSD_N
#define PSD_p mccPSD_Fermi1_PSD_p
#define PSD_p2 mccPSD_Fermi1_PSD_p2
{   /* Declarations of PSD_Fermi1=PSD_monitor() SETTING parameters. */
char* filename = mccPSD_Fermi1_filename;
MCNUM xmin = mccPSD_Fermi1_xmin;
MCNUM xmax = mccPSD_Fermi1_xmax;
MCNUM ymin = mccPSD_Fermi1_ymin;
MCNUM ymax = mccPSD_Fermi1_ymax;
MCNUM xwidth = mccPSD_Fermi1_xwidth;
MCNUM yheight = mccPSD_Fermi1_yheight;
MCNUM restore_neutron = mccPSD_Fermi1_restore_neutron;
int nowritefile = mccPSD_Fermi1_nowritefile;
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
#line 24733 "PSI_Focus.c"
}   /* End of PSD_Fermi1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPSD_Fermi1:
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

  /* TRACE Component FoChopper [31] */
  mccoordschange(mcposrFoChopper, mcrotrFoChopper,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component FoChopper (without coords transformations) */
  mcJumpTrace_FoChopper:
  SIG_MESSAGE("FoChopper (Trace)");
  mcDEBUG_COMP("FoChopper")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompFoChopper
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
#define mccompcurname  FoChopper
#define mccompcurtype  FermiChopper
#define mccompcurindex 31
#define FCVars mccFoChopper_FCVars
{   /* Declarations of FoChopper=FermiChopper() SETTING parameters. */
MCNUM phase = mccFoChopper_phase;
MCNUM radius = mccFoChopper_radius;
MCNUM nu = mccFoChopper_nu;
MCNUM w = mccFoChopper_w;
MCNUM nslit = mccFoChopper_nslit;
MCNUM R0 = mccFoChopper_R0;
MCNUM Qc = mccFoChopper_Qc;
MCNUM alpha = mccFoChopper_alpha;
MCNUM m = mccFoChopper_m;
MCNUM W = mccFoChopper_W;
MCNUM length = mccFoChopper_length;
MCNUM eff = mccFoChopper_eff;
MCNUM zero_time = mccFoChopper_zero_time;
MCNUM xwidth = mccFoChopper_xwidth;
MCNUM verbose = mccFoChopper_verbose;
MCNUM yheight = mccFoChopper_yheight;
MCNUM curvature = mccFoChopper_curvature;
MCNUM delay = mccFoChopper_delay;
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
#line 25286 "PSI_Focus.c"
}   /* End of FoChopper=FermiChopper() SETTING parameter declarations. */
#undef FCVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompFoChopper:
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

  /* TRACE Component PSD_Fermi2 [32] */
  mccoordschange(mcposrPSD_Fermi2, mcrotrPSD_Fermi2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PSD_Fermi2 (without coords transformations) */
  mcJumpTrace_PSD_Fermi2:
  SIG_MESSAGE("PSD_Fermi2 (Trace)");
  mcDEBUG_COMP("PSD_Fermi2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompPSD_Fermi2
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
#define mccompcurname  PSD_Fermi2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 32
#define nx mccPSD_Fermi2_nx
#define ny mccPSD_Fermi2_ny
#define PSD_N mccPSD_Fermi2_PSD_N
#define PSD_p mccPSD_Fermi2_PSD_p
#define PSD_p2 mccPSD_Fermi2_PSD_p2
{   /* Declarations of PSD_Fermi2=PSD_monitor() SETTING parameters. */
char* filename = mccPSD_Fermi2_filename;
MCNUM xmin = mccPSD_Fermi2_xmin;
MCNUM xmax = mccPSD_Fermi2_xmax;
MCNUM ymin = mccPSD_Fermi2_ymin;
MCNUM ymax = mccPSD_Fermi2_ymax;
MCNUM xwidth = mccPSD_Fermi2_xwidth;
MCNUM yheight = mccPSD_Fermi2_yheight;
MCNUM restore_neutron = mccPSD_Fermi2_restore_neutron;
int nowritefile = mccPSD_Fermi2_nowritefile;
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
#line 25425 "PSI_Focus.c"
}   /* End of PSD_Fermi2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPSD_Fermi2:
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

  /* TRACE Component DivMonfermi2 [33] */
  mccoordschange(mcposrDivMonfermi2, mcrotrDivMonfermi2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component DivMonfermi2 (without coords transformations) */
  mcJumpTrace_DivMonfermi2:
  SIG_MESSAGE("DivMonfermi2 (Trace)");
  mcDEBUG_COMP("DivMonfermi2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompDivMonfermi2
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
#define mccompcurname  DivMonfermi2
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 33
#define nh mccDivMonfermi2_nh
#define nv mccDivMonfermi2_nv
#define Div_N mccDivMonfermi2_Div_N
#define Div_p mccDivMonfermi2_Div_p
#define Div_p2 mccDivMonfermi2_Div_p2
{   /* Declarations of DivMonfermi2=Divergence_monitor() SETTING parameters. */
char* filename = mccDivMonfermi2_filename;
MCNUM xmin = mccDivMonfermi2_xmin;
MCNUM xmax = mccDivMonfermi2_xmax;
MCNUM ymin = mccDivMonfermi2_ymin;
MCNUM ymax = mccDivMonfermi2_ymax;
MCNUM xwidth = mccDivMonfermi2_xwidth;
MCNUM yheight = mccDivMonfermi2_yheight;
MCNUM maxdiv_h = mccDivMonfermi2_maxdiv_h;
MCNUM maxdiv_v = mccDivMonfermi2_maxdiv_v;
MCNUM restore_neutron = mccDivMonfermi2_restore_neutron;
MCNUM nx = mccDivMonfermi2_nx;
MCNUM ny = mccDivMonfermi2_ny;
MCNUM nz = mccDivMonfermi2_nz;
int nowritefile = mccDivMonfermi2_nowritefile;
#line 89 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
{
    int i,j;
    double h_div, v_div;
    double v, vn;

    PROP_Z0;
    if (x>xmin && x<xmax && y>ymin && y<ymax)
    {
      /* Find length of projection onto the [nx ny nz] axis */
      vn = scalar_prod(vx, vy, vz, nx, ny, nz);
      h_div = RAD2DEG*atan2(vx,vn);
      v_div = RAD2DEG*atan2(vy,vn);
      if (h_div < maxdiv_h && h_div > -maxdiv_h &&
          v_div < maxdiv_v && v_div > -maxdiv_v)
      {
        i = floor((h_div + maxdiv_h)*nh/(2.0*maxdiv_h));
        j = floor((v_div + maxdiv_v)*nv/(2.0*maxdiv_v));
        Div_N[i][j]++;
        Div_p[i][j] += p;
        Div_p2[i][j] += p*p;
        SCATTER;
      }
    }
    if (restore_neutron) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
    }
}
#line 25583 "PSI_Focus.c"
}   /* End of DivMonfermi2=Divergence_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompDivMonfermi2:
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

  /* TRACE Component FERMITOF1 [34] */
  mccoordschange(mcposrFERMITOF1, mcrotrFERMITOF1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component FERMITOF1 (without coords transformations) */
  mcJumpTrace_FERMITOF1:
  SIG_MESSAGE("FERMITOF1 (Trace)");
  mcDEBUG_COMP("FERMITOF1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompFERMITOF1
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
#define mccompcurname  FERMITOF1
#define mccompcurtype  Monitor_nD
#define mccompcurindex 34
#define user1 mccFERMITOF1_user1
#define user2 mccFERMITOF1_user2
#define user3 mccFERMITOF1_user3
#define DEFS mccFERMITOF1_DEFS
#define Vars mccFERMITOF1_Vars
#define detector mccFERMITOF1_detector
#define offdata mccFERMITOF1_offdata
{   /* Declarations of FERMITOF1=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccFERMITOF1_xwidth;
MCNUM yheight = mccFERMITOF1_yheight;
MCNUM zdepth = mccFERMITOF1_zdepth;
MCNUM xmin = mccFERMITOF1_xmin;
MCNUM xmax = mccFERMITOF1_xmax;
MCNUM ymin = mccFERMITOF1_ymin;
MCNUM ymax = mccFERMITOF1_ymax;
MCNUM zmin = mccFERMITOF1_zmin;
MCNUM zmax = mccFERMITOF1_zmax;
MCNUM bins = mccFERMITOF1_bins;
MCNUM min = mccFERMITOF1_min;
MCNUM max = mccFERMITOF1_max;
MCNUM restore_neutron = mccFERMITOF1_restore_neutron;
MCNUM radius = mccFERMITOF1_radius;
char* options = mccFERMITOF1_options;
char* filename = mccFERMITOF1_filename;
char* geometry = mccFERMITOF1_geometry;
char* username1 = mccFERMITOF1_username1;
char* username2 = mccFERMITOF1_username2;
char* username3 = mccFERMITOF1_username3;
int nowritefile = mccFERMITOF1_nowritefile;
#line 309 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
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
#line 25891 "PSI_Focus.c"
}   /* End of FERMITOF1=Monitor_nD() SETTING parameter declarations. */
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
  mcabsorbCompFERMITOF1:
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

  /* TRACE Component SAMPLE_SLIT [35] */
  mccoordschange(mcposrSAMPLE_SLIT, mcrotrSAMPLE_SLIT,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component SAMPLE_SLIT (without coords transformations) */
  mcJumpTrace_SAMPLE_SLIT:
  SIG_MESSAGE("SAMPLE_SLIT (Trace)");
  mcDEBUG_COMP("SAMPLE_SLIT")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompSAMPLE_SLIT
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
#define mccompcurname  SAMPLE_SLIT
#define mccompcurtype  Slit
#define mccompcurindex 35
{   /* Declarations of SAMPLE_SLIT=Slit() SETTING parameters. */
MCNUM xmin = mccSAMPLE_SLIT_xmin;
MCNUM xmax = mccSAMPLE_SLIT_xmax;
MCNUM ymin = mccSAMPLE_SLIT_ymin;
MCNUM ymax = mccSAMPLE_SLIT_ymax;
MCNUM radius = mccSAMPLE_SLIT_radius;
MCNUM xwidth = mccSAMPLE_SLIT_xwidth;
MCNUM yheight = mccSAMPLE_SLIT_yheight;
#line 71 "/usr/share/mcstas/2.5/optics/Slit.comp"
{
    PROP_Z0;
    if (((radius == 0) && (x<xmin || x>xmax || y<ymin || y>ymax))
    || ((radius != 0) && (x*x + y*y > radius*radius)))
      ABSORB;
    else
        SCATTER;
}
#line 26020 "PSI_Focus.c"
}   /* End of SAMPLE_SLIT=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompSAMPLE_SLIT:
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

  /* TRACE Component FERMITOF2 [36] */
  mccoordschange(mcposrFERMITOF2, mcrotrFERMITOF2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component FERMITOF2 (without coords transformations) */
  mcJumpTrace_FERMITOF2:
  SIG_MESSAGE("FERMITOF2 (Trace)");
  mcDEBUG_COMP("FERMITOF2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompFERMITOF2
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
#define mccompcurname  FERMITOF2
#define mccompcurtype  Monitor_nD
#define mccompcurindex 36
#define user1 mccFERMITOF2_user1
#define user2 mccFERMITOF2_user2
#define user3 mccFERMITOF2_user3
#define DEFS mccFERMITOF2_DEFS
#define Vars mccFERMITOF2_Vars
#define detector mccFERMITOF2_detector
#define offdata mccFERMITOF2_offdata
{   /* Declarations of FERMITOF2=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccFERMITOF2_xwidth;
MCNUM yheight = mccFERMITOF2_yheight;
MCNUM zdepth = mccFERMITOF2_zdepth;
MCNUM xmin = mccFERMITOF2_xmin;
MCNUM xmax = mccFERMITOF2_xmax;
MCNUM ymin = mccFERMITOF2_ymin;
MCNUM ymax = mccFERMITOF2_ymax;
MCNUM zmin = mccFERMITOF2_zmin;
MCNUM zmax = mccFERMITOF2_zmax;
MCNUM bins = mccFERMITOF2_bins;
MCNUM min = mccFERMITOF2_min;
MCNUM max = mccFERMITOF2_max;
MCNUM restore_neutron = mccFERMITOF2_restore_neutron;
MCNUM radius = mccFERMITOF2_radius;
char* options = mccFERMITOF2_options;
char* filename = mccFERMITOF2_filename;
char* geometry = mccFERMITOF2_geometry;
char* username1 = mccFERMITOF2_username1;
char* username2 = mccFERMITOF2_username2;
char* username3 = mccFERMITOF2_username3;
int nowritefile = mccFERMITOF2_nowritefile;
#line 309 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
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
#line 26323 "PSI_Focus.c"
}   /* End of FERMITOF2=Monitor_nD() SETTING parameter declarations. */
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
  mcabsorbCompFERMITOF2:
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

  /* TRACE Component PSD_SAMPLE [37] */
  mccoordschange(mcposrPSD_SAMPLE, mcrotrPSD_SAMPLE,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PSD_SAMPLE (without coords transformations) */
  mcJumpTrace_PSD_SAMPLE:
  SIG_MESSAGE("PSD_SAMPLE (Trace)");
  mcDEBUG_COMP("PSD_SAMPLE")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompPSD_SAMPLE
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
#define mccompcurname  PSD_SAMPLE
#define mccompcurtype  PSD_monitor
#define mccompcurindex 37
#define nx mccPSD_SAMPLE_nx
#define ny mccPSD_SAMPLE_ny
#define PSD_N mccPSD_SAMPLE_PSD_N
#define PSD_p mccPSD_SAMPLE_PSD_p
#define PSD_p2 mccPSD_SAMPLE_PSD_p2
{   /* Declarations of PSD_SAMPLE=PSD_monitor() SETTING parameters. */
char* filename = mccPSD_SAMPLE_filename;
MCNUM xmin = mccPSD_SAMPLE_xmin;
MCNUM xmax = mccPSD_SAMPLE_xmax;
MCNUM ymin = mccPSD_SAMPLE_ymin;
MCNUM ymax = mccPSD_SAMPLE_ymax;
MCNUM xwidth = mccPSD_SAMPLE_xwidth;
MCNUM yheight = mccPSD_SAMPLE_yheight;
MCNUM restore_neutron = mccPSD_SAMPLE_restore_neutron;
int nowritefile = mccPSD_SAMPLE_nowritefile;
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
#line 26468 "PSI_Focus.c"
}   /* End of PSD_SAMPLE=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPSD_SAMPLE:
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

  /* TRACE Component DivMon_Sample [38] */
  mccoordschange(mcposrDivMon_Sample, mcrotrDivMon_Sample,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component DivMon_Sample (without coords transformations) */
  mcJumpTrace_DivMon_Sample:
  SIG_MESSAGE("DivMon_Sample (Trace)");
  mcDEBUG_COMP("DivMon_Sample")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompDivMon_Sample
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
#define mccompcurname  DivMon_Sample
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 38
#define nh mccDivMon_Sample_nh
#define nv mccDivMon_Sample_nv
#define Div_N mccDivMon_Sample_Div_N
#define Div_p mccDivMon_Sample_Div_p
#define Div_p2 mccDivMon_Sample_Div_p2
{   /* Declarations of DivMon_Sample=Divergence_monitor() SETTING parameters. */
char* filename = mccDivMon_Sample_filename;
MCNUM xmin = mccDivMon_Sample_xmin;
MCNUM xmax = mccDivMon_Sample_xmax;
MCNUM ymin = mccDivMon_Sample_ymin;
MCNUM ymax = mccDivMon_Sample_ymax;
MCNUM xwidth = mccDivMon_Sample_xwidth;
MCNUM yheight = mccDivMon_Sample_yheight;
MCNUM maxdiv_h = mccDivMon_Sample_maxdiv_h;
MCNUM maxdiv_v = mccDivMon_Sample_maxdiv_v;
MCNUM restore_neutron = mccDivMon_Sample_restore_neutron;
MCNUM nx = mccDivMon_Sample_nx;
MCNUM ny = mccDivMon_Sample_ny;
MCNUM nz = mccDivMon_Sample_nz;
int nowritefile = mccDivMon_Sample_nowritefile;
#line 89 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
{
    int i,j;
    double h_div, v_div;
    double v, vn;

    PROP_Z0;
    if (x>xmin && x<xmax && y>ymin && y<ymax)
    {
      /* Find length of projection onto the [nx ny nz] axis */
      vn = scalar_prod(vx, vy, vz, nx, ny, nz);
      h_div = RAD2DEG*atan2(vx,vn);
      v_div = RAD2DEG*atan2(vy,vn);
      if (h_div < maxdiv_h && h_div > -maxdiv_h &&
          v_div < maxdiv_v && v_div > -maxdiv_v)
      {
        i = floor((h_div + maxdiv_h)*nh/(2.0*maxdiv_h));
        j = floor((v_div + maxdiv_v)*nv/(2.0*maxdiv_v));
        Div_N[i][j]++;
        Div_p[i][j] += p;
        Div_p2[i][j] += p*p;
        SCATTER;
      }
    }
    if (restore_neutron) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
    }
}
#line 26626 "PSI_Focus.c"
}   /* End of DivMon_Sample=Divergence_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompDivMon_Sample:
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

  /* TRACE Component EMON_SAMPLE [39] */
  mccoordschange(mcposrEMON_SAMPLE, mcrotrEMON_SAMPLE,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component EMON_SAMPLE (without coords transformations) */
  mcJumpTrace_EMON_SAMPLE:
  SIG_MESSAGE("EMON_SAMPLE (Trace)");
  mcDEBUG_COMP("EMON_SAMPLE")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompEMON_SAMPLE
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
#define mccompcurname  EMON_SAMPLE
#define mccompcurtype  E_monitor
#define mccompcurindex 39
#define nE mccEMON_SAMPLE_nE
#define E_N mccEMON_SAMPLE_E_N
#define E_p mccEMON_SAMPLE_E_p
#define E_p2 mccEMON_SAMPLE_E_p2
#define S_p mccEMON_SAMPLE_S_p
#define S_pE mccEMON_SAMPLE_S_pE
#define S_pE2 mccEMON_SAMPLE_S_pE2
{   /* Declarations of EMON_SAMPLE=E_monitor() SETTING parameters. */
char* filename = mccEMON_SAMPLE_filename;
MCNUM xmin = mccEMON_SAMPLE_xmin;
MCNUM xmax = mccEMON_SAMPLE_xmax;
MCNUM ymin = mccEMON_SAMPLE_ymin;
MCNUM ymax = mccEMON_SAMPLE_ymax;
MCNUM xwidth = mccEMON_SAMPLE_xwidth;
MCNUM yheight = mccEMON_SAMPLE_yheight;
MCNUM Emin = mccEMON_SAMPLE_Emin;
MCNUM Emax = mccEMON_SAMPLE_Emax;
MCNUM restore_neutron = mccEMON_SAMPLE_restore_neutron;
int nowritefile = mccEMON_SAMPLE_nowritefile;
#line 89 "/usr/share/mcstas/2.5/monitors/E_monitor.comp"
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
#line 26782 "PSI_Focus.c"
}   /* End of EMON_SAMPLE=E_monitor() SETTING parameter declarations. */
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
  mcabsorbCompEMON_SAMPLE:
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

  /* TRACE Component a3 [40] */
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
#define mccompcurname  a3
#define mccompcurtype  Arm
#define mccompcurindex 40
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompa3:
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

  /* TRACE Component Sample [41] */
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
  if (!mcSplit_Sample) {                   /* STORE only the first time */
    if (floor(10) > 1) p /= floor(10); /* adapt weight for SPLITed neutron */
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
  } else {
    RESTORE_NEUTRON(41,
      mcnlx,
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
  mcSplit_Sample++; /* SPLIT number */
  mcScattered=0;
  mcRestore=0;
  mcNCounter[41]++;
  mcPCounter[41] += p;
  mcP2Counter[41] += p*p;
#define mccompcurname  Sample
#define mccompcurtype  V_sample
#define mccompcurindex 41
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
#line 27185 "PSI_Focus.c"
}   /* End of Sample=V_sample() SETTING parameter declarations. */
#undef VarsV
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompSample:
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

  /* TRACE Component TOF_Det [42] */
  mccoordschange(mcposrTOF_Det, mcrotrTOF_Det,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component TOF_Det (without coords transformations) */
  mcJumpTrace_TOF_Det:
  SIG_MESSAGE("TOF_Det (Trace)");
  mcDEBUG_COMP("TOF_Det")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompTOF_Det
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
#define mccompcurname  TOF_Det
#define mccompcurtype  Monitor_nD
#define mccompcurindex 42
#define user1 mccTOF_Det_user1
#define user2 mccTOF_Det_user2
#define user3 mccTOF_Det_user3
#define DEFS mccTOF_Det_DEFS
#define Vars mccTOF_Det_Vars
#define detector mccTOF_Det_detector
#define offdata mccTOF_Det_offdata
{   /* Declarations of TOF_Det=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccTOF_Det_xwidth;
MCNUM yheight = mccTOF_Det_yheight;
MCNUM zdepth = mccTOF_Det_zdepth;
MCNUM xmin = mccTOF_Det_xmin;
MCNUM xmax = mccTOF_Det_xmax;
MCNUM ymin = mccTOF_Det_ymin;
MCNUM ymax = mccTOF_Det_ymax;
MCNUM zmin = mccTOF_Det_zmin;
MCNUM zmax = mccTOF_Det_zmax;
MCNUM bins = mccTOF_Det_bins;
MCNUM min = mccTOF_Det_min;
MCNUM max = mccTOF_Det_max;
MCNUM restore_neutron = mccTOF_Det_restore_neutron;
MCNUM radius = mccTOF_Det_radius;
char* options = mccTOF_Det_options;
char* filename = mccTOF_Det_filename;
char* geometry = mccTOF_Det_geometry;
char* username1 = mccTOF_Det_username1;
char* username2 = mccTOF_Det_username2;
char* username3 = mccTOF_Det_username3;
int nowritefile = mccTOF_Det_nowritefile;
#line 309 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
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
#line 27489 "PSI_Focus.c"
}   /* End of TOF_Det=Monitor_nD() SETTING parameter declarations. */
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
  mcabsorbCompTOF_Det:
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

  /* TRACE Component FoDet [43] */
  mccoordschange(mcposrFoDet, mcrotrFoDet,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component FoDet (without coords transformations) */
  mcJumpTrace_FoDet:
  SIG_MESSAGE("FoDet (Trace)");
  mcDEBUG_COMP("FoDet")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompFoDet
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
#define mccompcurname  FoDet
#define mccompcurtype  Monitor_nD
#define mccompcurindex 43
#define user1 mccFoDet_user1
#define user2 mccFoDet_user2
#define user3 mccFoDet_user3
#define DEFS mccFoDet_DEFS
#define Vars mccFoDet_Vars
#define detector mccFoDet_detector
#define offdata mccFoDet_offdata
{   /* Declarations of FoDet=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccFoDet_xwidth;
MCNUM yheight = mccFoDet_yheight;
MCNUM zdepth = mccFoDet_zdepth;
MCNUM xmin = mccFoDet_xmin;
MCNUM xmax = mccFoDet_xmax;
MCNUM ymin = mccFoDet_ymin;
MCNUM ymax = mccFoDet_ymax;
MCNUM zmin = mccFoDet_zmin;
MCNUM zmax = mccFoDet_zmax;
MCNUM bins = mccFoDet_bins;
MCNUM min = mccFoDet_min;
MCNUM max = mccFoDet_max;
MCNUM restore_neutron = mccFoDet_restore_neutron;
MCNUM radius = mccFoDet_radius;
char* options = mccFoDet_options;
char* filename = mccFoDet_filename;
char* geometry = mccFoDet_geometry;
char* username1 = mccFoDet_username1;
char* username2 = mccFoDet_username2;
char* username3 = mccFoDet_username3;
int nowritefile = mccFoDet_nowritefile;
#line 309 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
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
#line 27799 "PSI_Focus.c"
}   /* End of FoDet=Monitor_nD() SETTING parameter declarations. */
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
  mcabsorbCompFoDet:
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

  /* TRACE Component EMON_DET [44] */
  mccoordschange(mcposrEMON_DET, mcrotrEMON_DET,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component EMON_DET (without coords transformations) */
  mcJumpTrace_EMON_DET:
  SIG_MESSAGE("EMON_DET (Trace)");
  mcDEBUG_COMP("EMON_DET")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompEMON_DET
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
#define mccompcurname  EMON_DET
#define mccompcurtype  E_monitor
#define mccompcurindex 44
#define nE mccEMON_DET_nE
#define E_N mccEMON_DET_E_N
#define E_p mccEMON_DET_E_p
#define E_p2 mccEMON_DET_E_p2
#define S_p mccEMON_DET_S_p
#define S_pE mccEMON_DET_S_pE
#define S_pE2 mccEMON_DET_S_pE2
{   /* Declarations of EMON_DET=E_monitor() SETTING parameters. */
char* filename = mccEMON_DET_filename;
MCNUM xmin = mccEMON_DET_xmin;
MCNUM xmax = mccEMON_DET_xmax;
MCNUM ymin = mccEMON_DET_ymin;
MCNUM ymax = mccEMON_DET_ymax;
MCNUM xwidth = mccEMON_DET_xwidth;
MCNUM yheight = mccEMON_DET_yheight;
MCNUM Emin = mccEMON_DET_Emin;
MCNUM Emax = mccEMON_DET_Emax;
MCNUM restore_neutron = mccEMON_DET_restore_neutron;
int nowritefile = mccEMON_DET_nowritefile;
#line 89 "/usr/share/mcstas/2.5/monitors/E_monitor.comp"
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
#line 27957 "PSI_Focus.c"
}   /* End of EMON_DET=E_monitor() SETTING parameter declarations. */
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
  mcabsorbCompEMON_DET:
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

  mcabsorbAll:
  /* SPLIT loops in reverse order */
  if (mcSplit_Sample && mcSplit_Sample < (10)) {
    goto mcJumpTrace_Sample;
  }
    else mcSplit_Sample=0;
  if (mcSplit_mono && mcSplit_mono < (10)) {
    goto mcJumpTrace_mono;
  }
    else mcSplit_mono=0;

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
#line 28074 "PSI_Focus.c"
}   /* End of a1=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'lambdaGuideExit'. */
  SIG_MESSAGE("lambdaGuideExit (Save)");
#define mccompcurname  lambdaGuideExit
#define mccompcurtype  L_monitor
#define mccompcurindex 7
#define nL mcclambdaGuideExit_nL
#define L_N mcclambdaGuideExit_L_N
#define L_p mcclambdaGuideExit_L_p
#define L_p2 mcclambdaGuideExit_L_p2
{   /* Declarations of lambdaGuideExit=L_monitor() SETTING parameters. */
char* filename = mcclambdaGuideExit_filename;
MCNUM xmin = mcclambdaGuideExit_xmin;
MCNUM xmax = mcclambdaGuideExit_xmax;
MCNUM ymin = mcclambdaGuideExit_ymin;
MCNUM ymax = mcclambdaGuideExit_ymax;
MCNUM xwidth = mcclambdaGuideExit_xwidth;
MCNUM yheight = mcclambdaGuideExit_yheight;
MCNUM Lmin = mcclambdaGuideExit_Lmin;
MCNUM Lmax = mcclambdaGuideExit_Lmax;
MCNUM restore_neutron = mcclambdaGuideExit_restore_neutron;
int nowritefile = mcclambdaGuideExit_nowritefile;
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
#line 28117 "PSI_Focus.c"
}   /* End of lambdaGuideExit=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'DivMonGuideExit'. */
  SIG_MESSAGE("DivMonGuideExit (Save)");
#define mccompcurname  DivMonGuideExit
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 8
#define nh mccDivMonGuideExit_nh
#define nv mccDivMonGuideExit_nv
#define Div_N mccDivMonGuideExit_Div_N
#define Div_p mccDivMonGuideExit_Div_p
#define Div_p2 mccDivMonGuideExit_Div_p2
{   /* Declarations of DivMonGuideExit=Divergence_monitor() SETTING parameters. */
char* filename = mccDivMonGuideExit_filename;
MCNUM xmin = mccDivMonGuideExit_xmin;
MCNUM xmax = mccDivMonGuideExit_xmax;
MCNUM ymin = mccDivMonGuideExit_ymin;
MCNUM ymax = mccDivMonGuideExit_ymax;
MCNUM xwidth = mccDivMonGuideExit_xwidth;
MCNUM yheight = mccDivMonGuideExit_yheight;
MCNUM maxdiv_h = mccDivMonGuideExit_maxdiv_h;
MCNUM maxdiv_v = mccDivMonGuideExit_maxdiv_v;
MCNUM restore_neutron = mccDivMonGuideExit_restore_neutron;
MCNUM nx = mccDivMonGuideExit_nx;
MCNUM ny = mccDivMonGuideExit_ny;
MCNUM nz = mccDivMonGuideExit_nz;
int nowritefile = mccDivMonGuideExit_nowritefile;
#line 117 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_2D(
        "Divergence monitor",
        "X divergence [deg]",
        "Y divergence [deg]",
        -maxdiv_h, maxdiv_h, -maxdiv_v, maxdiv_v,
        nh, nv,
        &Div_N[0][0],&Div_p[0][0],&Div_p2[0][0],
        filename);
    }
}
#line 28165 "PSI_Focus.c"
}   /* End of DivMonGuideExit=Divergence_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'PSDGuideExit'. */
  SIG_MESSAGE("PSDGuideExit (Save)");
#define mccompcurname  PSDGuideExit
#define mccompcurtype  PSD_monitor
#define mccompcurindex 9
#define nx mccPSDGuideExit_nx
#define ny mccPSDGuideExit_ny
#define PSD_N mccPSDGuideExit_PSD_N
#define PSD_p mccPSDGuideExit_PSD_p
#define PSD_p2 mccPSDGuideExit_PSD_p2
{   /* Declarations of PSDGuideExit=PSD_monitor() SETTING parameters. */
char* filename = mccPSDGuideExit_filename;
MCNUM xmin = mccPSDGuideExit_xmin;
MCNUM xmax = mccPSDGuideExit_xmax;
MCNUM ymin = mccPSDGuideExit_ymin;
MCNUM ymax = mccPSDGuideExit_ymax;
MCNUM xwidth = mccPSDGuideExit_xwidth;
MCNUM yheight = mccPSDGuideExit_yheight;
MCNUM restore_neutron = mccPSDGuideExit_restore_neutron;
int nowritefile = mccPSDGuideExit_nowritefile;
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
#line 28209 "PSI_Focus.c"
}   /* End of PSDGuideExit=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'DISCTOF'. */
  SIG_MESSAGE("DISCTOF (Save)");
#define mccompcurname  DISCTOF
#define mccompcurtype  Monitor_nD
#define mccompcurindex 12
#define user1 mccDISCTOF_user1
#define user2 mccDISCTOF_user2
#define user3 mccDISCTOF_user3
#define DEFS mccDISCTOF_DEFS
#define Vars mccDISCTOF_Vars
#define detector mccDISCTOF_detector
#define offdata mccDISCTOF_offdata
{   /* Declarations of DISCTOF=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccDISCTOF_xwidth;
MCNUM yheight = mccDISCTOF_yheight;
MCNUM zdepth = mccDISCTOF_zdepth;
MCNUM xmin = mccDISCTOF_xmin;
MCNUM xmax = mccDISCTOF_xmax;
MCNUM ymin = mccDISCTOF_ymin;
MCNUM ymax = mccDISCTOF_ymax;
MCNUM zmin = mccDISCTOF_zmin;
MCNUM zmax = mccDISCTOF_zmax;
MCNUM bins = mccDISCTOF_bins;
MCNUM min = mccDISCTOF_min;
MCNUM max = mccDISCTOF_max;
MCNUM restore_neutron = mccDISCTOF_restore_neutron;
MCNUM radius = mccDISCTOF_radius;
char* options = mccDISCTOF_options;
char* filename = mccDISCTOF_filename;
char* geometry = mccDISCTOF_geometry;
char* username1 = mccDISCTOF_username1;
char* username2 = mccDISCTOF_username2;
char* username3 = mccDISCTOF_username3;
int nowritefile = mccDISCTOF_nowritefile;
#line 479 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
{
  /* save results, but do not free pointers */
  detector = Monitor_nD_Save(&DEFS, &Vars);
}
#line 28259 "PSI_Focus.c"
}   /* End of DISCTOF=Monitor_nD() SETTING parameter declarations. */
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

  /* User SAVE code for component 'PSDmon1Chopper'. */
  SIG_MESSAGE("PSDmon1Chopper (Save)");
#define mccompcurname  PSDmon1Chopper
#define mccompcurtype  PSD_monitor
#define mccompcurindex 13
#define nx mccPSDmon1Chopper_nx
#define ny mccPSDmon1Chopper_ny
#define PSD_N mccPSDmon1Chopper_PSD_N
#define PSD_p mccPSDmon1Chopper_PSD_p
#define PSD_p2 mccPSDmon1Chopper_PSD_p2
{   /* Declarations of PSDmon1Chopper=PSD_monitor() SETTING parameters. */
char* filename = mccPSDmon1Chopper_filename;
MCNUM xmin = mccPSDmon1Chopper_xmin;
MCNUM xmax = mccPSDmon1Chopper_xmax;
MCNUM ymin = mccPSDmon1Chopper_ymin;
MCNUM ymax = mccPSDmon1Chopper_ymax;
MCNUM xwidth = mccPSDmon1Chopper_xwidth;
MCNUM yheight = mccPSDmon1Chopper_yheight;
MCNUM restore_neutron = mccPSDmon1Chopper_restore_neutron;
int nowritefile = mccPSDmon1Chopper_nowritefile;
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
#line 28305 "PSI_Focus.c"
}   /* End of PSDmon1Chopper=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'PSDmonMono'. */
  SIG_MESSAGE("PSDmonMono (Save)");
#define mccompcurname  PSDmonMono
#define mccompcurtype  PSD_monitor
#define mccompcurindex 20
#define nx mccPSDmonMono_nx
#define ny mccPSDmonMono_ny
#define PSD_N mccPSDmonMono_PSD_N
#define PSD_p mccPSDmonMono_PSD_p
#define PSD_p2 mccPSDmonMono_PSD_p2
{   /* Declarations of PSDmonMono=PSD_monitor() SETTING parameters. */
char* filename = mccPSDmonMono_filename;
MCNUM xmin = mccPSDmonMono_xmin;
MCNUM xmax = mccPSDmonMono_xmax;
MCNUM ymin = mccPSDmonMono_ymin;
MCNUM ymax = mccPSDmonMono_ymax;
MCNUM xwidth = mccPSDmonMono_xwidth;
MCNUM yheight = mccPSDmonMono_yheight;
MCNUM restore_neutron = mccPSDmonMono_restore_neutron;
int nowritefile = mccPSDmonMono_nowritefile;
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
#line 28349 "PSI_Focus.c"
}   /* End of PSDmonMono=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'MONOTOF'. */
  SIG_MESSAGE("MONOTOF (Save)");
#define mccompcurname  MONOTOF
#define mccompcurtype  Monitor_nD
#define mccompcurindex 21
#define user1 mccMONOTOF_user1
#define user2 mccMONOTOF_user2
#define user3 mccMONOTOF_user3
#define DEFS mccMONOTOF_DEFS
#define Vars mccMONOTOF_Vars
#define detector mccMONOTOF_detector
#define offdata mccMONOTOF_offdata
{   /* Declarations of MONOTOF=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccMONOTOF_xwidth;
MCNUM yheight = mccMONOTOF_yheight;
MCNUM zdepth = mccMONOTOF_zdepth;
MCNUM xmin = mccMONOTOF_xmin;
MCNUM xmax = mccMONOTOF_xmax;
MCNUM ymin = mccMONOTOF_ymin;
MCNUM ymax = mccMONOTOF_ymax;
MCNUM zmin = mccMONOTOF_zmin;
MCNUM zmax = mccMONOTOF_zmax;
MCNUM bins = mccMONOTOF_bins;
MCNUM min = mccMONOTOF_min;
MCNUM max = mccMONOTOF_max;
MCNUM restore_neutron = mccMONOTOF_restore_neutron;
MCNUM radius = mccMONOTOF_radius;
char* options = mccMONOTOF_options;
char* filename = mccMONOTOF_filename;
char* geometry = mccMONOTOF_geometry;
char* username1 = mccMONOTOF_username1;
char* username2 = mccMONOTOF_username2;
char* username3 = mccMONOTOF_username3;
int nowritefile = mccMONOTOF_nowritefile;
#line 479 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
{
  /* save results, but do not free pointers */
  detector = Monitor_nD_Save(&DEFS, &Vars);
}
#line 28399 "PSI_Focus.c"
}   /* End of MONOTOF=Monitor_nD() SETTING parameter declarations. */
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

  /* User SAVE code for component 'DivMonMono'. */
  SIG_MESSAGE("DivMonMono (Save)");
#define mccompcurname  DivMonMono
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 22
#define nh mccDivMonMono_nh
#define nv mccDivMonMono_nv
#define Div_N mccDivMonMono_Div_N
#define Div_p mccDivMonMono_Div_p
#define Div_p2 mccDivMonMono_Div_p2
{   /* Declarations of DivMonMono=Divergence_monitor() SETTING parameters. */
char* filename = mccDivMonMono_filename;
MCNUM xmin = mccDivMonMono_xmin;
MCNUM xmax = mccDivMonMono_xmax;
MCNUM ymin = mccDivMonMono_ymin;
MCNUM ymax = mccDivMonMono_ymax;
MCNUM xwidth = mccDivMonMono_xwidth;
MCNUM yheight = mccDivMonMono_yheight;
MCNUM maxdiv_h = mccDivMonMono_maxdiv_h;
MCNUM maxdiv_v = mccDivMonMono_maxdiv_v;
MCNUM restore_neutron = mccDivMonMono_restore_neutron;
MCNUM nx = mccDivMonMono_nx;
MCNUM ny = mccDivMonMono_ny;
MCNUM nz = mccDivMonMono_nz;
int nowritefile = mccDivMonMono_nowritefile;
#line 117 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_2D(
        "Divergence monitor",
        "X divergence [deg]",
        "Y divergence [deg]",
        -maxdiv_h, maxdiv_h, -maxdiv_v, maxdiv_v,
        nh, nv,
        &Div_N[0][0],&Div_p[0][0],&Div_p2[0][0],
        filename);
    }
}
#line 28450 "PSI_Focus.c"
}   /* End of DivMonMono=Divergence_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'FERMITOF_before'. */
  SIG_MESSAGE("FERMITOF_before (Save)");
#define mccompcurname  FERMITOF_before
#define mccompcurtype  Monitor_nD
#define mccompcurindex 26
#define user1 mccFERMITOF_before_user1
#define user2 mccFERMITOF_before_user2
#define user3 mccFERMITOF_before_user3
#define DEFS mccFERMITOF_before_DEFS
#define Vars mccFERMITOF_before_Vars
#define detector mccFERMITOF_before_detector
#define offdata mccFERMITOF_before_offdata
{   /* Declarations of FERMITOF_before=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccFERMITOF_before_xwidth;
MCNUM yheight = mccFERMITOF_before_yheight;
MCNUM zdepth = mccFERMITOF_before_zdepth;
MCNUM xmin = mccFERMITOF_before_xmin;
MCNUM xmax = mccFERMITOF_before_xmax;
MCNUM ymin = mccFERMITOF_before_ymin;
MCNUM ymax = mccFERMITOF_before_ymax;
MCNUM zmin = mccFERMITOF_before_zmin;
MCNUM zmax = mccFERMITOF_before_zmax;
MCNUM bins = mccFERMITOF_before_bins;
MCNUM min = mccFERMITOF_before_min;
MCNUM max = mccFERMITOF_before_max;
MCNUM restore_neutron = mccFERMITOF_before_restore_neutron;
MCNUM radius = mccFERMITOF_before_radius;
char* options = mccFERMITOF_before_options;
char* filename = mccFERMITOF_before_filename;
char* geometry = mccFERMITOF_before_geometry;
char* username1 = mccFERMITOF_before_username1;
char* username2 = mccFERMITOF_before_username2;
char* username3 = mccFERMITOF_before_username3;
int nowritefile = mccFERMITOF_before_nowritefile;
#line 479 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
{
  /* save results, but do not free pointers */
  detector = Monitor_nD_Save(&DEFS, &Vars);
}
#line 28500 "PSI_Focus.c"
}   /* End of FERMITOF_before=Monitor_nD() SETTING parameter declarations. */
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

  /* User SAVE code for component 'lambdaFermi'. */
  SIG_MESSAGE("lambdaFermi (Save)");
#define mccompcurname  lambdaFermi
#define mccompcurtype  L_monitor
#define mccompcurindex 27
#define nL mcclambdaFermi_nL
#define L_N mcclambdaFermi_L_N
#define L_p mcclambdaFermi_L_p
#define L_p2 mcclambdaFermi_L_p2
{   /* Declarations of lambdaFermi=L_monitor() SETTING parameters. */
char* filename = mcclambdaFermi_filename;
MCNUM xmin = mcclambdaFermi_xmin;
MCNUM xmax = mcclambdaFermi_xmax;
MCNUM ymin = mcclambdaFermi_ymin;
MCNUM ymax = mcclambdaFermi_ymax;
MCNUM xwidth = mcclambdaFermi_xwidth;
MCNUM yheight = mcclambdaFermi_yheight;
MCNUM Lmin = mcclambdaFermi_Lmin;
MCNUM Lmax = mcclambdaFermi_Lmax;
MCNUM restore_neutron = mcclambdaFermi_restore_neutron;
int nowritefile = mcclambdaFermi_nowritefile;
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
#line 28546 "PSI_Focus.c"
}   /* End of lambdaFermi=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'EMON_Fermi'. */
  SIG_MESSAGE("EMON_Fermi (Save)");
#define mccompcurname  EMON_Fermi
#define mccompcurtype  E_monitor
#define mccompcurindex 28
#define nE mccEMON_Fermi_nE
#define E_N mccEMON_Fermi_E_N
#define E_p mccEMON_Fermi_E_p
#define E_p2 mccEMON_Fermi_E_p2
#define S_p mccEMON_Fermi_S_p
#define S_pE mccEMON_Fermi_S_pE
#define S_pE2 mccEMON_Fermi_S_pE2
{   /* Declarations of EMON_Fermi=E_monitor() SETTING parameters. */
char* filename = mccEMON_Fermi_filename;
MCNUM xmin = mccEMON_Fermi_xmin;
MCNUM xmax = mccEMON_Fermi_xmax;
MCNUM ymin = mccEMON_Fermi_ymin;
MCNUM ymax = mccEMON_Fermi_ymax;
MCNUM xwidth = mccEMON_Fermi_xwidth;
MCNUM yheight = mccEMON_Fermi_yheight;
MCNUM Emin = mccEMON_Fermi_Emin;
MCNUM Emax = mccEMON_Fermi_Emax;
MCNUM restore_neutron = mccEMON_Fermi_restore_neutron;
int nowritefile = mccEMON_Fermi_nowritefile;
#line 117 "/usr/share/mcstas/2.5/monitors/E_monitor.comp"
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
#line 28594 "PSI_Focus.c"
}   /* End of EMON_Fermi=E_monitor() SETTING parameter declarations. */
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

  /* User SAVE code for component 'DivMonfermi1'. */
  SIG_MESSAGE("DivMonfermi1 (Save)");
#define mccompcurname  DivMonfermi1
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 29
#define nh mccDivMonfermi1_nh
#define nv mccDivMonfermi1_nv
#define Div_N mccDivMonfermi1_Div_N
#define Div_p mccDivMonfermi1_Div_p
#define Div_p2 mccDivMonfermi1_Div_p2
{   /* Declarations of DivMonfermi1=Divergence_monitor() SETTING parameters. */
char* filename = mccDivMonfermi1_filename;
MCNUM xmin = mccDivMonfermi1_xmin;
MCNUM xmax = mccDivMonfermi1_xmax;
MCNUM ymin = mccDivMonfermi1_ymin;
MCNUM ymax = mccDivMonfermi1_ymax;
MCNUM xwidth = mccDivMonfermi1_xwidth;
MCNUM yheight = mccDivMonfermi1_yheight;
MCNUM maxdiv_h = mccDivMonfermi1_maxdiv_h;
MCNUM maxdiv_v = mccDivMonfermi1_maxdiv_v;
MCNUM restore_neutron = mccDivMonfermi1_restore_neutron;
MCNUM nx = mccDivMonfermi1_nx;
MCNUM ny = mccDivMonfermi1_ny;
MCNUM nz = mccDivMonfermi1_nz;
int nowritefile = mccDivMonfermi1_nowritefile;
#line 117 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_2D(
        "Divergence monitor",
        "X divergence [deg]",
        "Y divergence [deg]",
        -maxdiv_h, maxdiv_h, -maxdiv_v, maxdiv_v,
        nh, nv,
        &Div_N[0][0],&Div_p[0][0],&Div_p2[0][0],
        filename);
    }
}
#line 28645 "PSI_Focus.c"
}   /* End of DivMonfermi1=Divergence_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'PSD_Fermi1'. */
  SIG_MESSAGE("PSD_Fermi1 (Save)");
#define mccompcurname  PSD_Fermi1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 30
#define nx mccPSD_Fermi1_nx
#define ny mccPSD_Fermi1_ny
#define PSD_N mccPSD_Fermi1_PSD_N
#define PSD_p mccPSD_Fermi1_PSD_p
#define PSD_p2 mccPSD_Fermi1_PSD_p2
{   /* Declarations of PSD_Fermi1=PSD_monitor() SETTING parameters. */
char* filename = mccPSD_Fermi1_filename;
MCNUM xmin = mccPSD_Fermi1_xmin;
MCNUM xmax = mccPSD_Fermi1_xmax;
MCNUM ymin = mccPSD_Fermi1_ymin;
MCNUM ymax = mccPSD_Fermi1_ymax;
MCNUM xwidth = mccPSD_Fermi1_xwidth;
MCNUM yheight = mccPSD_Fermi1_yheight;
MCNUM restore_neutron = mccPSD_Fermi1_restore_neutron;
int nowritefile = mccPSD_Fermi1_nowritefile;
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
#line 28689 "PSI_Focus.c"
}   /* End of PSD_Fermi1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'FoChopper'. */
  SIG_MESSAGE("FoChopper (Save)");
#define mccompcurname  FoChopper
#define mccompcurtype  FermiChopper
#define mccompcurindex 31
#define FCVars mccFoChopper_FCVars
{   /* Declarations of FoChopper=FermiChopper() SETTING parameters. */
MCNUM phase = mccFoChopper_phase;
MCNUM radius = mccFoChopper_radius;
MCNUM nu = mccFoChopper_nu;
MCNUM w = mccFoChopper_w;
MCNUM nslit = mccFoChopper_nslit;
MCNUM R0 = mccFoChopper_R0;
MCNUM Qc = mccFoChopper_Qc;
MCNUM alpha = mccFoChopper_alpha;
MCNUM m = mccFoChopper_m;
MCNUM W = mccFoChopper_W;
MCNUM length = mccFoChopper_length;
MCNUM eff = mccFoChopper_eff;
MCNUM zero_time = mccFoChopper_zero_time;
MCNUM xwidth = mccFoChopper_xwidth;
MCNUM verbose = mccFoChopper_verbose;
MCNUM yheight = mccFoChopper_yheight;
MCNUM curvature = mccFoChopper_curvature;
MCNUM delay = mccFoChopper_delay;
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
#line 28802 "PSI_Focus.c"
}   /* End of FoChopper=FermiChopper() SETTING parameter declarations. */
#undef FCVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'PSD_Fermi2'. */
  SIG_MESSAGE("PSD_Fermi2 (Save)");
#define mccompcurname  PSD_Fermi2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 32
#define nx mccPSD_Fermi2_nx
#define ny mccPSD_Fermi2_ny
#define PSD_N mccPSD_Fermi2_PSD_N
#define PSD_p mccPSD_Fermi2_PSD_p
#define PSD_p2 mccPSD_Fermi2_PSD_p2
{   /* Declarations of PSD_Fermi2=PSD_monitor() SETTING parameters. */
char* filename = mccPSD_Fermi2_filename;
MCNUM xmin = mccPSD_Fermi2_xmin;
MCNUM xmax = mccPSD_Fermi2_xmax;
MCNUM ymin = mccPSD_Fermi2_ymin;
MCNUM ymax = mccPSD_Fermi2_ymax;
MCNUM xwidth = mccPSD_Fermi2_xwidth;
MCNUM yheight = mccPSD_Fermi2_yheight;
MCNUM restore_neutron = mccPSD_Fermi2_restore_neutron;
int nowritefile = mccPSD_Fermi2_nowritefile;
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
#line 28842 "PSI_Focus.c"
}   /* End of PSD_Fermi2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'DivMonfermi2'. */
  SIG_MESSAGE("DivMonfermi2 (Save)");
#define mccompcurname  DivMonfermi2
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 33
#define nh mccDivMonfermi2_nh
#define nv mccDivMonfermi2_nv
#define Div_N mccDivMonfermi2_Div_N
#define Div_p mccDivMonfermi2_Div_p
#define Div_p2 mccDivMonfermi2_Div_p2
{   /* Declarations of DivMonfermi2=Divergence_monitor() SETTING parameters. */
char* filename = mccDivMonfermi2_filename;
MCNUM xmin = mccDivMonfermi2_xmin;
MCNUM xmax = mccDivMonfermi2_xmax;
MCNUM ymin = mccDivMonfermi2_ymin;
MCNUM ymax = mccDivMonfermi2_ymax;
MCNUM xwidth = mccDivMonfermi2_xwidth;
MCNUM yheight = mccDivMonfermi2_yheight;
MCNUM maxdiv_h = mccDivMonfermi2_maxdiv_h;
MCNUM maxdiv_v = mccDivMonfermi2_maxdiv_v;
MCNUM restore_neutron = mccDivMonfermi2_restore_neutron;
MCNUM nx = mccDivMonfermi2_nx;
MCNUM ny = mccDivMonfermi2_ny;
MCNUM nz = mccDivMonfermi2_nz;
int nowritefile = mccDivMonfermi2_nowritefile;
#line 117 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_2D(
        "Divergence monitor",
        "X divergence [deg]",
        "Y divergence [deg]",
        -maxdiv_h, maxdiv_h, -maxdiv_v, maxdiv_v,
        nh, nv,
        &Div_N[0][0],&Div_p[0][0],&Div_p2[0][0],
        filename);
    }
}
#line 28891 "PSI_Focus.c"
}   /* End of DivMonfermi2=Divergence_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'FERMITOF1'. */
  SIG_MESSAGE("FERMITOF1 (Save)");
#define mccompcurname  FERMITOF1
#define mccompcurtype  Monitor_nD
#define mccompcurindex 34
#define user1 mccFERMITOF1_user1
#define user2 mccFERMITOF1_user2
#define user3 mccFERMITOF1_user3
#define DEFS mccFERMITOF1_DEFS
#define Vars mccFERMITOF1_Vars
#define detector mccFERMITOF1_detector
#define offdata mccFERMITOF1_offdata
{   /* Declarations of FERMITOF1=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccFERMITOF1_xwidth;
MCNUM yheight = mccFERMITOF1_yheight;
MCNUM zdepth = mccFERMITOF1_zdepth;
MCNUM xmin = mccFERMITOF1_xmin;
MCNUM xmax = mccFERMITOF1_xmax;
MCNUM ymin = mccFERMITOF1_ymin;
MCNUM ymax = mccFERMITOF1_ymax;
MCNUM zmin = mccFERMITOF1_zmin;
MCNUM zmax = mccFERMITOF1_zmax;
MCNUM bins = mccFERMITOF1_bins;
MCNUM min = mccFERMITOF1_min;
MCNUM max = mccFERMITOF1_max;
MCNUM restore_neutron = mccFERMITOF1_restore_neutron;
MCNUM radius = mccFERMITOF1_radius;
char* options = mccFERMITOF1_options;
char* filename = mccFERMITOF1_filename;
char* geometry = mccFERMITOF1_geometry;
char* username1 = mccFERMITOF1_username1;
char* username2 = mccFERMITOF1_username2;
char* username3 = mccFERMITOF1_username3;
int nowritefile = mccFERMITOF1_nowritefile;
#line 479 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
{
  /* save results, but do not free pointers */
  detector = Monitor_nD_Save(&DEFS, &Vars);
}
#line 28941 "PSI_Focus.c"
}   /* End of FERMITOF1=Monitor_nD() SETTING parameter declarations. */
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

  /* User SAVE code for component 'FERMITOF2'. */
  SIG_MESSAGE("FERMITOF2 (Save)");
#define mccompcurname  FERMITOF2
#define mccompcurtype  Monitor_nD
#define mccompcurindex 36
#define user1 mccFERMITOF2_user1
#define user2 mccFERMITOF2_user2
#define user3 mccFERMITOF2_user3
#define DEFS mccFERMITOF2_DEFS
#define Vars mccFERMITOF2_Vars
#define detector mccFERMITOF2_detector
#define offdata mccFERMITOF2_offdata
{   /* Declarations of FERMITOF2=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccFERMITOF2_xwidth;
MCNUM yheight = mccFERMITOF2_yheight;
MCNUM zdepth = mccFERMITOF2_zdepth;
MCNUM xmin = mccFERMITOF2_xmin;
MCNUM xmax = mccFERMITOF2_xmax;
MCNUM ymin = mccFERMITOF2_ymin;
MCNUM ymax = mccFERMITOF2_ymax;
MCNUM zmin = mccFERMITOF2_zmin;
MCNUM zmax = mccFERMITOF2_zmax;
MCNUM bins = mccFERMITOF2_bins;
MCNUM min = mccFERMITOF2_min;
MCNUM max = mccFERMITOF2_max;
MCNUM restore_neutron = mccFERMITOF2_restore_neutron;
MCNUM radius = mccFERMITOF2_radius;
char* options = mccFERMITOF2_options;
char* filename = mccFERMITOF2_filename;
char* geometry = mccFERMITOF2_geometry;
char* username1 = mccFERMITOF2_username1;
char* username2 = mccFERMITOF2_username2;
char* username3 = mccFERMITOF2_username3;
int nowritefile = mccFERMITOF2_nowritefile;
#line 479 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
{
  /* save results, but do not free pointers */
  detector = Monitor_nD_Save(&DEFS, &Vars);
}
#line 28993 "PSI_Focus.c"
}   /* End of FERMITOF2=Monitor_nD() SETTING parameter declarations. */
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

  /* User SAVE code for component 'PSD_SAMPLE'. */
  SIG_MESSAGE("PSD_SAMPLE (Save)");
#define mccompcurname  PSD_SAMPLE
#define mccompcurtype  PSD_monitor
#define mccompcurindex 37
#define nx mccPSD_SAMPLE_nx
#define ny mccPSD_SAMPLE_ny
#define PSD_N mccPSD_SAMPLE_PSD_N
#define PSD_p mccPSD_SAMPLE_PSD_p
#define PSD_p2 mccPSD_SAMPLE_PSD_p2
{   /* Declarations of PSD_SAMPLE=PSD_monitor() SETTING parameters. */
char* filename = mccPSD_SAMPLE_filename;
MCNUM xmin = mccPSD_SAMPLE_xmin;
MCNUM xmax = mccPSD_SAMPLE_xmax;
MCNUM ymin = mccPSD_SAMPLE_ymin;
MCNUM ymax = mccPSD_SAMPLE_ymax;
MCNUM xwidth = mccPSD_SAMPLE_xwidth;
MCNUM yheight = mccPSD_SAMPLE_yheight;
MCNUM restore_neutron = mccPSD_SAMPLE_restore_neutron;
int nowritefile = mccPSD_SAMPLE_nowritefile;
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
#line 29039 "PSI_Focus.c"
}   /* End of PSD_SAMPLE=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'DivMon_Sample'. */
  SIG_MESSAGE("DivMon_Sample (Save)");
#define mccompcurname  DivMon_Sample
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 38
#define nh mccDivMon_Sample_nh
#define nv mccDivMon_Sample_nv
#define Div_N mccDivMon_Sample_Div_N
#define Div_p mccDivMon_Sample_Div_p
#define Div_p2 mccDivMon_Sample_Div_p2
{   /* Declarations of DivMon_Sample=Divergence_monitor() SETTING parameters. */
char* filename = mccDivMon_Sample_filename;
MCNUM xmin = mccDivMon_Sample_xmin;
MCNUM xmax = mccDivMon_Sample_xmax;
MCNUM ymin = mccDivMon_Sample_ymin;
MCNUM ymax = mccDivMon_Sample_ymax;
MCNUM xwidth = mccDivMon_Sample_xwidth;
MCNUM yheight = mccDivMon_Sample_yheight;
MCNUM maxdiv_h = mccDivMon_Sample_maxdiv_h;
MCNUM maxdiv_v = mccDivMon_Sample_maxdiv_v;
MCNUM restore_neutron = mccDivMon_Sample_restore_neutron;
MCNUM nx = mccDivMon_Sample_nx;
MCNUM ny = mccDivMon_Sample_ny;
MCNUM nz = mccDivMon_Sample_nz;
int nowritefile = mccDivMon_Sample_nowritefile;
#line 117 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_2D(
        "Divergence monitor",
        "X divergence [deg]",
        "Y divergence [deg]",
        -maxdiv_h, maxdiv_h, -maxdiv_v, maxdiv_v,
        nh, nv,
        &Div_N[0][0],&Div_p[0][0],&Div_p2[0][0],
        filename);
    }
}
#line 29088 "PSI_Focus.c"
}   /* End of DivMon_Sample=Divergence_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'EMON_SAMPLE'. */
  SIG_MESSAGE("EMON_SAMPLE (Save)");
#define mccompcurname  EMON_SAMPLE
#define mccompcurtype  E_monitor
#define mccompcurindex 39
#define nE mccEMON_SAMPLE_nE
#define E_N mccEMON_SAMPLE_E_N
#define E_p mccEMON_SAMPLE_E_p
#define E_p2 mccEMON_SAMPLE_E_p2
#define S_p mccEMON_SAMPLE_S_p
#define S_pE mccEMON_SAMPLE_S_pE
#define S_pE2 mccEMON_SAMPLE_S_pE2
{   /* Declarations of EMON_SAMPLE=E_monitor() SETTING parameters. */
char* filename = mccEMON_SAMPLE_filename;
MCNUM xmin = mccEMON_SAMPLE_xmin;
MCNUM xmax = mccEMON_SAMPLE_xmax;
MCNUM ymin = mccEMON_SAMPLE_ymin;
MCNUM ymax = mccEMON_SAMPLE_ymax;
MCNUM xwidth = mccEMON_SAMPLE_xwidth;
MCNUM yheight = mccEMON_SAMPLE_yheight;
MCNUM Emin = mccEMON_SAMPLE_Emin;
MCNUM Emax = mccEMON_SAMPLE_Emax;
MCNUM restore_neutron = mccEMON_SAMPLE_restore_neutron;
int nowritefile = mccEMON_SAMPLE_nowritefile;
#line 117 "/usr/share/mcstas/2.5/monitors/E_monitor.comp"
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
#line 29137 "PSI_Focus.c"
}   /* End of EMON_SAMPLE=E_monitor() SETTING parameter declarations. */
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

  /* User SAVE code for component 'TOF_Det'. */
  SIG_MESSAGE("TOF_Det (Save)");
#define mccompcurname  TOF_Det
#define mccompcurtype  Monitor_nD
#define mccompcurindex 42
#define user1 mccTOF_Det_user1
#define user2 mccTOF_Det_user2
#define user3 mccTOF_Det_user3
#define DEFS mccTOF_Det_DEFS
#define Vars mccTOF_Det_Vars
#define detector mccTOF_Det_detector
#define offdata mccTOF_Det_offdata
{   /* Declarations of TOF_Det=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccTOF_Det_xwidth;
MCNUM yheight = mccTOF_Det_yheight;
MCNUM zdepth = mccTOF_Det_zdepth;
MCNUM xmin = mccTOF_Det_xmin;
MCNUM xmax = mccTOF_Det_xmax;
MCNUM ymin = mccTOF_Det_ymin;
MCNUM ymax = mccTOF_Det_ymax;
MCNUM zmin = mccTOF_Det_zmin;
MCNUM zmax = mccTOF_Det_zmax;
MCNUM bins = mccTOF_Det_bins;
MCNUM min = mccTOF_Det_min;
MCNUM max = mccTOF_Det_max;
MCNUM restore_neutron = mccTOF_Det_restore_neutron;
MCNUM radius = mccTOF_Det_radius;
char* options = mccTOF_Det_options;
char* filename = mccTOF_Det_filename;
char* geometry = mccTOF_Det_geometry;
char* username1 = mccTOF_Det_username1;
char* username2 = mccTOF_Det_username2;
char* username3 = mccTOF_Det_username3;
int nowritefile = mccTOF_Det_nowritefile;
#line 479 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
{
  /* save results, but do not free pointers */
  detector = Monitor_nD_Save(&DEFS, &Vars);
}
#line 29189 "PSI_Focus.c"
}   /* End of TOF_Det=Monitor_nD() SETTING parameter declarations. */
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

  /* User SAVE code for component 'FoDet'. */
  SIG_MESSAGE("FoDet (Save)");
#define mccompcurname  FoDet
#define mccompcurtype  Monitor_nD
#define mccompcurindex 43
#define user1 mccFoDet_user1
#define user2 mccFoDet_user2
#define user3 mccFoDet_user3
#define DEFS mccFoDet_DEFS
#define Vars mccFoDet_Vars
#define detector mccFoDet_detector
#define offdata mccFoDet_offdata
{   /* Declarations of FoDet=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccFoDet_xwidth;
MCNUM yheight = mccFoDet_yheight;
MCNUM zdepth = mccFoDet_zdepth;
MCNUM xmin = mccFoDet_xmin;
MCNUM xmax = mccFoDet_xmax;
MCNUM ymin = mccFoDet_ymin;
MCNUM ymax = mccFoDet_ymax;
MCNUM zmin = mccFoDet_zmin;
MCNUM zmax = mccFoDet_zmax;
MCNUM bins = mccFoDet_bins;
MCNUM min = mccFoDet_min;
MCNUM max = mccFoDet_max;
MCNUM restore_neutron = mccFoDet_restore_neutron;
MCNUM radius = mccFoDet_radius;
char* options = mccFoDet_options;
char* filename = mccFoDet_filename;
char* geometry = mccFoDet_geometry;
char* username1 = mccFoDet_username1;
char* username2 = mccFoDet_username2;
char* username3 = mccFoDet_username3;
int nowritefile = mccFoDet_nowritefile;
#line 479 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
{
  /* save results, but do not free pointers */
  detector = Monitor_nD_Save(&DEFS, &Vars);
}
#line 29241 "PSI_Focus.c"
}   /* End of FoDet=Monitor_nD() SETTING parameter declarations. */
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

  /* User SAVE code for component 'EMON_DET'. */
  SIG_MESSAGE("EMON_DET (Save)");
#define mccompcurname  EMON_DET
#define mccompcurtype  E_monitor
#define mccompcurindex 44
#define nE mccEMON_DET_nE
#define E_N mccEMON_DET_E_N
#define E_p mccEMON_DET_E_p
#define E_p2 mccEMON_DET_E_p2
#define S_p mccEMON_DET_S_p
#define S_pE mccEMON_DET_S_pE
#define S_pE2 mccEMON_DET_S_pE2
{   /* Declarations of EMON_DET=E_monitor() SETTING parameters. */
char* filename = mccEMON_DET_filename;
MCNUM xmin = mccEMON_DET_xmin;
MCNUM xmax = mccEMON_DET_xmax;
MCNUM ymin = mccEMON_DET_ymin;
MCNUM ymax = mccEMON_DET_ymax;
MCNUM xwidth = mccEMON_DET_xwidth;
MCNUM yheight = mccEMON_DET_yheight;
MCNUM Emin = mccEMON_DET_Emin;
MCNUM Emax = mccEMON_DET_Emax;
MCNUM restore_neutron = mccEMON_DET_restore_neutron;
int nowritefile = mccEMON_DET_nowritefile;
#line 117 "/usr/share/mcstas/2.5/monitors/E_monitor.comp"
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
#line 29292 "PSI_Focus.c"
}   /* End of EMON_DET=E_monitor() SETTING parameter declarations. */
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
#line 29339 "PSI_Focus.c"
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
  /* User FINALLY code for component 'csource'. */
  SIG_MESSAGE("csource (Finally)");
#define mccompcurname  csource
#define mccompcurtype  Source_gen
#define mccompcurindex 2
#define p_in mcccsource_p_in
#define lambda1 mcccsource_lambda1
#define lambda2 mcccsource_lambda2
#define lambda3 mcccsource_lambda3
#define pTable mcccsource_pTable
#define pTable_x mcccsource_pTable_x
#define pTable_y mcccsource_pTable_y
#define pTable_xmin mcccsource_pTable_xmin
#define pTable_xmax mcccsource_pTable_xmax
#define pTable_xsum mcccsource_pTable_xsum
#define pTable_ymin mcccsource_pTable_ymin
#define pTable_ymax mcccsource_pTable_ymax
#define pTable_ysum mcccsource_pTable_ysum
#define pTable_dxmin mcccsource_pTable_dxmin
#define pTable_dxmax mcccsource_pTable_dxmax
#define pTable_dymin mcccsource_pTable_dymin
#define pTable_dymax mcccsource_pTable_dymax
{   /* Declarations of csource=Source_gen() SETTING parameters. */
char* flux_file = mcccsource_flux_file;
char* xdiv_file = mcccsource_xdiv_file;
char* ydiv_file = mcccsource_ydiv_file;
MCNUM radius = mcccsource_radius;
MCNUM dist = mcccsource_dist;
MCNUM focus_xw = mcccsource_focus_xw;
MCNUM focus_yh = mcccsource_focus_yh;
MCNUM focus_aw = mcccsource_focus_aw;
MCNUM focus_ah = mcccsource_focus_ah;
MCNUM E0 = mcccsource_E0;
MCNUM dE = mcccsource_dE;
MCNUM lambda0 = mcccsource_lambda0;
MCNUM dlambda = mcccsource_dlambda;
MCNUM I1 = mcccsource_I1;
MCNUM yheight = mcccsource_yheight;
MCNUM xwidth = mcccsource_xwidth;
MCNUM verbose = mcccsource_verbose;
MCNUM T1 = mcccsource_T1;
MCNUM flux_file_perAA = mcccsource_flux_file_perAA;
MCNUM flux_file_log = mcccsource_flux_file_log;
MCNUM Lmin = mcccsource_Lmin;
MCNUM Lmax = mcccsource_Lmax;
MCNUM Emin = mcccsource_Emin;
MCNUM Emax = mcccsource_Emax;
MCNUM T2 = mcccsource_T2;
MCNUM I2 = mcccsource_I2;
MCNUM T3 = mcccsource_T3;
MCNUM I3 = mcccsource_I3;
MCNUM zdepth = mcccsource_zdepth;
int target_index = mcccsource_target_index;
#line 571 "/usr/share/mcstas/2.5/sources/Source_gen.comp"
{
  Table_Free(&pTable);
  Table_Free(&pTable_x);
  Table_Free(&pTable_y);
}
#line 29410 "PSI_Focus.c"
}   /* End of csource=Source_gen() SETTING parameter declarations. */
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

    if (!mcNCounter[2]) fprintf(stderr, "Warning: No neutron could reach Component[2] csource\n");
    if (mcAbsorbProp[2]) fprintf(stderr, "Warning: %g events were removed in Component[2] csource=Source_gen()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[2]);
    if (!mcNCounter[3]) fprintf(stderr, "Warning: No neutron could reach Component[3] guide1\n");
    if (mcAbsorbProp[3]) fprintf(stderr, "Warning: %g events were removed in Component[3] guide1=Guide()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[3]);
    if (!mcNCounter[4]) fprintf(stderr, "Warning: No neutron could reach Component[4] guide2\n");
    if (mcAbsorbProp[4]) fprintf(stderr, "Warning: %g events were removed in Component[4] guide2=Bender()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[4]);
    if (!mcNCounter[5]) fprintf(stderr, "Warning: No neutron could reach Component[5] bunker\n");
    if (mcAbsorbProp[5]) fprintf(stderr, "Warning: %g events were removed in Component[5] bunker=Guide()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[5]);
    if (!mcNCounter[6]) fprintf(stderr, "Warning: No neutron could reach Component[6] guide3\n");
    if (mcAbsorbProp[6]) fprintf(stderr, "Warning: %g events were removed in Component[6] guide3=Guide()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[6]);
    if (!mcNCounter[7]) fprintf(stderr, "Warning: No neutron could reach Component[7] lambdaGuideExit\n");
    if (mcAbsorbProp[7]) fprintf(stderr, "Warning: %g events were removed in Component[7] lambdaGuideExit=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[7]);
    if (!mcNCounter[8]) fprintf(stderr, "Warning: No neutron could reach Component[8] DivMonGuideExit\n");
    if (mcAbsorbProp[8]) fprintf(stderr, "Warning: %g events were removed in Component[8] DivMonGuideExit=Divergence_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[8]);
    if (!mcNCounter[9]) fprintf(stderr, "Warning: No neutron could reach Component[9] PSDGuideExit\n");
    if (mcAbsorbProp[9]) fprintf(stderr, "Warning: %g events were removed in Component[9] PSDGuideExit=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[9]);
    if (!mcNCounter[10]) fprintf(stderr, "Warning: No neutron could reach Component[10] FOCUSguide\n");
    if (mcAbsorbProp[10]) fprintf(stderr, "Warning: %g events were removed in Component[10] FOCUSguide=Guide_channeled()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[10]);
    if (!mcNCounter[11]) fprintf(stderr, "Warning: No neutron could reach Component[11] FirstChopper\n");
    if (mcAbsorbProp[11]) fprintf(stderr, "Warning: %g events were removed in Component[11] FirstChopper=DiskChopper()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[11]);
  /* User FINALLY code for component 'DISCTOF'. */
  SIG_MESSAGE("DISCTOF (Finally)");
#define mccompcurname  DISCTOF
#define mccompcurtype  Monitor_nD
#define mccompcurindex 12
#define user1 mccDISCTOF_user1
#define user2 mccDISCTOF_user2
#define user3 mccDISCTOF_user3
#define DEFS mccDISCTOF_DEFS
#define Vars mccDISCTOF_Vars
#define detector mccDISCTOF_detector
#define offdata mccDISCTOF_offdata
{   /* Declarations of DISCTOF=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccDISCTOF_xwidth;
MCNUM yheight = mccDISCTOF_yheight;
MCNUM zdepth = mccDISCTOF_zdepth;
MCNUM xmin = mccDISCTOF_xmin;
MCNUM xmax = mccDISCTOF_xmax;
MCNUM ymin = mccDISCTOF_ymin;
MCNUM ymax = mccDISCTOF_ymax;
MCNUM zmin = mccDISCTOF_zmin;
MCNUM zmax = mccDISCTOF_zmax;
MCNUM bins = mccDISCTOF_bins;
MCNUM min = mccDISCTOF_min;
MCNUM max = mccDISCTOF_max;
MCNUM restore_neutron = mccDISCTOF_restore_neutron;
MCNUM radius = mccDISCTOF_radius;
char* options = mccDISCTOF_options;
char* filename = mccDISCTOF_filename;
char* geometry = mccDISCTOF_geometry;
char* username1 = mccDISCTOF_username1;
char* username2 = mccDISCTOF_username2;
char* username3 = mccDISCTOF_username3;
int nowritefile = mccDISCTOF_nowritefile;
#line 485 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
{
  /* free pointers */
  if (!nowritefile) {
    Monitor_nD_Finally(&DEFS, &Vars);
  }
}
#line 29494 "PSI_Focus.c"
}   /* End of DISCTOF=Monitor_nD() SETTING parameter declarations. */
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

    if (!mcNCounter[12]) fprintf(stderr, "Warning: No neutron could reach Component[12] DISCTOF\n");
    if (mcAbsorbProp[12]) fprintf(stderr, "Warning: %g events were removed in Component[12] DISCTOF=Monitor_nD()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[12]);
    if (!mcNCounter[13]) fprintf(stderr, "Warning: No neutron could reach Component[13] PSDmon1Chopper\n");
    if (mcAbsorbProp[13]) fprintf(stderr, "Warning: %g events were removed in Component[13] PSDmon1Chopper=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[13]);
    if (!mcNCounter[14]) fprintf(stderr, "Warning: No neutron could reach Component[14] VacuumTube1entry\n");
    if (mcAbsorbProp[14]) fprintf(stderr, "Warning: %g events were removed in Component[14] VacuumTube1entry=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[14]);
    if (!mcNCounter[15]) fprintf(stderr, "Warning: No neutron could reach Component[15] VacuumTube1exit\n");
    if (mcAbsorbProp[15]) fprintf(stderr, "Warning: %g events were removed in Component[15] VacuumTube1exit=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[15]);
    if (!mcNCounter[16]) fprintf(stderr, "Warning: No neutron could reach Component[16] VacuumTube2entry\n");
    if (mcAbsorbProp[16]) fprintf(stderr, "Warning: %g events were removed in Component[16] VacuumTube2entry=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[16]);
    if (!mcNCounter[17]) fprintf(stderr, "Warning: No neutron could reach Component[17] VacuumTube2exit\n");
    if (mcAbsorbProp[17]) fprintf(stderr, "Warning: %g events were removed in Component[17] VacuumTube2exit=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[17]);
    if (!mcNCounter[18]) fprintf(stderr, "Warning: No neutron could reach Component[18] VacuumTube3entry\n");
    if (mcAbsorbProp[18]) fprintf(stderr, "Warning: %g events were removed in Component[18] VacuumTube3entry=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[18]);
    if (!mcNCounter[19]) fprintf(stderr, "Warning: No neutron could reach Component[19] VacuumTube3exit\n");
    if (mcAbsorbProp[19]) fprintf(stderr, "Warning: %g events were removed in Component[19] VacuumTube3exit=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[19]);
    if (!mcNCounter[20]) fprintf(stderr, "Warning: No neutron could reach Component[20] PSDmonMono\n");
    if (mcAbsorbProp[20]) fprintf(stderr, "Warning: %g events were removed in Component[20] PSDmonMono=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[20]);
  /* User FINALLY code for component 'MONOTOF'. */
  SIG_MESSAGE("MONOTOF (Finally)");
#define mccompcurname  MONOTOF
#define mccompcurtype  Monitor_nD
#define mccompcurindex 21
#define user1 mccMONOTOF_user1
#define user2 mccMONOTOF_user2
#define user3 mccMONOTOF_user3
#define DEFS mccMONOTOF_DEFS
#define Vars mccMONOTOF_Vars
#define detector mccMONOTOF_detector
#define offdata mccMONOTOF_offdata
{   /* Declarations of MONOTOF=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccMONOTOF_xwidth;
MCNUM yheight = mccMONOTOF_yheight;
MCNUM zdepth = mccMONOTOF_zdepth;
MCNUM xmin = mccMONOTOF_xmin;
MCNUM xmax = mccMONOTOF_xmax;
MCNUM ymin = mccMONOTOF_ymin;
MCNUM ymax = mccMONOTOF_ymax;
MCNUM zmin = mccMONOTOF_zmin;
MCNUM zmax = mccMONOTOF_zmax;
MCNUM bins = mccMONOTOF_bins;
MCNUM min = mccMONOTOF_min;
MCNUM max = mccMONOTOF_max;
MCNUM restore_neutron = mccMONOTOF_restore_neutron;
MCNUM radius = mccMONOTOF_radius;
char* options = mccMONOTOF_options;
char* filename = mccMONOTOF_filename;
char* geometry = mccMONOTOF_geometry;
char* username1 = mccMONOTOF_username1;
char* username2 = mccMONOTOF_username2;
char* username3 = mccMONOTOF_username3;
int nowritefile = mccMONOTOF_nowritefile;
#line 485 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
{
  /* free pointers */
  if (!nowritefile) {
    Monitor_nD_Finally(&DEFS, &Vars);
  }
}
#line 29566 "PSI_Focus.c"
}   /* End of MONOTOF=Monitor_nD() SETTING parameter declarations. */
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

    if (!mcNCounter[21]) fprintf(stderr, "Warning: No neutron could reach Component[21] MONOTOF\n");
    if (mcAbsorbProp[21]) fprintf(stderr, "Warning: %g events were removed in Component[21] MONOTOF=Monitor_nD()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[21]);
    if (!mcNCounter[22]) fprintf(stderr, "Warning: No neutron could reach Component[22] DivMonMono\n");
    if (mcAbsorbProp[22]) fprintf(stderr, "Warning: %g events were removed in Component[22] DivMonMono=Divergence_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[22]);
    if (!mcNCounter[23]) fprintf(stderr, "Warning: No neutron could reach Component[23] focus_mono\n");
    if (mcAbsorbProp[23]) fprintf(stderr, "Warning: %g events were removed in Component[23] focus_mono=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[23]);
    if (!mcNCounter[24]) fprintf(stderr, "Warning: No neutron could reach Component[24] mono\n");
    if (mcNCounter[24] < 1000*(10)) fprintf(stderr, 
"Warning: Number of events %g reaching SPLIT position Component[24] mono=Monochromator_2foc()\n"
"         is probably too low. Increase Ncount.\n", mcNCounter[24]);

    if (mcAbsorbProp[24]) fprintf(stderr, "Warning: %g events were removed in Component[24] mono=Monochromator_2foc()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[24]);
    if (!mcNCounter[25]) fprintf(stderr, "Warning: No neutron could reach Component[25] a2\n");
    if (mcAbsorbProp[25]) fprintf(stderr, "Warning: %g events were removed in Component[25] a2=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[25]);
  /* User FINALLY code for component 'FERMITOF_before'. */
  SIG_MESSAGE("FERMITOF_before (Finally)");
#define mccompcurname  FERMITOF_before
#define mccompcurtype  Monitor_nD
#define mccompcurindex 26
#define user1 mccFERMITOF_before_user1
#define user2 mccFERMITOF_before_user2
#define user3 mccFERMITOF_before_user3
#define DEFS mccFERMITOF_before_DEFS
#define Vars mccFERMITOF_before_Vars
#define detector mccFERMITOF_before_detector
#define offdata mccFERMITOF_before_offdata
{   /* Declarations of FERMITOF_before=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccFERMITOF_before_xwidth;
MCNUM yheight = mccFERMITOF_before_yheight;
MCNUM zdepth = mccFERMITOF_before_zdepth;
MCNUM xmin = mccFERMITOF_before_xmin;
MCNUM xmax = mccFERMITOF_before_xmax;
MCNUM ymin = mccFERMITOF_before_ymin;
MCNUM ymax = mccFERMITOF_before_ymax;
MCNUM zmin = mccFERMITOF_before_zmin;
MCNUM zmax = mccFERMITOF_before_zmax;
MCNUM bins = mccFERMITOF_before_bins;
MCNUM min = mccFERMITOF_before_min;
MCNUM max = mccFERMITOF_before_max;
MCNUM restore_neutron = mccFERMITOF_before_restore_neutron;
MCNUM radius = mccFERMITOF_before_radius;
char* options = mccFERMITOF_before_options;
char* filename = mccFERMITOF_before_filename;
char* geometry = mccFERMITOF_before_geometry;
char* username1 = mccFERMITOF_before_username1;
char* username2 = mccFERMITOF_before_username2;
char* username3 = mccFERMITOF_before_username3;
int nowritefile = mccFERMITOF_before_nowritefile;
#line 485 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
{
  /* free pointers */
  if (!nowritefile) {
    Monitor_nD_Finally(&DEFS, &Vars);
  }
}
#line 29631 "PSI_Focus.c"
}   /* End of FERMITOF_before=Monitor_nD() SETTING parameter declarations. */
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

    if (!mcNCounter[26]) fprintf(stderr, "Warning: No neutron could reach Component[26] FERMITOF_before\n");
    if (mcAbsorbProp[26]) fprintf(stderr, "Warning: %g events were removed in Component[26] FERMITOF_before=Monitor_nD()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[26]);
    if (!mcNCounter[27]) fprintf(stderr, "Warning: No neutron could reach Component[27] lambdaFermi\n");
    if (mcAbsorbProp[27]) fprintf(stderr, "Warning: %g events were removed in Component[27] lambdaFermi=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[27]);
    if (!mcNCounter[28]) fprintf(stderr, "Warning: No neutron could reach Component[28] EMON_Fermi\n");
    if (mcAbsorbProp[28]) fprintf(stderr, "Warning: %g events were removed in Component[28] EMON_Fermi=E_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[28]);
    if (!mcNCounter[29]) fprintf(stderr, "Warning: No neutron could reach Component[29] DivMonfermi1\n");
    if (mcAbsorbProp[29]) fprintf(stderr, "Warning: %g events were removed in Component[29] DivMonfermi1=Divergence_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[29]);
    if (!mcNCounter[30]) fprintf(stderr, "Warning: No neutron could reach Component[30] PSD_Fermi1\n");
    if (mcAbsorbProp[30]) fprintf(stderr, "Warning: %g events were removed in Component[30] PSD_Fermi1=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[30]);
    if (!mcNCounter[31]) fprintf(stderr, "Warning: No neutron could reach Component[31] FoChopper\n");
    if (mcAbsorbProp[31]) fprintf(stderr, "Warning: %g events were removed in Component[31] FoChopper=FermiChopper()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[31]);
    if (!mcNCounter[32]) fprintf(stderr, "Warning: No neutron could reach Component[32] PSD_Fermi2\n");
    if (mcAbsorbProp[32]) fprintf(stderr, "Warning: %g events were removed in Component[32] PSD_Fermi2=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[32]);
    if (!mcNCounter[33]) fprintf(stderr, "Warning: No neutron could reach Component[33] DivMonfermi2\n");
    if (mcAbsorbProp[33]) fprintf(stderr, "Warning: %g events were removed in Component[33] DivMonfermi2=Divergence_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[33]);
  /* User FINALLY code for component 'FERMITOF1'. */
  SIG_MESSAGE("FERMITOF1 (Finally)");
#define mccompcurname  FERMITOF1
#define mccompcurtype  Monitor_nD
#define mccompcurindex 34
#define user1 mccFERMITOF1_user1
#define user2 mccFERMITOF1_user2
#define user3 mccFERMITOF1_user3
#define DEFS mccFERMITOF1_DEFS
#define Vars mccFERMITOF1_Vars
#define detector mccFERMITOF1_detector
#define offdata mccFERMITOF1_offdata
{   /* Declarations of FERMITOF1=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccFERMITOF1_xwidth;
MCNUM yheight = mccFERMITOF1_yheight;
MCNUM zdepth = mccFERMITOF1_zdepth;
MCNUM xmin = mccFERMITOF1_xmin;
MCNUM xmax = mccFERMITOF1_xmax;
MCNUM ymin = mccFERMITOF1_ymin;
MCNUM ymax = mccFERMITOF1_ymax;
MCNUM zmin = mccFERMITOF1_zmin;
MCNUM zmax = mccFERMITOF1_zmax;
MCNUM bins = mccFERMITOF1_bins;
MCNUM min = mccFERMITOF1_min;
MCNUM max = mccFERMITOF1_max;
MCNUM restore_neutron = mccFERMITOF1_restore_neutron;
MCNUM radius = mccFERMITOF1_radius;
char* options = mccFERMITOF1_options;
char* filename = mccFERMITOF1_filename;
char* geometry = mccFERMITOF1_geometry;
char* username1 = mccFERMITOF1_username1;
char* username2 = mccFERMITOF1_username2;
char* username3 = mccFERMITOF1_username3;
int nowritefile = mccFERMITOF1_nowritefile;
#line 485 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
{
  /* free pointers */
  if (!nowritefile) {
    Monitor_nD_Finally(&DEFS, &Vars);
  }
}
#line 29701 "PSI_Focus.c"
}   /* End of FERMITOF1=Monitor_nD() SETTING parameter declarations. */
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

    if (!mcNCounter[34]) fprintf(stderr, "Warning: No neutron could reach Component[34] FERMITOF1\n");
    if (mcAbsorbProp[34]) fprintf(stderr, "Warning: %g events were removed in Component[34] FERMITOF1=Monitor_nD()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[34]);
    if (!mcNCounter[35]) fprintf(stderr, "Warning: No neutron could reach Component[35] SAMPLE_SLIT\n");
    if (mcAbsorbProp[35]) fprintf(stderr, "Warning: %g events were removed in Component[35] SAMPLE_SLIT=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[35]);
  /* User FINALLY code for component 'FERMITOF2'. */
  SIG_MESSAGE("FERMITOF2 (Finally)");
#define mccompcurname  FERMITOF2
#define mccompcurtype  Monitor_nD
#define mccompcurindex 36
#define user1 mccFERMITOF2_user1
#define user2 mccFERMITOF2_user2
#define user3 mccFERMITOF2_user3
#define DEFS mccFERMITOF2_DEFS
#define Vars mccFERMITOF2_Vars
#define detector mccFERMITOF2_detector
#define offdata mccFERMITOF2_offdata
{   /* Declarations of FERMITOF2=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccFERMITOF2_xwidth;
MCNUM yheight = mccFERMITOF2_yheight;
MCNUM zdepth = mccFERMITOF2_zdepth;
MCNUM xmin = mccFERMITOF2_xmin;
MCNUM xmax = mccFERMITOF2_xmax;
MCNUM ymin = mccFERMITOF2_ymin;
MCNUM ymax = mccFERMITOF2_ymax;
MCNUM zmin = mccFERMITOF2_zmin;
MCNUM zmax = mccFERMITOF2_zmax;
MCNUM bins = mccFERMITOF2_bins;
MCNUM min = mccFERMITOF2_min;
MCNUM max = mccFERMITOF2_max;
MCNUM restore_neutron = mccFERMITOF2_restore_neutron;
MCNUM radius = mccFERMITOF2_radius;
char* options = mccFERMITOF2_options;
char* filename = mccFERMITOF2_filename;
char* geometry = mccFERMITOF2_geometry;
char* username1 = mccFERMITOF2_username1;
char* username2 = mccFERMITOF2_username2;
char* username3 = mccFERMITOF2_username3;
int nowritefile = mccFERMITOF2_nowritefile;
#line 485 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
{
  /* free pointers */
  if (!nowritefile) {
    Monitor_nD_Finally(&DEFS, &Vars);
  }
}
#line 29759 "PSI_Focus.c"
}   /* End of FERMITOF2=Monitor_nD() SETTING parameter declarations. */
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

    if (!mcNCounter[36]) fprintf(stderr, "Warning: No neutron could reach Component[36] FERMITOF2\n");
    if (mcAbsorbProp[36]) fprintf(stderr, "Warning: %g events were removed in Component[36] FERMITOF2=Monitor_nD()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[36]);
    if (!mcNCounter[37]) fprintf(stderr, "Warning: No neutron could reach Component[37] PSD_SAMPLE\n");
    if (mcAbsorbProp[37]) fprintf(stderr, "Warning: %g events were removed in Component[37] PSD_SAMPLE=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[37]);
    if (!mcNCounter[38]) fprintf(stderr, "Warning: No neutron could reach Component[38] DivMon_Sample\n");
    if (mcAbsorbProp[38]) fprintf(stderr, "Warning: %g events were removed in Component[38] DivMon_Sample=Divergence_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[38]);
    if (!mcNCounter[39]) fprintf(stderr, "Warning: No neutron could reach Component[39] EMON_SAMPLE\n");
    if (mcAbsorbProp[39]) fprintf(stderr, "Warning: %g events were removed in Component[39] EMON_SAMPLE=E_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[39]);
    if (!mcNCounter[40]) fprintf(stderr, "Warning: No neutron could reach Component[40] a3\n");
    if (mcAbsorbProp[40]) fprintf(stderr, "Warning: %g events were removed in Component[40] a3=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[40]);
    if (!mcNCounter[41]) fprintf(stderr, "Warning: No neutron could reach Component[41] Sample\n");
    if (mcNCounter[41] < 1000*(10)) fprintf(stderr, 
"Warning: Number of events %g reaching SPLIT position Component[41] Sample=V_sample()\n"
"         is probably too low. Increase Ncount.\n", mcNCounter[41]);

    if (mcAbsorbProp[41]) fprintf(stderr, "Warning: %g events were removed in Component[41] Sample=V_sample()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[41]);
  /* User FINALLY code for component 'TOF_Det'. */
  SIG_MESSAGE("TOF_Det (Finally)");
#define mccompcurname  TOF_Det
#define mccompcurtype  Monitor_nD
#define mccompcurindex 42
#define user1 mccTOF_Det_user1
#define user2 mccTOF_Det_user2
#define user3 mccTOF_Det_user3
#define DEFS mccTOF_Det_DEFS
#define Vars mccTOF_Det_Vars
#define detector mccTOF_Det_detector
#define offdata mccTOF_Det_offdata
{   /* Declarations of TOF_Det=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccTOF_Det_xwidth;
MCNUM yheight = mccTOF_Det_yheight;
MCNUM zdepth = mccTOF_Det_zdepth;
MCNUM xmin = mccTOF_Det_xmin;
MCNUM xmax = mccTOF_Det_xmax;
MCNUM ymin = mccTOF_Det_ymin;
MCNUM ymax = mccTOF_Det_ymax;
MCNUM zmin = mccTOF_Det_zmin;
MCNUM zmax = mccTOF_Det_zmax;
MCNUM bins = mccTOF_Det_bins;
MCNUM min = mccTOF_Det_min;
MCNUM max = mccTOF_Det_max;
MCNUM restore_neutron = mccTOF_Det_restore_neutron;
MCNUM radius = mccTOF_Det_radius;
char* options = mccTOF_Det_options;
char* filename = mccTOF_Det_filename;
char* geometry = mccTOF_Det_geometry;
char* username1 = mccTOF_Det_username1;
char* username2 = mccTOF_Det_username2;
char* username3 = mccTOF_Det_username3;
int nowritefile = mccTOF_Det_nowritefile;
#line 485 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
{
  /* free pointers */
  if (!nowritefile) {
    Monitor_nD_Finally(&DEFS, &Vars);
  }
}
#line 29826 "PSI_Focus.c"
}   /* End of TOF_Det=Monitor_nD() SETTING parameter declarations. */
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

    if (!mcNCounter[42]) fprintf(stderr, "Warning: No neutron could reach Component[42] TOF_Det\n");
    if (mcAbsorbProp[42]) fprintf(stderr, "Warning: %g events were removed in Component[42] TOF_Det=Monitor_nD()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[42]);
  /* User FINALLY code for component 'FoDet'. */
  SIG_MESSAGE("FoDet (Finally)");
#define mccompcurname  FoDet
#define mccompcurtype  Monitor_nD
#define mccompcurindex 43
#define user1 mccFoDet_user1
#define user2 mccFoDet_user2
#define user3 mccFoDet_user3
#define DEFS mccFoDet_DEFS
#define Vars mccFoDet_Vars
#define detector mccFoDet_detector
#define offdata mccFoDet_offdata
{   /* Declarations of FoDet=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccFoDet_xwidth;
MCNUM yheight = mccFoDet_yheight;
MCNUM zdepth = mccFoDet_zdepth;
MCNUM xmin = mccFoDet_xmin;
MCNUM xmax = mccFoDet_xmax;
MCNUM ymin = mccFoDet_ymin;
MCNUM ymax = mccFoDet_ymax;
MCNUM zmin = mccFoDet_zmin;
MCNUM zmax = mccFoDet_zmax;
MCNUM bins = mccFoDet_bins;
MCNUM min = mccFoDet_min;
MCNUM max = mccFoDet_max;
MCNUM restore_neutron = mccFoDet_restore_neutron;
MCNUM radius = mccFoDet_radius;
char* options = mccFoDet_options;
char* filename = mccFoDet_filename;
char* geometry = mccFoDet_geometry;
char* username1 = mccFoDet_username1;
char* username2 = mccFoDet_username2;
char* username3 = mccFoDet_username3;
int nowritefile = mccFoDet_nowritefile;
#line 485 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
{
  /* free pointers */
  if (!nowritefile) {
    Monitor_nD_Finally(&DEFS, &Vars);
  }
}
#line 29882 "PSI_Focus.c"
}   /* End of FoDet=Monitor_nD() SETTING parameter declarations. */
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

    if (!mcNCounter[43]) fprintf(stderr, "Warning: No neutron could reach Component[43] FoDet\n");
    if (mcAbsorbProp[43]) fprintf(stderr, "Warning: %g events were removed in Component[43] FoDet=Monitor_nD()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[43]);
    if (!mcNCounter[44]) fprintf(stderr, "Warning: No neutron could reach Component[44] EMON_DET\n");
    if (mcAbsorbProp[44]) fprintf(stderr, "Warning: %g events were removed in Component[44] EMON_DET=E_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[44]);
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
#line 147 "/usr/share/mcstas/2.5/misc/Progress_bar.comp"
{
  
}
#line 29933 "PSI_Focus.c"
}   /* End of a1=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'csource'. */
  SIG_MESSAGE("csource (McDisplay)");
  printf("MCDISPLAY: component %s\n", "csource");
#define mccompcurname  csource
#define mccompcurtype  Source_gen
#define mccompcurindex 2
#define p_in mcccsource_p_in
#define lambda1 mcccsource_lambda1
#define lambda2 mcccsource_lambda2
#define lambda3 mcccsource_lambda3
#define pTable mcccsource_pTable
#define pTable_x mcccsource_pTable_x
#define pTable_y mcccsource_pTable_y
#define pTable_xmin mcccsource_pTable_xmin
#define pTable_xmax mcccsource_pTable_xmax
#define pTable_xsum mcccsource_pTable_xsum
#define pTable_ymin mcccsource_pTable_ymin
#define pTable_ymax mcccsource_pTable_ymax
#define pTable_ysum mcccsource_pTable_ysum
#define pTable_dxmin mcccsource_pTable_dxmin
#define pTable_dxmax mcccsource_pTable_dxmax
#define pTable_dymin mcccsource_pTable_dymin
#define pTable_dymax mcccsource_pTable_dymax
{   /* Declarations of csource=Source_gen() SETTING parameters. */
char* flux_file = mcccsource_flux_file;
char* xdiv_file = mcccsource_xdiv_file;
char* ydiv_file = mcccsource_ydiv_file;
MCNUM radius = mcccsource_radius;
MCNUM dist = mcccsource_dist;
MCNUM focus_xw = mcccsource_focus_xw;
MCNUM focus_yh = mcccsource_focus_yh;
MCNUM focus_aw = mcccsource_focus_aw;
MCNUM focus_ah = mcccsource_focus_ah;
MCNUM E0 = mcccsource_E0;
MCNUM dE = mcccsource_dE;
MCNUM lambda0 = mcccsource_lambda0;
MCNUM dlambda = mcccsource_dlambda;
MCNUM I1 = mcccsource_I1;
MCNUM yheight = mcccsource_yheight;
MCNUM xwidth = mcccsource_xwidth;
MCNUM verbose = mcccsource_verbose;
MCNUM T1 = mcccsource_T1;
MCNUM flux_file_perAA = mcccsource_flux_file_perAA;
MCNUM flux_file_log = mcccsource_flux_file_log;
MCNUM Lmin = mcccsource_Lmin;
MCNUM Lmax = mcccsource_Lmax;
MCNUM Emin = mcccsource_Emin;
MCNUM Emax = mcccsource_Emax;
MCNUM T2 = mcccsource_T2;
MCNUM I2 = mcccsource_I2;
MCNUM T3 = mcccsource_T3;
MCNUM I3 = mcccsource_I3;
MCNUM zdepth = mcccsource_zdepth;
int target_index = mcccsource_target_index;
#line 578 "/usr/share/mcstas/2.5/sources/Source_gen.comp"
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
#line 30046 "PSI_Focus.c"
}   /* End of csource=Source_gen() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'guide1'. */
  SIG_MESSAGE("guide1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "guide1");
#define mccompcurname  guide1
#define mccompcurtype  Guide
#define mccompcurindex 3
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
#line 202 "/usr/share/mcstas/2.5/optics/Guide.comp"
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
#line 30108 "PSI_Focus.c"
}   /* End of guide1=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'guide2'. */
  SIG_MESSAGE("guide2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "guide2");
#define mccompcurname  guide2
#define mccompcurtype  Bender
#define mccompcurindex 4
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
#line 269 "/usr/share/mcstas/2.5/optics/Bender.comp"
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
#line 30174 "PSI_Focus.c"
}   /* End of guide2=Bender() SETTING parameter declarations. */
#undef mWin
#undef bk
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'bunker'. */
  SIG_MESSAGE("bunker (McDisplay)");
  printf("MCDISPLAY: component %s\n", "bunker");
#define mccompcurname  bunker
#define mccompcurtype  Guide
#define mccompcurindex 5
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
#line 202 "/usr/share/mcstas/2.5/optics/Guide.comp"
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
#line 30221 "PSI_Focus.c"
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
#define mccompcurindex 6
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
#line 202 "/usr/share/mcstas/2.5/optics/Guide.comp"
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
#line 30267 "PSI_Focus.c"
}   /* End of guide3=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'lambdaGuideExit'. */
  SIG_MESSAGE("lambdaGuideExit (McDisplay)");
  printf("MCDISPLAY: component %s\n", "lambdaGuideExit");
#define mccompcurname  lambdaGuideExit
#define mccompcurtype  L_monitor
#define mccompcurindex 7
#define nL mcclambdaGuideExit_nL
#define L_N mcclambdaGuideExit_L_N
#define L_p mcclambdaGuideExit_L_p
#define L_p2 mcclambdaGuideExit_L_p2
{   /* Declarations of lambdaGuideExit=L_monitor() SETTING parameters. */
char* filename = mcclambdaGuideExit_filename;
MCNUM xmin = mcclambdaGuideExit_xmin;
MCNUM xmax = mcclambdaGuideExit_xmax;
MCNUM ymin = mcclambdaGuideExit_ymin;
MCNUM ymax = mcclambdaGuideExit_ymax;
MCNUM xwidth = mcclambdaGuideExit_xwidth;
MCNUM yheight = mcclambdaGuideExit_yheight;
MCNUM Lmin = mcclambdaGuideExit_Lmin;
MCNUM Lmax = mcclambdaGuideExit_Lmax;
MCNUM restore_neutron = mcclambdaGuideExit_restore_neutron;
int nowritefile = mcclambdaGuideExit_nowritefile;
#line 120 "/usr/share/mcstas/2.5/monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 30305 "PSI_Focus.c"
}   /* End of lambdaGuideExit=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'DivMonGuideExit'. */
  SIG_MESSAGE("DivMonGuideExit (McDisplay)");
  printf("MCDISPLAY: component %s\n", "DivMonGuideExit");
#define mccompcurname  DivMonGuideExit
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 8
#define nh mccDivMonGuideExit_nh
#define nv mccDivMonGuideExit_nv
#define Div_N mccDivMonGuideExit_Div_N
#define Div_p mccDivMonGuideExit_Div_p
#define Div_p2 mccDivMonGuideExit_Div_p2
{   /* Declarations of DivMonGuideExit=Divergence_monitor() SETTING parameters. */
char* filename = mccDivMonGuideExit_filename;
MCNUM xmin = mccDivMonGuideExit_xmin;
MCNUM xmax = mccDivMonGuideExit_xmax;
MCNUM ymin = mccDivMonGuideExit_ymin;
MCNUM ymax = mccDivMonGuideExit_ymax;
MCNUM xwidth = mccDivMonGuideExit_xwidth;
MCNUM yheight = mccDivMonGuideExit_yheight;
MCNUM maxdiv_h = mccDivMonGuideExit_maxdiv_h;
MCNUM maxdiv_v = mccDivMonGuideExit_maxdiv_v;
MCNUM restore_neutron = mccDivMonGuideExit_restore_neutron;
MCNUM nx = mccDivMonGuideExit_nx;
MCNUM ny = mccDivMonGuideExit_ny;
MCNUM nz = mccDivMonGuideExit_nz;
int nowritefile = mccDivMonGuideExit_nowritefile;
#line 131 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
{
    
    multiline(5, (double)xmin, (double)ymin, 0.0,
                 (double)xmax, (double)ymin, 0.0,
                 (double)xmax, (double)ymax, 0.0,
                 (double)xmin, (double)ymax, 0.0,
                 (double)xmin, (double)ymin, 0.0);
}
#line 30350 "PSI_Focus.c"
}   /* End of DivMonGuideExit=Divergence_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PSDGuideExit'. */
  SIG_MESSAGE("PSDGuideExit (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PSDGuideExit");
#define mccompcurname  PSDGuideExit
#define mccompcurtype  PSD_monitor
#define mccompcurindex 9
#define nx mccPSDGuideExit_nx
#define ny mccPSDGuideExit_ny
#define PSD_N mccPSDGuideExit_PSD_N
#define PSD_p mccPSDGuideExit_PSD_p
#define PSD_p2 mccPSDGuideExit_PSD_p2
{   /* Declarations of PSDGuideExit=PSD_monitor() SETTING parameters. */
char* filename = mccPSDGuideExit_filename;
MCNUM xmin = mccPSDGuideExit_xmin;
MCNUM xmax = mccPSDGuideExit_xmax;
MCNUM ymin = mccPSDGuideExit_ymin;
MCNUM ymax = mccPSDGuideExit_ymax;
MCNUM xwidth = mccPSDGuideExit_xwidth;
MCNUM yheight = mccPSDGuideExit_yheight;
MCNUM restore_neutron = mccPSDGuideExit_restore_neutron;
int nowritefile = mccPSDGuideExit_nowritefile;
#line 115 "/usr/share/mcstas/2.5/monitors/PSD_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 30391 "PSI_Focus.c"
}   /* End of PSDGuideExit=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'FOCUSguide'. */
  SIG_MESSAGE("FOCUSguide (McDisplay)");
  printf("MCDISPLAY: component %s\n", "FOCUSguide");
#define mccompcurname  FOCUSguide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 10
#define w1c mccFOCUSguide_w1c
#define w2c mccFOCUSguide_w2c
#define ww mccFOCUSguide_ww
#define hh mccFOCUSguide_hh
#define whalf mccFOCUSguide_whalf
#define hhalf mccFOCUSguide_hhalf
#define lwhalf mccFOCUSguide_lwhalf
#define lhhalf mccFOCUSguide_lhhalf
{   /* Declarations of FOCUSguide=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccFOCUSguide_w1;
MCNUM h1 = mccFOCUSguide_h1;
MCNUM w2 = mccFOCUSguide_w2;
MCNUM h2 = mccFOCUSguide_h2;
MCNUM l = mccFOCUSguide_l;
MCNUM R0 = mccFOCUSguide_R0;
MCNUM Qc = mccFOCUSguide_Qc;
MCNUM alpha = mccFOCUSguide_alpha;
MCNUM m = mccFOCUSguide_m;
MCNUM nslit = mccFOCUSguide_nslit;
MCNUM d = mccFOCUSguide_d;
MCNUM Qcx = mccFOCUSguide_Qcx;
MCNUM Qcy = mccFOCUSguide_Qcy;
MCNUM alphax = mccFOCUSguide_alphax;
MCNUM alphay = mccFOCUSguide_alphay;
MCNUM W = mccFOCUSguide_W;
MCNUM mx = mccFOCUSguide_mx;
MCNUM my = mccFOCUSguide_my;
MCNUM nu = mccFOCUSguide_nu;
MCNUM phase = mccFOCUSguide_phase;
#line 297 "/usr/share/mcstas/2.5/optics/Guide_channeled.comp"
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
#line 30468 "PSI_Focus.c"
}   /* End of FOCUSguide=Guide_channeled() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'FirstChopper'. */
  SIG_MESSAGE("FirstChopper (McDisplay)");
  printf("MCDISPLAY: component %s\n", "FirstChopper");
#define mccompcurname  FirstChopper
#define mccompcurtype  DiskChopper
#define mccompcurindex 11
#define Tg mccFirstChopper_Tg
#define To mccFirstChopper_To
#define delta_y mccFirstChopper_delta_y
#define height mccFirstChopper_height
#define omega mccFirstChopper_omega
{   /* Declarations of FirstChopper=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccFirstChopper_theta_0;
MCNUM radius = mccFirstChopper_radius;
MCNUM yheight = mccFirstChopper_yheight;
MCNUM nu = mccFirstChopper_nu;
MCNUM nslit = mccFirstChopper_nslit;
MCNUM jitter = mccFirstChopper_jitter;
MCNUM delay = mccFirstChopper_delay;
MCNUM isfirst = mccFirstChopper_isfirst;
MCNUM n_pulse = mccFirstChopper_n_pulse;
MCNUM abs_out = mccFirstChopper_abs_out;
MCNUM phase = mccFirstChopper_phase;
MCNUM xwidth = mccFirstChopper_xwidth;
MCNUM verbose = mccFirstChopper_verbose;
#line 168 "/usr/share/mcstas/2.5/optics/DiskChopper.comp"
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
#line 30534 "PSI_Focus.c"
}   /* End of FirstChopper=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'DISCTOF'. */
  SIG_MESSAGE("DISCTOF (McDisplay)");
  printf("MCDISPLAY: component %s\n", "DISCTOF");
#define mccompcurname  DISCTOF
#define mccompcurtype  Monitor_nD
#define mccompcurindex 12
#define user1 mccDISCTOF_user1
#define user2 mccDISCTOF_user2
#define user3 mccDISCTOF_user3
#define DEFS mccDISCTOF_DEFS
#define Vars mccDISCTOF_Vars
#define detector mccDISCTOF_detector
#define offdata mccDISCTOF_offdata
{   /* Declarations of DISCTOF=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccDISCTOF_xwidth;
MCNUM yheight = mccDISCTOF_yheight;
MCNUM zdepth = mccDISCTOF_zdepth;
MCNUM xmin = mccDISCTOF_xmin;
MCNUM xmax = mccDISCTOF_xmax;
MCNUM ymin = mccDISCTOF_ymin;
MCNUM ymax = mccDISCTOF_ymax;
MCNUM zmin = mccDISCTOF_zmin;
MCNUM zmax = mccDISCTOF_zmax;
MCNUM bins = mccDISCTOF_bins;
MCNUM min = mccDISCTOF_min;
MCNUM max = mccDISCTOF_max;
MCNUM restore_neutron = mccDISCTOF_restore_neutron;
MCNUM radius = mccDISCTOF_radius;
char* options = mccDISCTOF_options;
char* filename = mccDISCTOF_filename;
char* geometry = mccDISCTOF_geometry;
char* username1 = mccDISCTOF_username1;
char* username2 = mccDISCTOF_username2;
char* username3 = mccDISCTOF_username3;
int nowritefile = mccDISCTOF_nowritefile;
#line 493 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
{
  if (geometry && strlen(geometry) && strcmp(geometry,"0") && strcmp(geometry, "NULL"))
  {
    off_display(offdata);
  } else {
    Monitor_nD_McDisplay(&DEFS, &Vars);
  }
}
#line 30589 "PSI_Focus.c"
}   /* End of DISCTOF=Monitor_nD() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'PSDmon1Chopper'. */
  SIG_MESSAGE("PSDmon1Chopper (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PSDmon1Chopper");
#define mccompcurname  PSDmon1Chopper
#define mccompcurtype  PSD_monitor
#define mccompcurindex 13
#define nx mccPSDmon1Chopper_nx
#define ny mccPSDmon1Chopper_ny
#define PSD_N mccPSDmon1Chopper_PSD_N
#define PSD_p mccPSDmon1Chopper_PSD_p
#define PSD_p2 mccPSDmon1Chopper_PSD_p2
{   /* Declarations of PSDmon1Chopper=PSD_monitor() SETTING parameters. */
char* filename = mccPSDmon1Chopper_filename;
MCNUM xmin = mccPSDmon1Chopper_xmin;
MCNUM xmax = mccPSDmon1Chopper_xmax;
MCNUM ymin = mccPSDmon1Chopper_ymin;
MCNUM ymax = mccPSDmon1Chopper_ymax;
MCNUM xwidth = mccPSDmon1Chopper_xwidth;
MCNUM yheight = mccPSDmon1Chopper_yheight;
MCNUM restore_neutron = mccPSDmon1Chopper_restore_neutron;
int nowritefile = mccPSDmon1Chopper_nowritefile;
#line 115 "/usr/share/mcstas/2.5/monitors/PSD_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 30632 "PSI_Focus.c"
}   /* End of PSDmon1Chopper=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'VacuumTube1entry'. */
  SIG_MESSAGE("VacuumTube1entry (McDisplay)");
  printf("MCDISPLAY: component %s\n", "VacuumTube1entry");
#define mccompcurname  VacuumTube1entry
#define mccompcurtype  Slit
#define mccompcurindex 14
{   /* Declarations of VacuumTube1entry=Slit() SETTING parameters. */
MCNUM xmin = mccVacuumTube1entry_xmin;
MCNUM xmax = mccVacuumTube1entry_xmax;
MCNUM ymin = mccVacuumTube1entry_ymin;
MCNUM ymax = mccVacuumTube1entry_ymax;
MCNUM radius = mccVacuumTube1entry_radius;
MCNUM xwidth = mccVacuumTube1entry_xwidth;
MCNUM yheight = mccVacuumTube1entry_yheight;
#line 81 "/usr/share/mcstas/2.5/optics/Slit.comp"
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
#line 30680 "PSI_Focus.c"
}   /* End of VacuumTube1entry=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'VacuumTube1exit'. */
  SIG_MESSAGE("VacuumTube1exit (McDisplay)");
  printf("MCDISPLAY: component %s\n", "VacuumTube1exit");
#define mccompcurname  VacuumTube1exit
#define mccompcurtype  Slit
#define mccompcurindex 15
{   /* Declarations of VacuumTube1exit=Slit() SETTING parameters. */
MCNUM xmin = mccVacuumTube1exit_xmin;
MCNUM xmax = mccVacuumTube1exit_xmax;
MCNUM ymin = mccVacuumTube1exit_ymin;
MCNUM ymax = mccVacuumTube1exit_ymax;
MCNUM radius = mccVacuumTube1exit_radius;
MCNUM xwidth = mccVacuumTube1exit_xwidth;
MCNUM yheight = mccVacuumTube1exit_yheight;
#line 81 "/usr/share/mcstas/2.5/optics/Slit.comp"
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
#line 30723 "PSI_Focus.c"
}   /* End of VacuumTube1exit=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'VacuumTube2entry'. */
  SIG_MESSAGE("VacuumTube2entry (McDisplay)");
  printf("MCDISPLAY: component %s\n", "VacuumTube2entry");
#define mccompcurname  VacuumTube2entry
#define mccompcurtype  Slit
#define mccompcurindex 16
{   /* Declarations of VacuumTube2entry=Slit() SETTING parameters. */
MCNUM xmin = mccVacuumTube2entry_xmin;
MCNUM xmax = mccVacuumTube2entry_xmax;
MCNUM ymin = mccVacuumTube2entry_ymin;
MCNUM ymax = mccVacuumTube2entry_ymax;
MCNUM radius = mccVacuumTube2entry_radius;
MCNUM xwidth = mccVacuumTube2entry_xwidth;
MCNUM yheight = mccVacuumTube2entry_yheight;
#line 81 "/usr/share/mcstas/2.5/optics/Slit.comp"
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
#line 30766 "PSI_Focus.c"
}   /* End of VacuumTube2entry=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'VacuumTube2exit'. */
  SIG_MESSAGE("VacuumTube2exit (McDisplay)");
  printf("MCDISPLAY: component %s\n", "VacuumTube2exit");
#define mccompcurname  VacuumTube2exit
#define mccompcurtype  Slit
#define mccompcurindex 17
{   /* Declarations of VacuumTube2exit=Slit() SETTING parameters. */
MCNUM xmin = mccVacuumTube2exit_xmin;
MCNUM xmax = mccVacuumTube2exit_xmax;
MCNUM ymin = mccVacuumTube2exit_ymin;
MCNUM ymax = mccVacuumTube2exit_ymax;
MCNUM radius = mccVacuumTube2exit_radius;
MCNUM xwidth = mccVacuumTube2exit_xwidth;
MCNUM yheight = mccVacuumTube2exit_yheight;
#line 81 "/usr/share/mcstas/2.5/optics/Slit.comp"
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
#line 30809 "PSI_Focus.c"
}   /* End of VacuumTube2exit=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'VacuumTube3entry'. */
  SIG_MESSAGE("VacuumTube3entry (McDisplay)");
  printf("MCDISPLAY: component %s\n", "VacuumTube3entry");
#define mccompcurname  VacuumTube3entry
#define mccompcurtype  Slit
#define mccompcurindex 18
{   /* Declarations of VacuumTube3entry=Slit() SETTING parameters. */
MCNUM xmin = mccVacuumTube3entry_xmin;
MCNUM xmax = mccVacuumTube3entry_xmax;
MCNUM ymin = mccVacuumTube3entry_ymin;
MCNUM ymax = mccVacuumTube3entry_ymax;
MCNUM radius = mccVacuumTube3entry_radius;
MCNUM xwidth = mccVacuumTube3entry_xwidth;
MCNUM yheight = mccVacuumTube3entry_yheight;
#line 81 "/usr/share/mcstas/2.5/optics/Slit.comp"
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
#line 30852 "PSI_Focus.c"
}   /* End of VacuumTube3entry=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'VacuumTube3exit'. */
  SIG_MESSAGE("VacuumTube3exit (McDisplay)");
  printf("MCDISPLAY: component %s\n", "VacuumTube3exit");
#define mccompcurname  VacuumTube3exit
#define mccompcurtype  Slit
#define mccompcurindex 19
{   /* Declarations of VacuumTube3exit=Slit() SETTING parameters. */
MCNUM xmin = mccVacuumTube3exit_xmin;
MCNUM xmax = mccVacuumTube3exit_xmax;
MCNUM ymin = mccVacuumTube3exit_ymin;
MCNUM ymax = mccVacuumTube3exit_ymax;
MCNUM radius = mccVacuumTube3exit_radius;
MCNUM xwidth = mccVacuumTube3exit_xwidth;
MCNUM yheight = mccVacuumTube3exit_yheight;
#line 81 "/usr/share/mcstas/2.5/optics/Slit.comp"
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
#line 30895 "PSI_Focus.c"
}   /* End of VacuumTube3exit=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PSDmonMono'. */
  SIG_MESSAGE("PSDmonMono (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PSDmonMono");
#define mccompcurname  PSDmonMono
#define mccompcurtype  PSD_monitor
#define mccompcurindex 20
#define nx mccPSDmonMono_nx
#define ny mccPSDmonMono_ny
#define PSD_N mccPSDmonMono_PSD_N
#define PSD_p mccPSDmonMono_PSD_p
#define PSD_p2 mccPSDmonMono_PSD_p2
{   /* Declarations of PSDmonMono=PSD_monitor() SETTING parameters. */
char* filename = mccPSDmonMono_filename;
MCNUM xmin = mccPSDmonMono_xmin;
MCNUM xmax = mccPSDmonMono_xmax;
MCNUM ymin = mccPSDmonMono_ymin;
MCNUM ymax = mccPSDmonMono_ymax;
MCNUM xwidth = mccPSDmonMono_xwidth;
MCNUM yheight = mccPSDmonMono_yheight;
MCNUM restore_neutron = mccPSDmonMono_restore_neutron;
int nowritefile = mccPSDmonMono_nowritefile;
#line 115 "/usr/share/mcstas/2.5/monitors/PSD_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 30931 "PSI_Focus.c"
}   /* End of PSDmonMono=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'MONOTOF'. */
  SIG_MESSAGE("MONOTOF (McDisplay)");
  printf("MCDISPLAY: component %s\n", "MONOTOF");
#define mccompcurname  MONOTOF
#define mccompcurtype  Monitor_nD
#define mccompcurindex 21
#define user1 mccMONOTOF_user1
#define user2 mccMONOTOF_user2
#define user3 mccMONOTOF_user3
#define DEFS mccMONOTOF_DEFS
#define Vars mccMONOTOF_Vars
#define detector mccMONOTOF_detector
#define offdata mccMONOTOF_offdata
{   /* Declarations of MONOTOF=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccMONOTOF_xwidth;
MCNUM yheight = mccMONOTOF_yheight;
MCNUM zdepth = mccMONOTOF_zdepth;
MCNUM xmin = mccMONOTOF_xmin;
MCNUM xmax = mccMONOTOF_xmax;
MCNUM ymin = mccMONOTOF_ymin;
MCNUM ymax = mccMONOTOF_ymax;
MCNUM zmin = mccMONOTOF_zmin;
MCNUM zmax = mccMONOTOF_zmax;
MCNUM bins = mccMONOTOF_bins;
MCNUM min = mccMONOTOF_min;
MCNUM max = mccMONOTOF_max;
MCNUM restore_neutron = mccMONOTOF_restore_neutron;
MCNUM radius = mccMONOTOF_radius;
char* options = mccMONOTOF_options;
char* filename = mccMONOTOF_filename;
char* geometry = mccMONOTOF_geometry;
char* username1 = mccMONOTOF_username1;
char* username2 = mccMONOTOF_username2;
char* username3 = mccMONOTOF_username3;
int nowritefile = mccMONOTOF_nowritefile;
#line 493 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
{
  if (geometry && strlen(geometry) && strcmp(geometry,"0") && strcmp(geometry, "NULL"))
  {
    off_display(offdata);
  } else {
    Monitor_nD_McDisplay(&DEFS, &Vars);
  }
}
#line 30986 "PSI_Focus.c"
}   /* End of MONOTOF=Monitor_nD() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'DivMonMono'. */
  SIG_MESSAGE("DivMonMono (McDisplay)");
  printf("MCDISPLAY: component %s\n", "DivMonMono");
#define mccompcurname  DivMonMono
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 22
#define nh mccDivMonMono_nh
#define nv mccDivMonMono_nv
#define Div_N mccDivMonMono_Div_N
#define Div_p mccDivMonMono_Div_p
#define Div_p2 mccDivMonMono_Div_p2
{   /* Declarations of DivMonMono=Divergence_monitor() SETTING parameters. */
char* filename = mccDivMonMono_filename;
MCNUM xmin = mccDivMonMono_xmin;
MCNUM xmax = mccDivMonMono_xmax;
MCNUM ymin = mccDivMonMono_ymin;
MCNUM ymax = mccDivMonMono_ymax;
MCNUM xwidth = mccDivMonMono_xwidth;
MCNUM yheight = mccDivMonMono_yheight;
MCNUM maxdiv_h = mccDivMonMono_maxdiv_h;
MCNUM maxdiv_v = mccDivMonMono_maxdiv_v;
MCNUM restore_neutron = mccDivMonMono_restore_neutron;
MCNUM nx = mccDivMonMono_nx;
MCNUM ny = mccDivMonMono_ny;
MCNUM nz = mccDivMonMono_nz;
int nowritefile = mccDivMonMono_nowritefile;
#line 131 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
{
    
    multiline(5, (double)xmin, (double)ymin, 0.0,
                 (double)xmax, (double)ymin, 0.0,
                 (double)xmax, (double)ymax, 0.0,
                 (double)xmin, (double)ymax, 0.0,
                 (double)xmin, (double)ymin, 0.0);
}
#line 31034 "PSI_Focus.c"
}   /* End of DivMonMono=Divergence_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'focus_mono'. */
  SIG_MESSAGE("focus_mono (McDisplay)");
  printf("MCDISPLAY: component %s\n", "focus_mono");
#define mccompcurname  focus_mono
#define mccompcurtype  Arm
#define mccompcurindex 23
#line 40 "/usr/share/mcstas/2.5/optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 31059 "PSI_Focus.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'mono'. */
  SIG_MESSAGE("mono (McDisplay)");
  printf("MCDISPLAY: component %s\n", "mono");
#define mccompcurname  mono
#define mccompcurtype  Monochromator_2foc
#define mccompcurindex 24
#define mos_y mccmono_mos_y
#define mos_z mccmono_mos_z
#define mono_Q mccmono_mono_Q
#define SlabWidth mccmono_SlabWidth
#define SlabHeight mccmono_SlabHeight
#define rTable mccmono_rTable
{   /* Declarations of mono=Monochromator_2foc() SETTING parameters. */
char* reflect = mccmono_reflect;
MCNUM zwidth = mccmono_zwidth;
MCNUM yheight = mccmono_yheight;
MCNUM gap = mccmono_gap;
MCNUM NH = mccmono_NH;
MCNUM NV = mccmono_NV;
MCNUM mosaich = mccmono_mosaich;
MCNUM mosaicv = mccmono_mosaicv;
MCNUM r0 = mccmono_r0;
MCNUM Q = mccmono_Q;
MCNUM RV = mccmono_RV;
MCNUM RH = mccmono_RH;
MCNUM DM = mccmono_DM;
MCNUM mosaic = mccmono_mosaic;
MCNUM width = mccmono_width;
MCNUM height = mccmono_height;
MCNUM verbose = mccmono_verbose;
#line 270 "/usr/share/mcstas/2.5/contrib/Monochromator_2foc.comp"
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
#line 31129 "PSI_Focus.c"
}   /* End of mono=Monochromator_2foc() SETTING parameter declarations. */
#undef rTable
#undef SlabHeight
#undef SlabWidth
#undef mono_Q
#undef mos_z
#undef mos_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'a2'. */
  SIG_MESSAGE("a2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "a2");
#define mccompcurname  a2
#define mccompcurtype  Arm
#define mccompcurindex 25
#line 40 "/usr/share/mcstas/2.5/optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 31155 "PSI_Focus.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'FERMITOF_before'. */
  SIG_MESSAGE("FERMITOF_before (McDisplay)");
  printf("MCDISPLAY: component %s\n", "FERMITOF_before");
#define mccompcurname  FERMITOF_before
#define mccompcurtype  Monitor_nD
#define mccompcurindex 26
#define user1 mccFERMITOF_before_user1
#define user2 mccFERMITOF_before_user2
#define user3 mccFERMITOF_before_user3
#define DEFS mccFERMITOF_before_DEFS
#define Vars mccFERMITOF_before_Vars
#define detector mccFERMITOF_before_detector
#define offdata mccFERMITOF_before_offdata
{   /* Declarations of FERMITOF_before=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccFERMITOF_before_xwidth;
MCNUM yheight = mccFERMITOF_before_yheight;
MCNUM zdepth = mccFERMITOF_before_zdepth;
MCNUM xmin = mccFERMITOF_before_xmin;
MCNUM xmax = mccFERMITOF_before_xmax;
MCNUM ymin = mccFERMITOF_before_ymin;
MCNUM ymax = mccFERMITOF_before_ymax;
MCNUM zmin = mccFERMITOF_before_zmin;
MCNUM zmax = mccFERMITOF_before_zmax;
MCNUM bins = mccFERMITOF_before_bins;
MCNUM min = mccFERMITOF_before_min;
MCNUM max = mccFERMITOF_before_max;
MCNUM restore_neutron = mccFERMITOF_before_restore_neutron;
MCNUM radius = mccFERMITOF_before_radius;
char* options = mccFERMITOF_before_options;
char* filename = mccFERMITOF_before_filename;
char* geometry = mccFERMITOF_before_geometry;
char* username1 = mccFERMITOF_before_username1;
char* username2 = mccFERMITOF_before_username2;
char* username3 = mccFERMITOF_before_username3;
int nowritefile = mccFERMITOF_before_nowritefile;
#line 493 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
{
  if (geometry && strlen(geometry) && strcmp(geometry,"0") && strcmp(geometry, "NULL"))
  {
    off_display(offdata);
  } else {
    Monitor_nD_McDisplay(&DEFS, &Vars);
  }
}
#line 31204 "PSI_Focus.c"
}   /* End of FERMITOF_before=Monitor_nD() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'lambdaFermi'. */
  SIG_MESSAGE("lambdaFermi (McDisplay)");
  printf("MCDISPLAY: component %s\n", "lambdaFermi");
#define mccompcurname  lambdaFermi
#define mccompcurtype  L_monitor
#define mccompcurindex 27
#define nL mcclambdaFermi_nL
#define L_N mcclambdaFermi_L_N
#define L_p mcclambdaFermi_L_p
#define L_p2 mcclambdaFermi_L_p2
{   /* Declarations of lambdaFermi=L_monitor() SETTING parameters. */
char* filename = mcclambdaFermi_filename;
MCNUM xmin = mcclambdaFermi_xmin;
MCNUM xmax = mcclambdaFermi_xmax;
MCNUM ymin = mcclambdaFermi_ymin;
MCNUM ymax = mcclambdaFermi_ymax;
MCNUM xwidth = mcclambdaFermi_xwidth;
MCNUM yheight = mcclambdaFermi_yheight;
MCNUM Lmin = mcclambdaFermi_Lmin;
MCNUM Lmax = mcclambdaFermi_Lmax;
MCNUM restore_neutron = mcclambdaFermi_restore_neutron;
int nowritefile = mcclambdaFermi_nowritefile;
#line 120 "/usr/share/mcstas/2.5/monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 31248 "PSI_Focus.c"
}   /* End of lambdaFermi=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'EMON_Fermi'. */
  SIG_MESSAGE("EMON_Fermi (McDisplay)");
  printf("MCDISPLAY: component %s\n", "EMON_Fermi");
#define mccompcurname  EMON_Fermi
#define mccompcurtype  E_monitor
#define mccompcurindex 28
#define nE mccEMON_Fermi_nE
#define E_N mccEMON_Fermi_E_N
#define E_p mccEMON_Fermi_E_p
#define E_p2 mccEMON_Fermi_E_p2
#define S_p mccEMON_Fermi_S_p
#define S_pE mccEMON_Fermi_S_pE
#define S_pE2 mccEMON_Fermi_S_pE2
{   /* Declarations of EMON_Fermi=E_monitor() SETTING parameters. */
char* filename = mccEMON_Fermi_filename;
MCNUM xmin = mccEMON_Fermi_xmin;
MCNUM xmax = mccEMON_Fermi_xmax;
MCNUM ymin = mccEMON_Fermi_ymin;
MCNUM ymax = mccEMON_Fermi_ymax;
MCNUM xwidth = mccEMON_Fermi_xwidth;
MCNUM yheight = mccEMON_Fermi_yheight;
MCNUM Emin = mccEMON_Fermi_Emin;
MCNUM Emax = mccEMON_Fermi_Emax;
MCNUM restore_neutron = mccEMON_Fermi_restore_neutron;
int nowritefile = mccEMON_Fermi_nowritefile;
#line 132 "/usr/share/mcstas/2.5/monitors/E_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 31292 "PSI_Focus.c"
}   /* End of EMON_Fermi=E_monitor() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'DivMonfermi1'. */
  SIG_MESSAGE("DivMonfermi1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "DivMonfermi1");
#define mccompcurname  DivMonfermi1
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 29
#define nh mccDivMonfermi1_nh
#define nv mccDivMonfermi1_nv
#define Div_N mccDivMonfermi1_Div_N
#define Div_p mccDivMonfermi1_Div_p
#define Div_p2 mccDivMonfermi1_Div_p2
{   /* Declarations of DivMonfermi1=Divergence_monitor() SETTING parameters. */
char* filename = mccDivMonfermi1_filename;
MCNUM xmin = mccDivMonfermi1_xmin;
MCNUM xmax = mccDivMonfermi1_xmax;
MCNUM ymin = mccDivMonfermi1_ymin;
MCNUM ymax = mccDivMonfermi1_ymax;
MCNUM xwidth = mccDivMonfermi1_xwidth;
MCNUM yheight = mccDivMonfermi1_yheight;
MCNUM maxdiv_h = mccDivMonfermi1_maxdiv_h;
MCNUM maxdiv_v = mccDivMonfermi1_maxdiv_v;
MCNUM restore_neutron = mccDivMonfermi1_restore_neutron;
MCNUM nx = mccDivMonfermi1_nx;
MCNUM ny = mccDivMonfermi1_ny;
MCNUM nz = mccDivMonfermi1_nz;
int nowritefile = mccDivMonfermi1_nowritefile;
#line 131 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
{
    
    multiline(5, (double)xmin, (double)ymin, 0.0,
                 (double)xmax, (double)ymin, 0.0,
                 (double)xmax, (double)ymax, 0.0,
                 (double)xmin, (double)ymax, 0.0,
                 (double)xmin, (double)ymin, 0.0);
}
#line 31340 "PSI_Focus.c"
}   /* End of DivMonfermi1=Divergence_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PSD_Fermi1'. */
  SIG_MESSAGE("PSD_Fermi1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PSD_Fermi1");
#define mccompcurname  PSD_Fermi1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 30
#define nx mccPSD_Fermi1_nx
#define ny mccPSD_Fermi1_ny
#define PSD_N mccPSD_Fermi1_PSD_N
#define PSD_p mccPSD_Fermi1_PSD_p
#define PSD_p2 mccPSD_Fermi1_PSD_p2
{   /* Declarations of PSD_Fermi1=PSD_monitor() SETTING parameters. */
char* filename = mccPSD_Fermi1_filename;
MCNUM xmin = mccPSD_Fermi1_xmin;
MCNUM xmax = mccPSD_Fermi1_xmax;
MCNUM ymin = mccPSD_Fermi1_ymin;
MCNUM ymax = mccPSD_Fermi1_ymax;
MCNUM xwidth = mccPSD_Fermi1_xwidth;
MCNUM yheight = mccPSD_Fermi1_yheight;
MCNUM restore_neutron = mccPSD_Fermi1_restore_neutron;
int nowritefile = mccPSD_Fermi1_nowritefile;
#line 115 "/usr/share/mcstas/2.5/monitors/PSD_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 31381 "PSI_Focus.c"
}   /* End of PSD_Fermi1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'FoChopper'. */
  SIG_MESSAGE("FoChopper (McDisplay)");
  printf("MCDISPLAY: component %s\n", "FoChopper");
#define mccompcurname  FoChopper
#define mccompcurtype  FermiChopper
#define mccompcurindex 31
#define FCVars mccFoChopper_FCVars
{   /* Declarations of FoChopper=FermiChopper() SETTING parameters. */
MCNUM phase = mccFoChopper_phase;
MCNUM radius = mccFoChopper_radius;
MCNUM nu = mccFoChopper_nu;
MCNUM w = mccFoChopper_w;
MCNUM nslit = mccFoChopper_nslit;
MCNUM R0 = mccFoChopper_R0;
MCNUM Qc = mccFoChopper_Qc;
MCNUM alpha = mccFoChopper_alpha;
MCNUM m = mccFoChopper_m;
MCNUM W = mccFoChopper_W;
MCNUM length = mccFoChopper_length;
MCNUM eff = mccFoChopper_eff;
MCNUM zero_time = mccFoChopper_zero_time;
MCNUM xwidth = mccFoChopper_xwidth;
MCNUM verbose = mccFoChopper_verbose;
MCNUM yheight = mccFoChopper_yheight;
MCNUM curvature = mccFoChopper_curvature;
MCNUM delay = mccFoChopper_delay;
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
#line 31487 "PSI_Focus.c"
}   /* End of FoChopper=FermiChopper() SETTING parameter declarations. */
#undef FCVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PSD_Fermi2'. */
  SIG_MESSAGE("PSD_Fermi2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PSD_Fermi2");
#define mccompcurname  PSD_Fermi2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 32
#define nx mccPSD_Fermi2_nx
#define ny mccPSD_Fermi2_ny
#define PSD_N mccPSD_Fermi2_PSD_N
#define PSD_p mccPSD_Fermi2_PSD_p
#define PSD_p2 mccPSD_Fermi2_PSD_p2
{   /* Declarations of PSD_Fermi2=PSD_monitor() SETTING parameters. */
char* filename = mccPSD_Fermi2_filename;
MCNUM xmin = mccPSD_Fermi2_xmin;
MCNUM xmax = mccPSD_Fermi2_xmax;
MCNUM ymin = mccPSD_Fermi2_ymin;
MCNUM ymax = mccPSD_Fermi2_ymax;
MCNUM xwidth = mccPSD_Fermi2_xwidth;
MCNUM yheight = mccPSD_Fermi2_yheight;
MCNUM restore_neutron = mccPSD_Fermi2_restore_neutron;
int nowritefile = mccPSD_Fermi2_nowritefile;
#line 115 "/usr/share/mcstas/2.5/monitors/PSD_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 31524 "PSI_Focus.c"
}   /* End of PSD_Fermi2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'DivMonfermi2'. */
  SIG_MESSAGE("DivMonfermi2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "DivMonfermi2");
#define mccompcurname  DivMonfermi2
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 33
#define nh mccDivMonfermi2_nh
#define nv mccDivMonfermi2_nv
#define Div_N mccDivMonfermi2_Div_N
#define Div_p mccDivMonfermi2_Div_p
#define Div_p2 mccDivMonfermi2_Div_p2
{   /* Declarations of DivMonfermi2=Divergence_monitor() SETTING parameters. */
char* filename = mccDivMonfermi2_filename;
MCNUM xmin = mccDivMonfermi2_xmin;
MCNUM xmax = mccDivMonfermi2_xmax;
MCNUM ymin = mccDivMonfermi2_ymin;
MCNUM ymax = mccDivMonfermi2_ymax;
MCNUM xwidth = mccDivMonfermi2_xwidth;
MCNUM yheight = mccDivMonfermi2_yheight;
MCNUM maxdiv_h = mccDivMonfermi2_maxdiv_h;
MCNUM maxdiv_v = mccDivMonfermi2_maxdiv_v;
MCNUM restore_neutron = mccDivMonfermi2_restore_neutron;
MCNUM nx = mccDivMonfermi2_nx;
MCNUM ny = mccDivMonfermi2_ny;
MCNUM nz = mccDivMonfermi2_nz;
int nowritefile = mccDivMonfermi2_nowritefile;
#line 131 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
{
    
    multiline(5, (double)xmin, (double)ymin, 0.0,
                 (double)xmax, (double)ymin, 0.0,
                 (double)xmax, (double)ymax, 0.0,
                 (double)xmin, (double)ymax, 0.0,
                 (double)xmin, (double)ymin, 0.0);
}
#line 31570 "PSI_Focus.c"
}   /* End of DivMonfermi2=Divergence_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'FERMITOF1'. */
  SIG_MESSAGE("FERMITOF1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "FERMITOF1");
#define mccompcurname  FERMITOF1
#define mccompcurtype  Monitor_nD
#define mccompcurindex 34
#define user1 mccFERMITOF1_user1
#define user2 mccFERMITOF1_user2
#define user3 mccFERMITOF1_user3
#define DEFS mccFERMITOF1_DEFS
#define Vars mccFERMITOF1_Vars
#define detector mccFERMITOF1_detector
#define offdata mccFERMITOF1_offdata
{   /* Declarations of FERMITOF1=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccFERMITOF1_xwidth;
MCNUM yheight = mccFERMITOF1_yheight;
MCNUM zdepth = mccFERMITOF1_zdepth;
MCNUM xmin = mccFERMITOF1_xmin;
MCNUM xmax = mccFERMITOF1_xmax;
MCNUM ymin = mccFERMITOF1_ymin;
MCNUM ymax = mccFERMITOF1_ymax;
MCNUM zmin = mccFERMITOF1_zmin;
MCNUM zmax = mccFERMITOF1_zmax;
MCNUM bins = mccFERMITOF1_bins;
MCNUM min = mccFERMITOF1_min;
MCNUM max = mccFERMITOF1_max;
MCNUM restore_neutron = mccFERMITOF1_restore_neutron;
MCNUM radius = mccFERMITOF1_radius;
char* options = mccFERMITOF1_options;
char* filename = mccFERMITOF1_filename;
char* geometry = mccFERMITOF1_geometry;
char* username1 = mccFERMITOF1_username1;
char* username2 = mccFERMITOF1_username2;
char* username3 = mccFERMITOF1_username3;
int nowritefile = mccFERMITOF1_nowritefile;
#line 493 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
{
  if (geometry && strlen(geometry) && strcmp(geometry,"0") && strcmp(geometry, "NULL"))
  {
    off_display(offdata);
  } else {
    Monitor_nD_McDisplay(&DEFS, &Vars);
  }
}
#line 31625 "PSI_Focus.c"
}   /* End of FERMITOF1=Monitor_nD() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'SAMPLE_SLIT'. */
  SIG_MESSAGE("SAMPLE_SLIT (McDisplay)");
  printf("MCDISPLAY: component %s\n", "SAMPLE_SLIT");
#define mccompcurname  SAMPLE_SLIT
#define mccompcurtype  Slit
#define mccompcurindex 35
{   /* Declarations of SAMPLE_SLIT=Slit() SETTING parameters. */
MCNUM xmin = mccSAMPLE_SLIT_xmin;
MCNUM xmax = mccSAMPLE_SLIT_xmax;
MCNUM ymin = mccSAMPLE_SLIT_ymin;
MCNUM ymax = mccSAMPLE_SLIT_ymax;
MCNUM radius = mccSAMPLE_SLIT_radius;
MCNUM xwidth = mccSAMPLE_SLIT_xwidth;
MCNUM yheight = mccSAMPLE_SLIT_yheight;
#line 81 "/usr/share/mcstas/2.5/optics/Slit.comp"
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
#line 31675 "PSI_Focus.c"
}   /* End of SAMPLE_SLIT=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'FERMITOF2'. */
  SIG_MESSAGE("FERMITOF2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "FERMITOF2");
#define mccompcurname  FERMITOF2
#define mccompcurtype  Monitor_nD
#define mccompcurindex 36
#define user1 mccFERMITOF2_user1
#define user2 mccFERMITOF2_user2
#define user3 mccFERMITOF2_user3
#define DEFS mccFERMITOF2_DEFS
#define Vars mccFERMITOF2_Vars
#define detector mccFERMITOF2_detector
#define offdata mccFERMITOF2_offdata
{   /* Declarations of FERMITOF2=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccFERMITOF2_xwidth;
MCNUM yheight = mccFERMITOF2_yheight;
MCNUM zdepth = mccFERMITOF2_zdepth;
MCNUM xmin = mccFERMITOF2_xmin;
MCNUM xmax = mccFERMITOF2_xmax;
MCNUM ymin = mccFERMITOF2_ymin;
MCNUM ymax = mccFERMITOF2_ymax;
MCNUM zmin = mccFERMITOF2_zmin;
MCNUM zmax = mccFERMITOF2_zmax;
MCNUM bins = mccFERMITOF2_bins;
MCNUM min = mccFERMITOF2_min;
MCNUM max = mccFERMITOF2_max;
MCNUM restore_neutron = mccFERMITOF2_restore_neutron;
MCNUM radius = mccFERMITOF2_radius;
char* options = mccFERMITOF2_options;
char* filename = mccFERMITOF2_filename;
char* geometry = mccFERMITOF2_geometry;
char* username1 = mccFERMITOF2_username1;
char* username2 = mccFERMITOF2_username2;
char* username3 = mccFERMITOF2_username3;
int nowritefile = mccFERMITOF2_nowritefile;
#line 493 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
{
  if (geometry && strlen(geometry) && strcmp(geometry,"0") && strcmp(geometry, "NULL"))
  {
    off_display(offdata);
  } else {
    Monitor_nD_McDisplay(&DEFS, &Vars);
  }
}
#line 31725 "PSI_Focus.c"
}   /* End of FERMITOF2=Monitor_nD() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'PSD_SAMPLE'. */
  SIG_MESSAGE("PSD_SAMPLE (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PSD_SAMPLE");
#define mccompcurname  PSD_SAMPLE
#define mccompcurtype  PSD_monitor
#define mccompcurindex 37
#define nx mccPSD_SAMPLE_nx
#define ny mccPSD_SAMPLE_ny
#define PSD_N mccPSD_SAMPLE_PSD_N
#define PSD_p mccPSD_SAMPLE_PSD_p
#define PSD_p2 mccPSD_SAMPLE_PSD_p2
{   /* Declarations of PSD_SAMPLE=PSD_monitor() SETTING parameters. */
char* filename = mccPSD_SAMPLE_filename;
MCNUM xmin = mccPSD_SAMPLE_xmin;
MCNUM xmax = mccPSD_SAMPLE_xmax;
MCNUM ymin = mccPSD_SAMPLE_ymin;
MCNUM ymax = mccPSD_SAMPLE_ymax;
MCNUM xwidth = mccPSD_SAMPLE_xwidth;
MCNUM yheight = mccPSD_SAMPLE_yheight;
MCNUM restore_neutron = mccPSD_SAMPLE_restore_neutron;
int nowritefile = mccPSD_SAMPLE_nowritefile;
#line 115 "/usr/share/mcstas/2.5/monitors/PSD_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 31768 "PSI_Focus.c"
}   /* End of PSD_SAMPLE=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'DivMon_Sample'. */
  SIG_MESSAGE("DivMon_Sample (McDisplay)");
  printf("MCDISPLAY: component %s\n", "DivMon_Sample");
#define mccompcurname  DivMon_Sample
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 38
#define nh mccDivMon_Sample_nh
#define nv mccDivMon_Sample_nv
#define Div_N mccDivMon_Sample_Div_N
#define Div_p mccDivMon_Sample_Div_p
#define Div_p2 mccDivMon_Sample_Div_p2
{   /* Declarations of DivMon_Sample=Divergence_monitor() SETTING parameters. */
char* filename = mccDivMon_Sample_filename;
MCNUM xmin = mccDivMon_Sample_xmin;
MCNUM xmax = mccDivMon_Sample_xmax;
MCNUM ymin = mccDivMon_Sample_ymin;
MCNUM ymax = mccDivMon_Sample_ymax;
MCNUM xwidth = mccDivMon_Sample_xwidth;
MCNUM yheight = mccDivMon_Sample_yheight;
MCNUM maxdiv_h = mccDivMon_Sample_maxdiv_h;
MCNUM maxdiv_v = mccDivMon_Sample_maxdiv_v;
MCNUM restore_neutron = mccDivMon_Sample_restore_neutron;
MCNUM nx = mccDivMon_Sample_nx;
MCNUM ny = mccDivMon_Sample_ny;
MCNUM nz = mccDivMon_Sample_nz;
int nowritefile = mccDivMon_Sample_nowritefile;
#line 131 "/usr/share/mcstas/2.5/monitors/Divergence_monitor.comp"
{
    
    multiline(5, (double)xmin, (double)ymin, 0.0,
                 (double)xmax, (double)ymin, 0.0,
                 (double)xmax, (double)ymax, 0.0,
                 (double)xmin, (double)ymax, 0.0,
                 (double)xmin, (double)ymin, 0.0);
}
#line 31814 "PSI_Focus.c"
}   /* End of DivMon_Sample=Divergence_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'EMON_SAMPLE'. */
  SIG_MESSAGE("EMON_SAMPLE (McDisplay)");
  printf("MCDISPLAY: component %s\n", "EMON_SAMPLE");
#define mccompcurname  EMON_SAMPLE
#define mccompcurtype  E_monitor
#define mccompcurindex 39
#define nE mccEMON_SAMPLE_nE
#define E_N mccEMON_SAMPLE_E_N
#define E_p mccEMON_SAMPLE_E_p
#define E_p2 mccEMON_SAMPLE_E_p2
#define S_p mccEMON_SAMPLE_S_p
#define S_pE mccEMON_SAMPLE_S_pE
#define S_pE2 mccEMON_SAMPLE_S_pE2
{   /* Declarations of EMON_SAMPLE=E_monitor() SETTING parameters. */
char* filename = mccEMON_SAMPLE_filename;
MCNUM xmin = mccEMON_SAMPLE_xmin;
MCNUM xmax = mccEMON_SAMPLE_xmax;
MCNUM ymin = mccEMON_SAMPLE_ymin;
MCNUM ymax = mccEMON_SAMPLE_ymax;
MCNUM xwidth = mccEMON_SAMPLE_xwidth;
MCNUM yheight = mccEMON_SAMPLE_yheight;
MCNUM Emin = mccEMON_SAMPLE_Emin;
MCNUM Emax = mccEMON_SAMPLE_Emax;
MCNUM restore_neutron = mccEMON_SAMPLE_restore_neutron;
int nowritefile = mccEMON_SAMPLE_nowritefile;
#line 132 "/usr/share/mcstas/2.5/monitors/E_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 31859 "PSI_Focus.c"
}   /* End of EMON_SAMPLE=E_monitor() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'a3'. */
  SIG_MESSAGE("a3 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "a3");
#define mccompcurname  a3
#define mccompcurtype  Arm
#define mccompcurindex 40
#line 40 "/usr/share/mcstas/2.5/optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 31886 "PSI_Focus.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Sample'. */
  SIG_MESSAGE("Sample (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Sample");
#define mccompcurname  Sample
#define mccompcurtype  V_sample
#define mccompcurindex 41
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
#line 31976 "PSI_Focus.c"
}   /* End of Sample=V_sample() SETTING parameter declarations. */
#undef VarsV
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'TOF_Det'. */
  SIG_MESSAGE("TOF_Det (McDisplay)");
  printf("MCDISPLAY: component %s\n", "TOF_Det");
#define mccompcurname  TOF_Det
#define mccompcurtype  Monitor_nD
#define mccompcurindex 42
#define user1 mccTOF_Det_user1
#define user2 mccTOF_Det_user2
#define user3 mccTOF_Det_user3
#define DEFS mccTOF_Det_DEFS
#define Vars mccTOF_Det_Vars
#define detector mccTOF_Det_detector
#define offdata mccTOF_Det_offdata
{   /* Declarations of TOF_Det=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccTOF_Det_xwidth;
MCNUM yheight = mccTOF_Det_yheight;
MCNUM zdepth = mccTOF_Det_zdepth;
MCNUM xmin = mccTOF_Det_xmin;
MCNUM xmax = mccTOF_Det_xmax;
MCNUM ymin = mccTOF_Det_ymin;
MCNUM ymax = mccTOF_Det_ymax;
MCNUM zmin = mccTOF_Det_zmin;
MCNUM zmax = mccTOF_Det_zmax;
MCNUM bins = mccTOF_Det_bins;
MCNUM min = mccTOF_Det_min;
MCNUM max = mccTOF_Det_max;
MCNUM restore_neutron = mccTOF_Det_restore_neutron;
MCNUM radius = mccTOF_Det_radius;
char* options = mccTOF_Det_options;
char* filename = mccTOF_Det_filename;
char* geometry = mccTOF_Det_geometry;
char* username1 = mccTOF_Det_username1;
char* username2 = mccTOF_Det_username2;
char* username3 = mccTOF_Det_username3;
int nowritefile = mccTOF_Det_nowritefile;
#line 493 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
{
  if (geometry && strlen(geometry) && strcmp(geometry,"0") && strcmp(geometry, "NULL"))
  {
    off_display(offdata);
  } else {
    Monitor_nD_McDisplay(&DEFS, &Vars);
  }
}
#line 32027 "PSI_Focus.c"
}   /* End of TOF_Det=Monitor_nD() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'FoDet'. */
  SIG_MESSAGE("FoDet (McDisplay)");
  printf("MCDISPLAY: component %s\n", "FoDet");
#define mccompcurname  FoDet
#define mccompcurtype  Monitor_nD
#define mccompcurindex 43
#define user1 mccFoDet_user1
#define user2 mccFoDet_user2
#define user3 mccFoDet_user3
#define DEFS mccFoDet_DEFS
#define Vars mccFoDet_Vars
#define detector mccFoDet_detector
#define offdata mccFoDet_offdata
{   /* Declarations of FoDet=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccFoDet_xwidth;
MCNUM yheight = mccFoDet_yheight;
MCNUM zdepth = mccFoDet_zdepth;
MCNUM xmin = mccFoDet_xmin;
MCNUM xmax = mccFoDet_xmax;
MCNUM ymin = mccFoDet_ymin;
MCNUM ymax = mccFoDet_ymax;
MCNUM zmin = mccFoDet_zmin;
MCNUM zmax = mccFoDet_zmax;
MCNUM bins = mccFoDet_bins;
MCNUM min = mccFoDet_min;
MCNUM max = mccFoDet_max;
MCNUM restore_neutron = mccFoDet_restore_neutron;
MCNUM radius = mccFoDet_radius;
char* options = mccFoDet_options;
char* filename = mccFoDet_filename;
char* geometry = mccFoDet_geometry;
char* username1 = mccFoDet_username1;
char* username2 = mccFoDet_username2;
char* username3 = mccFoDet_username3;
int nowritefile = mccFoDet_nowritefile;
#line 493 "/usr/share/mcstas/2.5/monitors/Monitor_nD.comp"
{
  if (geometry && strlen(geometry) && strcmp(geometry,"0") && strcmp(geometry, "NULL"))
  {
    off_display(offdata);
  } else {
    Monitor_nD_McDisplay(&DEFS, &Vars);
  }
}
#line 32084 "PSI_Focus.c"
}   /* End of FoDet=Monitor_nD() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'EMON_DET'. */
  SIG_MESSAGE("EMON_DET (McDisplay)");
  printf("MCDISPLAY: component %s\n", "EMON_DET");
#define mccompcurname  EMON_DET
#define mccompcurtype  E_monitor
#define mccompcurindex 44
#define nE mccEMON_DET_nE
#define E_N mccEMON_DET_E_N
#define E_p mccEMON_DET_E_p
#define E_p2 mccEMON_DET_E_p2
#define S_p mccEMON_DET_S_p
#define S_pE mccEMON_DET_S_pE
#define S_pE2 mccEMON_DET_S_pE2
{   /* Declarations of EMON_DET=E_monitor() SETTING parameters. */
char* filename = mccEMON_DET_filename;
MCNUM xmin = mccEMON_DET_xmin;
MCNUM xmax = mccEMON_DET_xmax;
MCNUM ymin = mccEMON_DET_ymin;
MCNUM ymax = mccEMON_DET_ymax;
MCNUM xwidth = mccEMON_DET_xwidth;
MCNUM yheight = mccEMON_DET_yheight;
MCNUM Emin = mccEMON_DET_Emin;
MCNUM Emax = mccEMON_DET_Emax;
MCNUM restore_neutron = mccEMON_DET_restore_neutron;
int nowritefile = mccEMON_DET_nowritefile;
#line 132 "/usr/share/mcstas/2.5/monitors/E_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 32131 "PSI_Focus.c"
}   /* End of EMON_DET=E_monitor() SETTING parameter declarations. */
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
/* end of generated C code PSI_Focus.c */
