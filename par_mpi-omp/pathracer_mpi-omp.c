#define _POSIX_C_SOURCE 199506L /* Pour rand_r */
#define _XOPEN_SOURCE
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>  /* pour mkdir    */
#include <unistd.h>    /* pour getuid   */
#include <sys/types.h> /* pour getpwuid */
#include <pwd.h>       /* pour getpwuid */
#include <string.h>	   /* pour memcpy */
#include "mpi.h"
#include <omp.h>
#include "list_res.c"

enum Refl_t {DIFF, SPEC, REFR};   /* types de matériaux (DIFFuse, SPECular, REFRactive) */

struct Sphere {
  double radius;
  double position[3];
  double emission[3];     /* couleur émise (=source de lumière) */
  double color[3];        /* couleur de l'objet RGB (diffusion, refraction, ...) */
  enum Refl_t refl;       /* type de reflection */
  double max_reflexivity;
};

static const int KILL_DEPTH = 7;
static const int SPLIT_DEPTH = 4;

/* la scène est composée uniquement de spheres */
struct Sphere spheres[] = {
  // radius position,                         emission,     color,              material
  {1e5,  { 1e5+1,  40.8,       81.6},      {},           {.75,  .25,  .25},  DIFF, -1}, // Left
  {1e5,  {-1e5+99, 40.8,       81.6},      {},           {.25,  .25,  .75},  DIFF, -1}, // Right
  {1e5,  {50,      40.8,       1e5},       {},           {.75,  .75,  .75},  DIFF, -1}, // Back
  {1e5,  {50,      40.8,      -1e5 + 170}, {},           {},                 DIFF, -1}, // Front
  {1e5,  {50,      1e5,        81.6},      {},           {0.75, .75,  .75},  DIFF, -1}, // Bottom
  {1e5,  {50,     -1e5 + 81.6, 81.6},      {},           {0.75, .75,  .75},  DIFF, -1}, // Top
  {16.5, {40,      16.5,       47},        {},           {.999, .999, .999}, SPEC, -1}, // Mirror
  {16.5, {73,      46.5,       88},        {},           {.999, .999, .999}, REFR, -1}, // Glass
  {10,   {15,      45,         112},       {},           {.999, .999, .999}, DIFF, -1}, // white ball
  {15,   {16,      16,         130},       {},           {.999, .999, 0},    REFR, -1}, // big yellow glass
  {7.5,  {40,      8,          120},        {},           {.999, .999, 0   }, REFR, -1}, // small yellow glass middle
  {8.5,  {60,      9,          110},        {},           {.999, .999, 0   }, REFR, -1}, // small yellow glass right
  {10,   {80,      12,         92},        {},           {0, .999, 0},       DIFF, -1}, // green ball
  {600,  {50,      681.33,     81.6},      {12, 12, 12}, {},                 DIFF, -1},  // Light
  {5,    {50,      75,         81.6},      {},           {0, .682, .999}, DIFF, -1}, // occlusion, mirror
};


/********** micro BLAS LEVEL-1 + quelques fonctions non-standard **************/
static inline void copy(const double *x, double *y)
{
  for (int i = 0; i < 3; i++)
    y[i] = x[i];
}

static inline void zero(double *x)
{
  for (int i = 0; i < 3; i++)
    x[i] = 0;
}

static inline void axpy(double alpha, const double *x, double *y)
{
  for (int i = 0; i < 3; i++)
    y[i] += alpha * x[i];
}

static inline void scal(double alpha, double *x)
{
  for (int i = 0; i < 3; i++)
    x[i] *= alpha;
}

static inline double dot(const double *a, const double *b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static inline double nrm2(const double *a)
{
  return sqrt(dot(a, a));
}

/********* fonction non-standard *************/
static inline void mul(const double *x, const double *y, double *z)
{
  for (int i = 0; i < 3; i++)
    z[i] = x[i] * y[i];
}

static inline void normalize(double *x)
{
  scal(1 / nrm2(x), x);
}

/* produit vectoriel */
static inline void cross(const double *a, const double *b, double *c)
{
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

/****** tronque *************/
static inline void clamp(double *x)
{
  int i;
  for (i = 0; i < 3; i++) {
    if (x[i] < 0)
      x[i] = 0;
    if (x[i] > 1)
      x[i] = 1;
  }
}

/******************************* calcul des intersections rayon / sphere *************************************/

/**
 * Calcul l'interssection entre un rayon et une sphere donnée
 * returns distance, 0 if nohit
 */
double sphere_intersect(const struct Sphere *s, const double *ray_origin, const double *ray_direction)
{
  double op[3];
  // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
  copy(s->position, op);
  axpy(-1, ray_origin, op);
  double eps = 1e-4;
  double b = dot(op, ray_direction);
  double discriminant = b * b - dot(op, op) + s->radius * s->radius;
  if (discriminant < 0)
    return 0;   /* pas d'intersection */
  else
    discriminant = sqrt(discriminant);
  /* détermine la plus petite solution positive (i.e. point d'intersection le plus proche, mais devant nous) */
  double t = b - discriminant;
  if (t > eps) {
    return t;
  } else {
    t = b + discriminant;
    if (t > eps)
      return t;
    else
      return 0;  /* cas bizarre, racine double, etc. */
  }
}

/**
 * détermine si le rayon intersecte l'une des sphere; si oui renvoie true et fixe t, id
 * (Calcule l'interssection entre un rayon et TOUS les spheres)
 */
bool intersect(const double *ray_origin, const double *ray_direction, double *t, int *id)
{
  int n = sizeof(spheres) / sizeof(struct Sphere);
  double inf = 1e20;
  *t = inf;
  // Parallélisation impossible ici car i+1 utilise *t de i
  for (int i = 0; i < n; i++) {
    double d = sphere_intersect(&spheres[i], ray_origin, ray_direction);
    if ((d > 0) && (d < *t)) {
      *t = d;
      *id = i;
    }
  }
  return *t < inf;
}

/**
 * calcule (dans out) la lumiance reçue par la camera sur le rayon donné
 */
void radiance(const double *ray_origin, const double *ray_direction, int depth, unsigned int *PRNG_state, double *out)
{
  int id = 0;                             // id de la sphère intersectée par le rayon
  double t;                               // distance à l'intersection
  if (!intersect(ray_origin, ray_direction, &t, &id)) {
    zero(out);    // if miss, return black
    return;
  }
  const struct Sphere *obj = &spheres[id];

  /* point d'intersection du rayon et de la sphère */
  double x[3];
  copy(ray_origin, x);
  axpy(t, ray_direction, x);

  /* vecteur normal à la sphere, au point d'intersection */
  double n[3];
  copy(x, n);
  axpy(-1, obj->position, n);
  normalize(n);

  /* vecteur normal, orienté dans le sens opposé au rayon
     (vers l'extérieur si le rayon entre, vers l'intérieur s'il sort) */
  double nl[3];
  copy(n, nl);
  if (dot(n, ray_direction) > 0)
    scal(-1, nl);

  /* couleur de la sphere */
  double f[3];
  copy(obj->color, f);
  double p = obj->max_reflexivity;

  /* processus aléatoire : au-delà d'une certaine profondeur,
     décide aléatoirement d'arrêter la récusion. Plus l'objet est
     clair, plus le processus a de chance de continuer. */
  depth++;
  if (depth > KILL_DEPTH) {
    if ((double)rand_r(PRNG_state)/RAND_MAX < p) {
      scal(1 / p, f);
    } else {
      copy(obj->emission, out);
      return;
    }
  }

  /* Cas de la réflection DIFFuse (= non-brillante).
     On récupère la luminance en provenance de l'ensemble de l'univers.
     Pour cela : (processus de monte-carlo) on choisit une direction
     aléatoire dans un certain cone, et on récupère la luminance en
     provenance de cette direction. */
  if (obj->refl == DIFF) {
    double r1 = 2 * M_PI * (double)rand_r(PRNG_state)/RAND_MAX;  /* angle aléatoire */
    double r2 = (double)rand_r(PRNG_state)/RAND_MAX;             /* distance au centre aléatoire */
    double r2s = sqrt(r2);

    double w[3];   /* vecteur normal */
    copy(nl, w);

    double u[3];   /* u est orthogonal à w */
    double uw[3] = {0, 0, 0};
    if (fabs(w[0]) > .1)
      uw[1] = 1;
    else
      uw[0] = 1;
    cross(uw, w, u);
    normalize(u);

    double v[3];   /* v est orthogonal à u et w */
    cross(w, u, v);

    double d[3];   /* d est le vecteur incident aléatoire, selon la bonne distribution */
    zero(d);
    axpy(cos(r1) * r2s, u, d);
    axpy(sin(r1) * r2s, v, d);
    axpy(sqrt(1 - r2), w, d);
    normalize(d);

    /* calcule récursivement la luminance du rayon incident */
    double rec[3];
    radiance(x, d, depth, PRNG_state, rec);

    /* pondère par la couleur de la sphère, prend en compte l'emissivité */
    mul(f, rec, out);
    axpy(1, obj->emission, out);
    return;
  }

  /* dans les deux autres cas (réflection parfaite / refraction), on considère le rayon
     réfléchi par la spère */

  double reflected_dir[3];
  copy(ray_direction, reflected_dir);
  axpy(-2 * dot(n, ray_direction), n, reflected_dir);

  /* cas de la reflection SPEculaire parfaire (==mirroir) */
  if (obj->refl == SPEC) {
    double rec[3];
    /* calcule récursivement la luminance du rayon réflechi */
    radiance(x, reflected_dir, depth, PRNG_state, rec);
    /* pondère par la couleur de la sphère, prend en compte l'emissivité */
    mul(f, rec, out);
    axpy(1, obj->emission, out);
    return;
  }

  /* cas des surfaces diélectriques (==verre). Combinaison de réflection et de réfraction. */
  bool into = dot(n, nl) > 0;      /* vient-il de l'extérieur ? */
  double nc = 1;                   /* indice de réfraction de l'air */
  double nt = 1.5;                 /* indice de réfraction du verre */
  double nnt = into ? (nc / nt) : (nt / nc);
  double ddn = dot(ray_direction, nl);

  /* si le rayon essaye de sortir de l'objet en verre avec un angle incident trop faible,
     il rebondit entièrement */
  double cos2t = 1 - nnt * nnt * (1 - ddn * ddn);
  if (cos2t < 0) {
    double rec[3];
    /* calcule seulement le rayon réfléchi */
    radiance(x, reflected_dir, depth, PRNG_state, rec);
    mul(f, rec, out);
    axpy(1, obj->emission, out);
    return;
  }

  /* calcule la direction du rayon réfracté */
  double tdir[3];
  zero(tdir);
  axpy(nnt, ray_direction, tdir);
  axpy(-(into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)), n, tdir);

  /* calcul de la réflectance (==fraction de la lumière réfléchie) */
  double a = nt - nc;
  double b = nt + nc;
  double R0 = a * a / (b * b);
  double c = 1 - (into ? -ddn : dot(tdir, n));
  double Re = R0 + (1 - R0) * c * c * c * c * c;   /* réflectance */
  double Tr = 1 - Re;                              /* transmittance */

  /* au-dela d'une certaine profondeur, on choisit aléatoirement si
     on calcule le rayon réfléchi ou bien le rayon réfracté. En dessous du
     seuil, on calcule les deux. */
  double rec[3];
  if (depth > SPLIT_DEPTH) {
    double P = .25 + .5 * Re;             /* probabilité de réflection */
    if ((double)rand_r(PRNG_state)/RAND_MAX < P) {
      radiance(x, reflected_dir, depth, PRNG_state, rec);
      double RP = Re / P;
      scal(RP, rec);
    } else {
      radiance(x, tdir, depth, PRNG_state, rec);
      double TP = Tr / (1 - P);
      scal(TP, rec);
    }
  } else {
    double rec_re[3], rec_tr[3];
    radiance(x, reflected_dir, depth, PRNG_state, rec_re);
    radiance(x, tdir, depth, PRNG_state, rec_tr);
    zero(rec);
    axpy(Re, rec_re, rec);
    axpy(Tr, rec_tr, rec);
  }
  /* pondère, prend en compte la luminance */
  mul(f, rec, out);
  axpy(1, obj->emission, out);
  return;
}

double wtime()
{
  struct timeval ts;
  gettimeofday(&ts, NULL);
  return (double)ts.tv_sec + ts.tv_usec / 1E6;
}

int toInt(double x)
{
  return pow(x, 1 / 2.2) * 255 + .5;   /* gamma correction = 2.2 */
}

/**
   Nombre de bloc restant au process
*/
int nblocs_left(int* blocs, int nblocs, int ind) {
  int res=0;
  for(int i=ind; i<nblocs && blocs[i]!=(-1); i++, res++) {}
  return res;
}

void print_tab(int *tab, int l) {
  fprintf(stderr, "{ ");
#pragma omp parallel for ordered
  for(int i=0; i<l; i++) fprintf(stderr, "%d ", tab[i]);
  fprintf(stderr, "}\n");
}


#define TAG_DATA 73 /* Tag pour recevoir les données de l'image */
#define TAG_REQ 37 /* Tag pour une demande d'envoi de nouveaux blocs */
#define TAG_END 0 /* Tag pour annoncer la fin de ses blocs (reponse négative à TAG_REQ) */

int main(int argc, char **argv)
{
  /* Moyen cas test */
  int w = 2000;
  int h = 1500;
  int samples = 500;

  int cur_h; /* hauteur courant dans l'image intégrale (valeur de i dans l'image intégrale)*/
  int nlines; /* Nombre de ligne que contiendra chaque bloc de travail */
  int flag=0, flag_end=0;
  bool end_sent=false;

  double *image; /* image intégrale du resultat finale */
  int *blocs; /* Tableau de numeros des blocs sur lesquels devra travailler le process */
  int ind=0; /* indice du bloc sur lequel le process travaille dans le tableau */
  int nblocs; /* Nombre de blocs */
  int buf=1, rest;

  list_res_t* list_res=allocate_list(); /* liste des resultats des images traitées localement */

  /* Variables mpi */
  int rank, np, tag, root=0; // rang, nbre processus, tag, processus root
  MPI_Status status;
  MPI_Request request;

  /* Chronometrage */
  double debut, fin;

  /* debut du chronometrage */
  debut = wtime();

  /* Initialisation MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  nlines=5; /* Valeur par défaut */

  if (argc > 1) {
    int im=atoi(argv[1]);
    if(im==1) {
      /* Petit cas test (small, quick and dirty): */
      w = 320;
      h = 200;
      samples = 200;
    } else if(im==2) {
      /* Gros cas test (big, slow and pretty): */
      w = 3840;
      h = 2160;
      samples = 5000;
    }
    if (argc > 2) {
      samples = atoi(argv[2]) / 4;
    } if(argc > 3) {
      nlines = atoi(argv[3]);
      if (argc == 5) root = atoi(argv[4]);
      else if(argc > 5 && rank==root) fprintf(stderr, "mpirun ... pathracer_mpi-omp <dim_image : 1 petit cas | 2 : gros cas> <samples : int> <nlines : int> <root : int>\n");
    }
  }

  static const double CST = 0.5135;  /* ceci défini l'angle de vue */
  double camera_position[3] = {50, 52, 295.6};
  double camera_direction[3] = {0, -0.042612, -1};
  normalize(camera_direction);

  if(rank==root) fprintf(stderr, "hauteur=%d, largeur=%d, nombre de lignes par bloc=%d, nombre de processus=%d\n", h, w, nlines, np);

  /* hypothèses de départ */
  if(np<2 || h%nlines != 0 || np*nlines > h) {
    if(rank==root) fprintf(stderr, "Nombre de process invalide ou nombre de ligne invalide. Nombre de ligne doit diviser la hauteur\n");
    exit(1);
  }

  /* incréments pour passer d'un pixel à l'autre */
  double cx[3] = {w * CST / h, 0, 0};
  double cy[3];
  cross(cx, camera_direction, cy);  /* cy est orthogonal à cx ET à la direction dans laquelle regarde la caméra */
  normalize(cy);
  scal(CST, cy);

  nblocs=floor(h/(np*nlines));
  rest=(h/nlines) % np;
  if(rank < rest) nblocs++; /* Les processus de rang < rest ont un bloc en plus */

  /* précalcule la norme infinie des couleurs */
  /* Est-ce nécessaire de paralléliser ? Nbre de sphere négligeable
   * Seul le process root le fera */
  if(rank==root) {
    int n = sizeof(spheres) / sizeof(struct Sphere);
    int i;
#pragma omp parallel for schedule(runtime)
    for (i = 0; i < n; i++) {
      double *f = spheres[i].color;
      if ((f[0] > f[1]) && (f[0] > f[2]))
	spheres[i].max_reflexivity = f[0];
      else {
	if (f[1] > f[2])
	  spheres[i].max_reflexivity = f[1];
	else
	  spheres[i].max_reflexivity = f[2];
      }
    }
  }

  /* Les autres process attendent que root finissent de faire les calculs préalabes */
  MPI_Barrier(MPI_COMM_WORLD);

  /* Allocation du tableau des numeros des blocs */
  blocs=(int*)malloc(sizeof(int)*nblocs);

  /* Répartition préalable (statique) des tâches */
  {
    int i;
#pragma omp parallel for schedule(runtime)
    for(i=0; i<nblocs; i++){
      blocs[i]=np * i + rank;
    }
  }
  fprintf(stderr, "\n[%d] : nblocs=%d, Répartition initiale des blocs de tâche : ", rank, nblocs);
  print_tab(blocs, nblocs);

  /* root alloue l'image intégrale pour le résultat final
   * root travaillera directement sur l'image finale
   * */
  if(rank==root) {
    image = malloc(3 * w * h * sizeof(*image));
    if (image == NULL) {
      perror("Impossible d'allouer l'image");
      exit(1);
    }
  } else image=NULL;

  /* boucle principale */
  /* Tant qu'il reste des blocs de travail (le travail n'est pas totalement achevé) */
  while(!flag_end) {

    /* Tant qu'on a encore des blocs de travail */
    while(blocs[ind] != -1) {

      /* Les process non root allouent leur tampon image */
      image_res_t* ima_res;
      if(rank != root) ima_res = create_image_res(blocs[ind], w, nlines);
      else ima_res = create_image_res(blocs[ind], 0, 0); /* pour la répartition finale */
      blocs[ind]=-1;
      ind=(ind+1)%nblocs; /* incrémentation vers le prochain bloc de travail */

      for (int i = 0; i < nlines; i++) {
	cur_h=i+h-ima_res->nbloc*nlines-nlines;
	unsigned int PRNG_state=time(NULL);
	for (unsigned short j = 0; j < w; j++) {
	  /* calcule la luminance d'un pixel, avec sur-échantillonnage 2x2 */
	  double pixel_radiance[3] = {0, 0, 0};
	  for (int sub_i = 0; sub_i < 2; sub_i++) {
	    for (int sub_j = 0; sub_j < 2; sub_j++) {
	      double subpixel_radiance[3] = {0};
	      /* simulation de monte-carlo : on effectue plein de lancers de rayons et on moyenne */
	      //PARALLELISABLE omp
#pragma omp parallel for schedule(runtime)
	      for (int s = 0; s < samples; s++) {
		/* tire un rayon aléatoire dans une zone de la caméra qui correspond à peu près au pixel à calculer */
		double r1 = 2 * (double)rand_r(&PRNG_state)/RAND_MAX;
		double dx = (r1 < 1) ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
		double r2 = 2 * (double)rand_r(&PRNG_state)/RAND_MAX;
		double dy = (r2 < 1) ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
		double ray_direction[3];
		copy(camera_direction, ray_direction);
		axpy(((sub_i + .5 + dy) / 2 + cur_h) / h - .5, cy, ray_direction); /**/
		axpy(((sub_j + .5 + dx) / 2 + j) / w - .5, cx, ray_direction);
		normalize(ray_direction);

		double ray_origin[3];
		copy(camera_position, ray_origin);
		axpy(140, ray_direction, ray_origin);

		/* estime la lumiance qui arrive sur la caméra par ce rayon */
		double sample_radiance[3];
		radiance(ray_origin, ray_direction, 0, &PRNG_state, sample_radiance);
		/* fait la moyenne sur tous les rayons */
		axpy(1. / samples, sample_radiance, subpixel_radiance);
	      }
	      clamp(subpixel_radiance);
	      /* fait la moyenne sur les 4 sous-pixels */
	      axpy(0.25, subpixel_radiance, pixel_radiance);
	    }
	  }
	  if(rank==root) copy(pixel_radiance, image + 3 * ((h - 1 - cur_h) * w + j)); // <-- retournement vertical
	  else copy(pixel_radiance, ima_res->image + 3 * ((nlines - 1 - i) * w + j)); // <-- retournement vertical

	  /* Vérifie si un autre processus n'est pas en attente
	   * i.e un processus a fini son lot de travail et demande
	   * au processus actuel s'il lui reste du travail
	   */
	  MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);

	  if(flag) { /* si oui */
	    tag=status.MPI_TAG;
	    if(tag==TAG_REQ) { // Demande de travail
	      fprintf(stderr, "[%d] : attente vérifiée pour %d\n", rank, status.MPI_SOURCE);
	      MPI_Recv(&buf, 1, MPI_INT, status.MPI_SOURCE, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	      int nl=nblocs_left(blocs, nblocs, ind);
	      if(nl==0) {/* il ne reste plus de bloc pour ce process */
		fprintf(stderr, "[%d] : passage de la requête de %d à %d\n", rank, buf, (rank+1)%np);
		MPI_Send(&buf, 1, MPI_INT, (rank+1)%np, TAG_REQ, MPI_COMM_WORLD);
	      }
	      else {
		/* Envoi des blocs */
		int count=ceil(nl/2.0);
		MPI_Send(blocs+ind, count, MPI_INT, buf, TAG_DATA, MPI_COMM_WORLD);
		fprintf(stderr, "[%d] : %d blocs envoyés à %d\n", rank, count, buf);
		for(int i=0; i<count; i++) {
		  blocs[ind]=-1;
		  ind=(ind+1)%nblocs; /* incrémentation vers le prochain bloc de travail */
		}
	      }
	    } else if(tag==TAG_END) { // Fin de tous les travaux
	      // le process fait suivre le msg
	      MPI_Recv(&buf, 1, MPI_INT, status.MPI_SOURCE, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	      MPI_Isend(&buf, 1, MPI_INT, (rank+1)%np, tag, MPI_COMM_WORLD, &request);
	      fprintf(stderr, "[%d] : Tag end transmis à %d\n", rank, (rank+1)%np);
	      end_sent=true;
	      flag_end=true;
	    }
	    flag=0;
	  }
	} /* fin for j */
      } /* fin for i */

      /* ajout du resultat à la liste */
      push_ima(list_res, ima_res);

    }/* Fin des travaux (boucle while 2) */
    

    /* Entrée dans l'ANNEAU LOGIQUE */
    if(!end_sent) {
      MPI_Send(&rank, 1, MPI_INT, (rank+1)%np, TAG_REQ, MPI_COMM_WORLD);
      fprintf(stderr, "[%d] : entrée en entente\n", rank);
      while(true) {
	MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	tag=status.MPI_TAG;
	if(tag == TAG_REQ) {
	  MPI_Recv(&buf, 1, MPI_INT, status.MPI_SOURCE, tag, MPI_COMM_WORLD, &status);
	  if(buf == rank) { // Le msg envoyé a fait le tour de l'anneau -> terminaison
	    MPI_Isend(&rank, 1, MPI_INT, (rank+1)%np, TAG_END, MPI_COMM_WORLD, &request);
	    fprintf(stderr, "[%d] : tag end envoyé à %d\n", rank, (rank+1)%np);
	    end_sent=true;
	    break;
	  }
	  else { // Nouvelle demande de travail d'un autre proc -> fait suivre
	    MPI_Isend(&buf, 1, MPI_INT, (rank+1)%np, TAG_REQ, MPI_COMM_WORLD, &request);
	    fprintf(stderr, "[%d] : passage de la requête de %d à %d\n", rank, buf, (rank+1)%np);
	  }
	} else if(tag==TAG_DATA){
	  int count;
	  MPI_Get_count(&status, MPI_INT, &count);
	  /* Si on arrive ici, le tableau 'blocs' entier vaut -1. On stcoke donc au debut */
	  MPI_Recv(blocs, count, MPI_INT, status.MPI_SOURCE, TAG_DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  fprintf(stderr, "[%d] : %d blocs reçu de %d\n", rank, count, status.MPI_SOURCE);
	  ind=0;
	  end_sent=false;
	  break;
	} else {
	  MPI_Recv(&buf, 1, MPI_INT, status.MPI_SOURCE, TAG_END, MPI_COMM_WORLD, &status);
	  fprintf(stderr, "[%d] : Tag end reçu de %d\n", rank, status.MPI_SOURCE);
	  MPI_Isend(&buf, 1, MPI_INT, (rank+1)%np, TAG_END, MPI_COMM_WORLD, &request);
	  fprintf(stderr, "[%d] : Tag end transmis à %d\n", rank, (rank+1)%np);
	  end_sent=true;
	  flag_end=true;
	  break;
	}
      }
    }

    while(end_sent && !flag_end) {
      fprintf(stderr, "[%d] : En attente...\n", rank);
      MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      tag=status.MPI_TAG;
      if(tag == TAG_REQ) {
	MPI_Recv(&buf, 1, MPI_INT, status.MPI_SOURCE, tag, MPI_COMM_WORLD, &status);
      } else if(tag == TAG_END){ // TAG_END
	MPI_Recv(&buf, 1, MPI_INT, status.MPI_SOURCE, tag, MPI_COMM_WORLD, &status);
	fprintf(stderr, "[%d] : Tag end reçu de %d. Terminaison...\n", rank, status.MPI_SOURCE);
	flag_end=true;
	break;
      }
    }
  } /* fin boucle principale */

  /* Répartition finale des blocs de tâche */
  fprintf(stderr, "\n[%d] : Répartition finale des blocs de tâche : ", rank);
  print_list(list_res);

  /* Regroupement des données calculées */
  fprintf(stderr, "%d dans Gather\n", rank);
  image_res_t* ima;
  int nbb; /* Numero du bloc */
  if(rank==root) {
    int i;
    /* root reçoit les images des autres processeurs */
#pragma omp parallel for schedule(runtime) private(tag, status, nbb) ordered
    for(i=0; i<np; i++) {
      if(i != rank) {
	fprintf(stderr, "root receving from %d...\n", i);
	MPI_Recv(&nbb, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	tag=status.MPI_TAG;
	while(tag != TAG_END) {
	  if(tag==TAG_DATA) {
#pragma omp ordered
	    MPI_Recv(image+3*w*nbb*nlines, 3*w*nlines, MPI_DOUBLE, i, TAG_DATA, MPI_COMM_WORLD, &status);
	  }
	  MPI_Recv(&nbb, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	  tag=status.MPI_TAG;
	}
      }
    }
  } else {
    while(!list_empty(list_res)) {
      ima=pop_ima(list_res);
      MPI_Ssend(&ima->nbloc, 1, MPI_INT, root, TAG_DATA, MPI_COMM_WORLD);
      MPI_Ssend(ima->image, 3*w*nlines, MPI_DOUBLE, root, TAG_DATA, MPI_COMM_WORLD);
    }
    nbb=-1;
    MPI_Send(&nbb, 1, MPI_INT, root, TAG_END, MPI_COMM_WORLD);
  }

  /* libération des tampons */
  free_list(list_res);
  free(blocs);

  /* fin du chronometrage */
  fin = wtime();
  if(rank==root) fprintf( stderr, "\nTemps total de calcul : %g sec\n", fin - debut);

  /* root stocke l'image dans un fichier au format NetPbm */
  if(rank==root){
    fprintf(stderr, "\nSAUVEGARDE EN COURS...\n");
    struct passwd *pass;
    char nom_sortie[100] = "";
    char nom_rep[30] = "";

    pass = getpwuid(getuid());
    sprintf(nom_rep, "/tmp/%s", pass->pw_name);
    mkdir(nom_rep, S_IRWXU);
    sprintf(nom_sortie, "%s/par-mpi-image-omp.ppm", nom_rep);

    FILE *f = fopen(nom_sortie, "w");
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w * h; i++) {
      fprintf(f,"%d %d %d ", toInt(image[3 * i]), toInt(image[3 * i + 1]), toInt(image[3 * i + 2]));
    }
    fclose(f);
    fprintf(stderr, "\nSauvegardé.\nOpen file %s for result\n", nom_sortie);
    free(image);
  }

  MPI_Finalize();
		
}

