/**
 * La liste des résultats obtenus sur chaque portion travaillé par un process :
 * cette liste est simplement chaînées.
*/

#ifndef _LIST_RES_H
#define _LIST_RES_H

#include <stdbool.h>

/**
   Image resultat
*/
typedef struct image_res_s{
  int nbloc;
  double *image;
  struct image_res_s *suivant;
} image_res_t;


/**
 Liste chainée des résultats
*/
typedef struct list_res_s{
  image_res_t *first;
  int length;
} list_res_t;


image_res_t* create_image_res(int n, int w, int h);
void free_image_res(image_res_t* ima);
list_res_t* allocate_list();
bool list_empty(list_res_t* l);
void push_ima(list_res_t* l, image_res_t* ima);
image_res_t* pop_ima(list_res_t* l);
void free_list(list_res_t* l);
void print_list(list_res_t* l);


#endif
