#include <stdio.h>
#include <stdlib.h>
#include "list_res.h"


image_res_t* create_image_res(int n, int w, int h) {
	image_res_t * res = (image_res_t *)malloc(sizeof(image_res_t));
	res->image=(double*)malloc(sizeof(res->image) * 3 * w * h);
	if (res->image == NULL) {
		perror("Impossible d'allouer l'image");
		exit(1);
	}
	res->nbloc=n;
	res->suivant = NULL;
	return res;
}

void free_image_res(image_res_t* ima) {
	free(ima->image);
	free(ima);
}

list_res_t* allocate_list(){
	list_res_t *res = (list_res_t*)malloc(sizeof(list_res_t));
	if (res == NULL) {
		perror("Impossible d'allouer la liste");
		exit(1);
	}
	res->first = NULL;
	res->length = 0;
	return res;
}

bool list_empty(list_res_t *l){
  return l->length==0;
}

/**
 * Ajout au debut de la liste
 */
void push_ima(list_res_t *l, image_res_t *new_ima){
	new_ima->suivant=l->first;
	l->first=new_ima;
	l->length++;
}

/**
 * Extraction au debut
 */
image_res_t* pop_ima(list_res_t *l){
	image_res_t* res = l->first;
	if(!list_empty(l)){
		l->first=res->suivant;
		res->suivant=NULL;
		l->length--;
	}
	return res;
}

void free_list(list_res_t *l){
  while(!list_empty(l))
    free_image_res(pop_ima(l));
  free(l);
}

void print_list(list_res_t* l) {
	image_res_t* ima;
	fprintf(stderr, "{");
	for(ima=l->first; ima != NULL; ima=ima->suivant) {
		fprintf(stderr, "%d ", ima->nbloc);
	}
	fprintf(stderr, "}\n");
}

