#define _CRT_SECURE_NO_WARNINGS
#if !defined MATRIX_H
#define MATRIX_H

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h> // Per avere uint8_t
#include <stdbool.h>
#include <math.h>
struct matrix {
	size_t rows;
	size_t cols;
	double *data;
};

struct bmatrix {
	size_t rows, cols;
	bool *data;
};

enum comparisons {
	LT, LE, EQ,
	NE, GE, GT
};

struct image {
	size_t rows, cols;
	uint8_t *data;
};
struct forest {
	size_t rows, cols;
	char *data;
};



//altre funzioni
extern struct image *image_doublesize(const struct image *img);
extern struct image *aggiungi_cornice(const struct image *img);
extern void propagate_fire(const struct forest *f);


//funzioni gestione matrice
extern void print_matrix(const struct matrix *A);
extern void prod_per_scalare(struct matrix *mat, double k);
extern struct matrix *mat_copy(struct matrix *dst, struct matrix const *src);
extern struct matrix *mat_constructor(struct matrix *m, int rows, int cols);
extern void mat_destructor(struct matrix *m);
extern struct matrix *new_mat(int rows, int cols);
extern void delete_mat(struct matrix *m);
extern struct matrix *mat_create_copy(struct matrix *src);
extern double *diag(const struct matrix *matr, size_t n);
extern double det3x3(const struct matrix *matr);
extern void matrix_write(const struct matrix *matr, FILE *f);
extern int matrix_read(struct matrix *matr, FILE *f);
extern int mat_isupper(const struct matrix *matr);
extern void mat_swaprows(struct matrix *mat, size_t r1, size_t r2);
extern void mat_swapcols(struct matrix *mat, size_t c1, size_t c2);
extern struct matrix *mat_transpose(const struct matrix *mat);
extern struct matrix *mat_replica(const struct matrix *a, int dir); //NON FUNZIONA
extern struct matrix *mat_rendiquadrata(const struct matrix *a);
extern struct matrix *mat_creatediag(const double *values, size_t n);
extern struct matrix *matrix_flip_v(const struct matrix *m);
extern struct matrix *matrix_flip_h(const struct matrix *m);
extern double *bordo_esterno(const struct matrix *m, size_t *new_size);//funziona in parte
extern struct matrix *rotate_v(const struct matrix *m, int n); //da implementare
extern struct matrix **leggi_matrici(const char *filename, size_t *size); //da implementare
extern struct matrix *mat_permute_rows(const struct matrix *m, const size_t *p);
extern struct bmatrix *mat_boolean(const struct matrix *m, double rhs, enum comparisons cmp);
extern struct matrix *matrix_quadruplica(const struct matrix *m);
#endif //MATRIX_H