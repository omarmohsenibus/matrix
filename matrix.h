#if !defined MATRIX_H
#define MATRIX_H

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>

struct matrix {
    size_t rows;
    size_t cols;
	void *data;
};

extern void matrix_write(const struct matrix *matr, FILE *f);
extern int matrix_read(struct matrix *matr, FILE *f);
extern struct matrix *mat_rendiquadrata(const struct matrix *a);
extern struct matrix *mat_permute_rows(const struct matrix *m, const size_t *p);
extern struct matrix *mat_scale(const struct matrix *m, double x);
extern struct matrix *scambia_diagonali(const struct matrix *m);
extern struct matrix **leggi_matrici(const char *filename, size_t *size);
extern struct matrix *prod_kronecker(const struct matrix *a, const struct matrix *b);
extern double *bordo_esterno(const struct matrix *m, size_t *new_size);
extern struct matrix *matrix_flip_h(const struct matrix *m);
extern struct matrix *matrix_flip_v(const struct matrix *m);
extern struct matrix *mat_creatediag(const double *values, size_t n);
extern struct matrix *mat_rendiquadrata(const struct matrix *a);
extern struct matrix *mat_replica(const struct matrix *a, int dir);
extern struct matrix *mat_sommadiretta(const struct matrix *a, const struct matrix *b);
extern struct matrix *mat_copy(const struct matrix *mat);
extern struct matrix *mat_transpose(const struct matrix *mat);
extern void mat_swaprows(struct matrix *mat, size_t r1, size_t r2);
extern int mat_isupper(const struct matrix *matr);


#endif //MATRIX_H

