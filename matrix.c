#include "matrix.h"

void print_matrix(const struct matrix *a){
	print_matrix("Matrice %dx%d:\n", a->rows, a->cols);
	print_matrix("---------------\n");
	for(size_t r=0; r < a->rows; r++){
		for(size_t c=0; c < a->cols; c++){
			printf("%d", a->data[r*a->cols+c]);
		}
		printf("\n");
	}
	print_matrix("---------------\n");
}

void prod_per_scalare(struct matrix *mat, double k){
    for (int r = 0; r < mat->rows; ++r) {
        for (int c = 0; c < mat->cols; ++c) {
            mat->data[r * mat->cols + c] *= k;
        }
    }
}

struct matrix *mat_copy(struct matrix *dst, struct matrix const *src){
    int size = src->rows * src->cols;
    for (int i = 0; i < size; ++i) {
        dst->data[i] = src->data[i];
    }
    return dst;
}

struct matrix *mat_constructor(struct matrix *m, int rows, int cols){
    m->rows = rows;
    m->cols = cols;
    m->data = malloc(rows*cols * sizeof(double));
    return m;
}

void mat_destructor(struct matrix *m){
    free(m->data);
}

struct matrix *new_mat(int rows, int cols){
    return mat_constructor(malloc(sizeof(struct matrix)), 
        rows, cols);
}

void delete_mat(struct matrix *m){
    mat_destructor(m);
    free(m);
}

struct matrix *mat_create_copy(struct matrix *src){
    return mat_copy(new_mat(src->rows, src->cols), src);
}

double *diag (const struct matrix *matr){
	double *result = NULL;
	if(matr != NULL){
		result = malloc(matr->rows*sizeof(double));
		size_t i=0;
		for(size_t r = 0; r < matr->rows; r++){
			for(size_t c = 0; c < matr->cols; c++){
				if(r==c)	result[i] = matr->data[r*matr->cols+c];
			}
		} 
	}
	return result;
}

double det3x3 (const struct matrix *matr){
	double result = NAN;
	if(matr != NULL){
		result = matr->data[0]*matr->data[4]*matr->data[8] -
				 matr->data[2]*matr->data[4]*matr->data[6];
	}
	return result;
}

void matrix_write(const struct matrix *matr, FILE *f) {
	if (matr != NULL) {
		if (matr->data != NULL) {
			if (f != NULL) {
				fprintf(f, "%u\n%u\n", matr->rows, matr->cols);
				for (size_t r = 0; r < matr->rows; r++) {
					for (size_t c = 0; c < matr->cols; c++) {
						fprintf(f, "%f", matr->data[r * matr->cols + c]);
						if(c + 1 != matr->cols)	fputc('\t', f);
					}
					//a capo
					fputc('\n', f);
				}
			}
		}
	}
}

int matrix_read(struct matrix *matr, FILE *f) {
	int result = 0;
	if (f != NULL) {
		int letti = fscanf(f, "%u%u", &matr->rows, &matr->cols);
		if (letti == 2) {
			matr->data = malloc(matr->rows *matr->cols * sizeof(double));
			size_t i = 0;
			while (1) {
				double value = 0;
				int letti = fscanf(f, "%lf", &value);
				if (letti != 1) {
					if (feof(f)) {
						result = 1;
						break;
					}
					if (ferror(f)) {
						free(matr->data);
						matr->data = NULL;
						break;
					}
				}
				else {
					matr->data[i] = value;
					i++;
				}
			}
		}
	}
	return result;
}

int mat_isupper(const struct matrix *matr) {
	int result = 0;
	if (matr->rows == matr->cols) {
		result = 1;
		for (size_t r = 1; r < matr->rows; r++) {
			for (size_t c = 0; c < matr->cols; c++) {
				if (r > c) {
					if (matr->data[r * matr->cols + c] != 0) {
						result = 0;
						break;
					}
				}
			}
			if (result == 0)	break;
		}
	}
	return result;
}

void mat_swaprows(struct matrix *mat, size_t r1, size_t r2) {
	if (mat != NULL) {
		if (mat->data != NULL) {
			if (r1 < mat->rows && r2 < mat->rows) {
				double *tmp_row = malloc(mat->cols * sizeof(double));
				for (size_t c = 0, i = 0; c < mat->cols; c++, i++) {
					tmp_row[i] = mat->data[r1*mat->cols + c];	//copia tmp del vettore riga
					mat->data[r1*mat->cols + c] = mat->data[r2*mat->cols + c];
					mat->data[r2*mat->cols + c] = tmp_row[i];
				}
				free(tmp_row);
			}
		}
	}
}

void mat_swapcols(struct matrix *mat, size_t c1, size_t c2) {
	if (mat != NULL) {
		if (mat->data != NULL) {
			if (c1 < mat->cols && c2 < mat->cols) {
				double *tmp_col = malloc(mat->rows * sizeof(double));
				for (size_t r = 0, i = 0; r < mat->rows; r++, i++) {
					tmp_col[i] = mat->data[r * mat->cols + c1];	//copia tmp del vettore riga
					mat->data[r * mat->cols + c1] = mat->data[r * mat->cols + c2];
					mat->data[r * mat->cols + c2] = tmp_col[i];
				}
				free(tmp_col);
			}
		}
	}
}

struct matrix *mat_transpose(const struct matrix *mat) {
	struct matrix *result = NULL;
	if (mat != NULL) {
		if (mat->data != NULL) {
			result = malloc(sizeof(struct matrix));
			result->rows = mat->rows;
			result->cols = mat->cols;
			result->data = malloc(mat->rows*mat->cols * sizeof(double));

			for (size_t c = 0; c < result->cols; c++) {
				for (size_t r = 0; r < result->rows; r++) {
					result->data[c * result->cols + r] = mat->data[r*result->rows + c];
				}
			}
		}
	}
	return result;
}

struct matrix *mat_replica(const struct matrix *a, int dir){//non funziona
	//dir = 0 -> replica orizzontalmente
	//dir != 0 -> replica verticalmente
	
	struct matrix *result = NULL;
	
	if(a != NULL){
		if(a->data != NULL){
			result = malloc(sizeof(struct matrix));
			if(dir == 0){
				result->rows = a->rows;
				result->cols = a->cols * 2;
			} else {
				result->rows = a->rows * 2;
				result->cols = a->cols;
			}
			
			result->data = malloc(result->rows*result->cols*sizeof(double));
			
			//horizzontal
			for(size_t i = 0; i < 2; i++){
				for(size_t r = 0; r < result->rows; r++){
					for(size_t c = i*a->cols; c < result->cols; c++){
						result->data[r * a->cols + c] = a->data[r * a->cols + c];
					}
				}
			}
			
			//vertical
			
			for(size_t i = 0; i < 2; i++){
				for(size_t r = 0; r < result->rows; r++){
					for(size_t c = i*a->cols; c < result->cols; c++){
						result->data[r * a->cols + c] = a->data[r * a->cols + c];
					}
				}
			}
		}
	}
	
	return result;
}

struct matrix *mat_rendiquadrata(const struct matrix *a) {
	struct matrix *result = NULL;
	if (a != NULL) {
		if (a->data != NULL) {
			result = malloc(sizeof(struct matrix));
			if (a->rows < a->cols) {
				result->rows = a->cols;
				result->cols = a->cols;
			}
			else {
				result->rows = a->rows;
				result->cols = a->rows;
			}

			result->data = calloc(result->rows*result->cols, sizeof(double));

			if (a->rows < a->cols) {
				for (size_t r = 0; r < a->rows; r++) {
					for (size_t c = 0; c < a->cols; c++) {
						result->data[r * a->cols + c] = a->data[r * a->cols + c];
					}
				}
			}
			else {
				for (size_t c = 0; c < a->cols; c++) {
					for (size_t r = 0; r < a->rows; r++) {
						result->data[r * result->cols + c] = a->data[r * a->cols + c];
					}
				}
			}
		}
	}
	return result;
}

struct matrix *mat_creatediag(const double *values, size_t n) {
	struct matrix *result = NULL;
	result = malloc(sizeof(struct matrix));
	result->rows = n;
	result->cols = n;
	result->data = calloc(n * n, sizeof(double));
	if (n != 0) {
		for (size_t r = 0, i = 0; r < result->rows; r++, i++) {
			for (size_t c = 0; c < result->cols; c++) {
				if (r == c) {
					result->data[r*result->cols + c] = values[i];
				}
			}
		}
	}
	return result;
}

struct image *image_doublesize(const struct image *img){//non funziona
	struct image *result = NULL;
	
	return NULL;
}

struct matrix *matrix_flip_v(const struct matrix *m) {
	struct matrix *result = NULL;
	if (m != NULL) {
		result = malloc(sizeof(struct matrix));
		result->rows = m->rows;
		result->cols = m->cols;
		result->data = malloc(m->cols*m->rows * sizeof(double));
		mat_copy(result, m);
		//flip vertically
		for (size_t r = 0; r < result->rows / 2; ++r) {
			size_t  k = result->rows - 1 - r;
			for (size_t c = 0; c < result->cols; ++c) {
				double temp = result->data[r * result->cols + c];
				result->data[r * result->cols + c] = result->data[k * result->cols + c];
				result->data[k * result->cols + c] = temp;
			}
		}
	}
	return result;
}

struct matrix *matrix_flip_h(const struct matrix *m) {
	struct matrix *result = NULL;
	if (m != NULL) {
		result = malloc(sizeof(struct matrix));
		result->rows = m->rows;
		result->cols = m->cols;
		result->data = malloc(m->cols*m->rows * sizeof(double));
		mat_copy(result, m);
		//flip h
		for (size_t c = 0; c < result->cols / 2; c++) {
			size_t  k = result->cols - 1 - c;
			for (size_t r = 0; r < result->rows; r++) {
				double temp = result->data[r * result->cols + c];
				result->data[r * result->cols + c] = result->data[r*result->cols + k];
				result->data[r * result->cols + k] = temp;
			}
		}
	}
	return result;
}

double *bordo_esterno(const struct matrix *m, size_t *new_size) {
	double *result = NULL;

	if (m != NULL) {
		result = malloc(m->rows*m->cols*sizeof(double));

		size_t i = 0;
		for (size_t r = 0; r < m->rows; r++) {
			for (size_t c = 0; c < m->cols; c++) {
				if (r == 0 || c == 0 || r == m->rows - 1 || c == m->cols - 1) {
					result[i] = m->data[r*m->cols + c];
					i++;
				}
			}
		}
		result = realloc(result, i*sizeof(double));
		*new_size = i;
	}
	return result;
}

struct matrix *rotate_v(const struct matrix *m, int n){
	return NULL;
}

struct matrix *mat_permute_rows(const struct matrix *m, const size_t *p) {
	struct matrix *result = malloc(sizeof(struct matrix));

	result->rows = m->rows;
	result->cols = m->cols;
	result->data = malloc(m->cols*m->rows * sizeof(double));

	size_t r = 0;
	for (size_t i = 0; i < m->rows; i++) {
		size_t new_row_index = p[i];
		for (size_t c = 0; c < m->cols; c++) {
			result->data[r*m->cols + c] = m->data[new_row_index*m->cols + c];
		}
		r++;
	}


	return result;

}

bool cmp_eval(const double mat_el, double rhs, enum comparisons cmp) {
	switch (cmp){
		case LT:
			return mat_el < rhs;
		break; 
		case EQ:
			return mat_el == rhs;
		break;
		case LE:
			return mat_el <= rhs;
		break; 
		case NE:
			return mat_el != rhs;
		break; 
		case GE:
			return mat_el >= rhs;
		break; 
		case GT:
			return mat_el > rhs;
			break;
		default: break;
	}

	return false;
}

struct bmatrix *mat_boolean(const struct matrix *m, double rhs, enum comparisons cmp) {
	struct bmatrix *result = NULL;
	if (m != NULL) {
		result = malloc(sizeof(struct matrix));

		result->rows = m->rows;
		result->cols = m->cols;

		result->data = malloc(m->cols*m->rows * sizeof(bool));

		for (size_t r = 0; r < m->rows; r++) {
			for (size_t c = 0; c < m->cols; c++) {
				result->data[r*m->cols + c] = cmp_eval(m->data[r*m->cols + c], rhs, cmp);
			}
		}

	}
	return result;
}

struct matrix *matrix_quadruplica(const struct matrix *m) {
	struct matrix *result = NULL;
	if (m != NULL) {
		result = malloc(sizeof(struct matrix));
		result->rows = m->rows * 2;
		result->cols = m->cols * 2;

		result->data = malloc(result->rows*result->cols * sizeof(double));

		size_t r1 = 0;
		for (size_t r = 0; r < result->rows; r++) {
			size_t c1 = 0;
			if (r == m->rows)	r1 = 0;
			for (size_t c = 0; c < result->cols; c++) {
				if (c == m->cols)	c1 = 0;
				result->data[r*result->cols + c] = m->data[r1*m->cols + c1];
				c1++;
			}
			r1++;
		}
	}
	return result;
}

struct image *aggiungi_cornice(const struct image *img) {
	struct image *result = NULL;

	if (img != NULL) {
		result = malloc(sizeof(struct image));

		result->rows = img->rows + 2;
		result->cols = img->cols + 2;

		result->data = calloc(result->rows*result->cols, sizeof(uint8_t));

		for (size_t r = 0; r < img->rows; r++) {
			for (size_t c = 0; c < img->cols;c++) {
				result->data[(r + 1)*result->cols + (c + 1)] = img->data[r*img->cols + c];
			}
		}

	}

	return result;
}











