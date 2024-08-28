#include "s21_matrix.h"

int not_real_matrix(matrix_t *A) {
  return (A == NULL || A->rows < 1 || A->columns < 1 || A->matrix == NULL) ? 1
                                                                           : 0;
}

int not_same_size(matrix_t *A, matrix_t *B) {
  return (A->rows != B->rows || A->columns != B->columns) ? 1 : 0;
}

int not_squere(matrix_t *A) { return (A->columns != A->rows) ? 1 : 0; }

double minor_m(int x, int y, matrix_t *A) {
  double res = 1;
  matrix_t M;
  if (!s21_create_matrix(A->rows - 1, A->columns - 1, &M)) {
    int r = 0, c;
    for (int i = 0; i < A->rows; i++) {
      if (i == x) continue;
      c = 0;
      for (int j = 0; j < A->columns; j++) {
        if (j == y) continue;
        M.matrix[r][c] = A->matrix[i][j];
        c++;
      }
      r++;
    }
    s21_determinant(&M, &res);
    s21_remove_matrix(&M);
  }
  return res;
}

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int res = 0;
  if (rows < 1 || columns < 1 || result == NULL)
    res = 1;
  else {
    result->rows = rows;
    result->columns = columns;
    result->matrix =
        calloc(1, rows * sizeof(double *) + rows * columns * (sizeof(double)));
    if (result->matrix == NULL)
      res = 1;
    else
      for (int i = 0; i < rows; i++)
        result->matrix[i] = (void *)result->matrix + rows * sizeof(double *) +
                            i * columns * sizeof(double);
  }
  return res;
}

void s21_remove_matrix(matrix_t *A) {
  if (A->matrix != NULL) free(A->matrix);
  A->matrix = NULL;
  A->rows = 0;
  A->columns = 0;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  if (not_real_matrix(A) || not_real_matrix(B) || not_same_size(A, B))
    return FAILURE;
  for (int i = 0; i < A->rows; i++)
    for (int j = 0; j < A->columns; j++)
      if (fabs(A->matrix[i][j] - B->matrix[i][j]) > ACCURACY) return FAILURE;
  return SUCCESS;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (not_real_matrix(A) || not_real_matrix(B)) return 1;
  if (not_same_size(A, B)) return 2;
  s21_create_matrix(A->rows, A->columns, result);
  for (int i = 0; i < A->rows; i++)
    for (int j = 0; j < A->columns; j++)
      result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
  return 0;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (not_real_matrix(A) || not_real_matrix(B)) return 1;
  if (not_same_size(A, B)) return 2;
  s21_create_matrix(A->rows, A->columns, result);
  for (int i = 0; i < A->rows; i++)
    for (int j = 0; j < A->columns; j++)
      result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
  return 0;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  if (not_real_matrix(A)) return 1;
  s21_create_matrix(A->rows, A->columns, result);
  for (int i = 0; i < A->rows; i++)
    for (int j = 0; j < A->columns; j++)
      result->matrix[i][j] = number * A->matrix[i][j];
  return 0;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (not_real_matrix(A) || not_real_matrix(B)) return 1;
  if (A->columns != B->rows) return 2;
  s21_create_matrix(A->rows, B->columns, result);
  for (int i = 0; i < A->rows; i++)
    for (int j = 0; j < B->columns; j++)
      for (int m = 0; m < A->columns; m++)
        result->matrix[i][j] =
            result->matrix[i][j] + A->matrix[i][m] * B->matrix[m][j];
  return 0;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  if (not_real_matrix(A)) return 1;
  s21_create_matrix(A->columns, A->rows, result);
  for (int i = 0; i < A->rows; i++)
    for (int j = 0; j < A->columns; j++) result->matrix[j][i] = A->matrix[i][j];
  return 0;
}

int s21_determinant(matrix_t *A, double *result) {
  if (not_real_matrix(A)) return 1;
  if (not_squere(A)) return 2;
  matrix_t B;
  s21_create_matrix(A->rows, A->columns, &B);
  for (int i = 0; i < A->rows; i++)
    for (int j = 0; j < A->columns; j++) B.matrix[i][j] = A->matrix[i][j];
  int sign = 1, row_num, flag;
  double temp;
  *result = 1;
  for (int k = 0; k < B.columns - 1; k++) {
    flag = 0;
    for (int i = k; i < B.rows; i++)
      if (B.matrix[i][k] != 0) {
        row_num = i;
        flag = 1;
        break;
      }
    if (flag) {
      if (row_num != k) {
        sign *= -1;
        for (int j = 0; j < B.columns; j++) {
          temp = B.matrix[k][j];
          B.matrix[k][j] = B.matrix[row_num][j];
          B.matrix[row_num][j] = temp;
        }
      }
      double first_in_row;
      for (int i = k + 1; i < B.rows; i++) {
        first_in_row = B.matrix[i][k] / B.matrix[k][k];
        for (int j = k; j < B.columns; j++) {
          B.matrix[i][j] -= B.matrix[k][j] * first_in_row;
        }
      }
    } else
      *result = 0;
  }
  for (int i = 0; i < B.rows; i++) *result *= B.matrix[i][i];
  *result *= sign;
  s21_remove_matrix(&B);
  return 0;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  if (not_real_matrix(A)) return 1;
  if (not_squere(A)) return 2;
  s21_create_matrix(A->rows, A->columns, result);
  for (int i = 0; i < A->rows; i++)
    for (int j = 0; j < A->columns; j++)
      result->matrix[i][j] = minor_m(i, j, A) * pow(-1, i + j);
  return 0;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  if (not_real_matrix(A)) return 1;
  double det;
  s21_determinant(A, &det);
  if (not_squere(A) || !det) return 2;
  matrix_t B;
  s21_calc_complements(A, &B);
  s21_transpose(&B, result);
  for (int i = 0; i < A->rows; i++)
    for (int j = 0; j < A->columns; j++) result->matrix[i][j] /= det;
  s21_remove_matrix(&B);
  return 0;
}