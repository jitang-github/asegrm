/* as-eGRM - Ancestry-specific Genetic Relationship Matrix
 * Chiang Lab (https://chianglab.usc.edu/) - University of Southern California
 * Copyright (C) 2024 Ji Tang, Charleston W.K. Chiang
 *
 * This program is licensed for academic research use only
 * unless otherwise stated. Contact jitang@usc.edu, charleston.chiang@med.usc.edu for
 * commercial licensing options.
 * 
 * For use in any publications please cite: Ji Tang, Charleston W.K. Chiang (2025). A genealogy-based approach for revealing ancestry-specific structures in admixed populations. The American Journal of Human Genetics, Volume 112, Issue 8, 1906 - 1922.
 *
 * Acknowledgement: A part of code is adapted from egrm (https://github.com/Ephraim-usc/egrm, Copyright (C) 2022 Caoqi Fan, Nicholas Mancuso, Charleston W.K. Chiang)
*/
#define PY_SSIZE_T_CLEAN
#include <string.h>
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#include "numpy/arrayobject.h"

typedef double DTYPE;
typedef unsigned long ITYPE;

typedef struct matrix {
  ITYPE n;
  DTYPE *data;
} matrix;

matrix* new_matrix(ITYPE n)
{
  matrix* mat;
  mat = (matrix *)malloc(sizeof(matrix));
  mat->n = n;
  mat->data = calloc(n*n,sizeof(DTYPE));
  return mat;
}

void destroy_matrix(matrix* mat)
{
  free(mat->data);
  free(mat);
}

void print_matrix(matrix* mat)
{
  ITYPE n = mat->n;
  printf("matrix of dimension %lu\n", n);
  
  ITYPE i, j;
  for (i=0; i<n; i++)
  {
    for (j=0; j<n; j++)
      printf("%lf ", mat->data[i*n + j]);
    printf("\n");
  }
}

void print_vector(DTYPE* vec, ITYPE* n)
{
  printf("vector of dimension %lu\n", n);

  ITYPE i;
  for (i=0; i<n; i++)
    printf("%lf ", vec[i]);
  printf("\n");
}

void add_square(matrix* mat, ITYPE* idx, ITYPE len_idx, DTYPE q)
{
  ITYPE* x = idx;
  ITYPE* y = idx;
  ITYPE* end = idx + len_idx;
  DTYPE* data = mat->data;
  ITYPE n = mat->n;
  
  for (; x < end; x++)
    for (y=idx; y <= x; y++)
    {
      data[*x * n + *y] += q;
    }
}

void tri_to_sq(matrix* mat, ITYPE* idx, ITYPE len_idx){
  ITYPE* x;
  ITYPE* y;
  ITYPE* end = idx + len_idx;
  DTYPE* data = mat->data;
  ITYPE n = mat->n;
  for (x=idx+1; x < end; x++)
    for (y=idx; y < x; y++)
    {
      data[*y * n + *x] = data[*x * n + *y];
    }
}

// subtract mean along axis for the elements indexed by idx
void sub_mean_by_idx(matrix* mat, ITYPE* idx, ITYPE len_idx, ITYPE axis)
{
  ITYPE* i;
  ITYPE* x;
  ITYPE* y;
  ITYPE* end = idx + len_idx;
  DTYPE* data = mat->data;
  ITYPE n = mat->n;
  DTYPE mean;
  DTYPE sum;

  if (axis==0){
    for (i=idx; i < end; i++)
    {
      sum = 0.0;
      for (x=idx; x<end; x++)
      {
        sum += data[*x * n + *i];
      };
      mean = sum/len_idx;
      for (x=idx; x<end; x++){
        data[*x * n + *i] -= mean;
      };
    }
  };

  if (axis==1){
    for (i=idx; i < end; i++)
    {
      sum = 0.0;
      for (y=idx; y < end; y++)
      {
        sum += data[*i * n + *y];
      };
      mean = sum/len_idx;
      for (y=idx; y<end; y++)
      {
        data[*i * n + *y] -= mean;
      };
    };
  };
}

// use q to divide the elements indexed by idx
void divide_by_idx(matrix* mat, ITYPE* idx, ITYPE len_idx, DTYPE q)
{
  ITYPE* x = idx;
  ITYPE* y = idx;
  ITYPE* end = idx + len_idx;
  DTYPE* data = mat->data;
  ITYPE n = mat->n;

  for (; x < end; x++)
    for (y=idx; y <= x; y++)
    {
      data[*x * n + *y] /= q;
    }
}

void set_square(matrix* mat, ITYPE* idx, ITYPE len_idx, DTYPE q)
{
  ITYPE* x = idx;
  ITYPE* y = idx;
  ITYPE* end = idx + len_idx;
  DTYPE* data = mat->data;
  ITYPE n = mat->n;
  
  for (; x < end; x++)
    for (y=idx; y <= x; y++)
    {
      data[*x * n + *y] = q;
    }
}

void add(matrix* mat_1, matrix* mat_2)
{
  ITYPE n = mat_1->n;
  DTYPE* p = mat_1->data;
  DTYPE* q = mat_2->data;
  ITYPE i;
  ITYPE j;

  for (i = 0; i < n*n; i++)
  {
    *(p+i) += *(q+i);
  }

  /* It is wrong because the input mats are non-symmetric in some cases
  for (i=0; i<n; i++){
    for (j=0; j<=i; j++){
      p[i*n + j] += q[i*n + j];
     };
  };
  */
}

void set_zeros(matrix* mat_1)
{
  ITYPE n = mat_1->n;
  DTYPE* data = mat_1->data;
  
  memset(data, 0, n*n * sizeof(DTYPE));
}

void set_ones(matrix* mat_1)
{
  ITYPE n = mat_1->n;
  DTYPE* data = mat_1->data;
  ITYPE i;
  ITYPE j;

  //memset(data, 1.0, n*n * sizeof(DTYPE));
  for (i=0; i<n; i++)
  {
    for (j=0; j<=i; j++)
      data[i*n+j] = 1;
  }
}

void set_values(matrix* mat_1, DTYPE q)
{
  ITYPE n = mat_1->n;
  DTYPE* data = mat_1->data;
  ITYPE i;
  ITYPE j;

  //memset(data, 1.0, n*n * sizeof(DTYPE));
  for (i=0; i<n; i++)
  {
    for (j=0; j<=i; j++)
      data[i*n+j] = q;
  }
}

void set_zeros_by_idx(matrix* mat, ITYPE* idx, ITYPE len_idx)
{
  ITYPE* i = idx;
  ITYPE* end = idx + len_idx;
  DTYPE* data = mat->data;
  ITYPE n = mat->n;
  ITYPE r; // for traversing all rows of the triangle
  ITYPE c; // for traversing all columns of the triangle

  for (; i < end; i++)
  {
    c = *i;
    for (r=*i; r < n; r++)
    {
      data[r*n + c] = 0;
    };
    r = *i*n;
    for (c=0; c <= *i; c++)
    {
      data[r + c] = 0;
    }
  }
}

static void parse_py_int_seq(PyObject *py_int_seq, ITYPE** pr, ITYPE* len)
{
  *len = (ITYPE)PySequence_Fast_GET_SIZE(py_int_seq);
  *pr = (ITYPE *) malloc(sizeof(ITYPE) * *len);
  ITYPE i;
  for (i = 0; i < *len; i++) 
  {
    PyObject *item = PySequence_Fast_GET_ITEM(py_int_seq, i);
    (*pr)[i] = (ITYPE)PyLong_AsLong(item);
  }
}

static PyObject* py_new_matrix(PyObject* self, PyObject* args)
{
  ITYPE n;
  PyArg_ParseTuple(args, "k", &n);

  matrix* mat = new_matrix((ITYPE)n);
  PyObject* py_mat = PyCapsule_New((void *)mat, "matrix._matrix_C_API", NULL);

  return py_mat;
}

static PyObject* py_destroy_matrix(PyObject* self, PyObject* args)
{
  PyObject* py_mat;
  PyArg_UnpackTuple(args, NULL, 1, 1, &py_mat);
  
  matrix* mat = (matrix *)PyCapsule_GetPointer(py_mat, "matrix._matrix_C_API");
  destroy_matrix(mat);
  
  Py_RETURN_NONE;
}

static PyObject* py_add_square(PyObject* self, PyObject* args)
{
  PyObject* py_mat;
  PyObject* py_idx;
  PyObject* py_q;
  PyArg_UnpackTuple(args, NULL, 3, 3, &py_mat, &py_idx, &py_q);
  
  matrix* mat = (matrix *)PyCapsule_GetPointer(py_mat, "matrix._matrix_C_API");
  DTYPE q = (DTYPE)PyFloat_AS_DOUBLE(py_q);
  
  PyObject* py_int_seq;
  py_int_seq = PySequence_Fast(py_idx, NULL);
  
  ITYPE* idx; ITYPE len_idx;
  parse_py_int_seq(py_int_seq, &idx, &len_idx);
  add_square(mat, idx, len_idx, q);
  free(idx);
  
  Py_DECREF(py_int_seq);
  Py_RETURN_NONE;
}

static PyObject* py_tri_to_sq(PyObject* self, PyObject* args)
{
  PyObject* py_mat;
  PyObject* py_idx;
  PyArg_UnpackTuple(args, NULL, 2, 2, &py_mat, &py_idx);

  matrix* mat = (matrix *)PyCapsule_GetPointer(py_mat, "matrix._matrix_C_API");

  PyObject* py_int_seq;
  py_int_seq = PySequence_Fast(py_idx, NULL);

  ITYPE* idx; ITYPE len_idx;
  parse_py_int_seq(py_int_seq, &idx, &len_idx);
  tri_to_sq(mat, idx, len_idx);
  free(idx);

  Py_DECREF(py_int_seq);
  Py_RETURN_NONE;
}

static PyObject* py_sub_mean_by_idx(PyObject* self, PyObject* args)
{
  PyObject* py_mat;
  PyObject* py_idx;
  PyObject* py_axis; // 0 or 1
  PyArg_UnpackTuple(args, NULL, 3, 3, &py_mat, &py_idx, &py_axis);

  matrix* mat = (matrix *)PyCapsule_GetPointer(py_mat, "matrix._matrix_C_API");
  ITYPE axis = (ITYPE)PyLong_AsLong(py_axis);

  PyObject* py_int_seq;
  py_int_seq = PySequence_Fast(py_idx, NULL);

  ITYPE* idx; ITYPE len_idx;
  parse_py_int_seq(py_int_seq, &idx, &len_idx);
  sub_mean_by_idx(mat, idx, len_idx, axis);
  free(idx);

  Py_DECREF(py_int_seq);
  Py_RETURN_NONE;
}

static PyObject* py_divide_by_idx(PyObject* self, PyObject* args)
{
  PyObject* py_mat;
  PyObject* py_idx;
  PyObject* py_q;
  PyArg_UnpackTuple(args, NULL, 3, 3, &py_mat, &py_idx, &py_q);

  matrix* mat = (matrix *)PyCapsule_GetPointer(py_mat, "matrix._matrix_C_API");
  DTYPE q = (DTYPE)PyFloat_AS_DOUBLE(py_q);

  PyObject* py_int_seq;
  py_int_seq = PySequence_Fast(py_idx, NULL);

  ITYPE* idx; ITYPE len_idx;
  parse_py_int_seq(py_int_seq, &idx, &len_idx);
  divide_by_idx(mat, idx, len_idx, q);
  free(idx);

  Py_DECREF(py_int_seq);
  Py_RETURN_NONE;
}

static PyObject* py_set_square(PyObject* self, PyObject* args)
{
  PyObject* py_mat;
  PyObject* py_idx;
  PyObject* py_q;
  PyArg_UnpackTuple(args, NULL, 3, 3, &py_mat, &py_idx, &py_q);
  
  matrix* mat = (matrix *)PyCapsule_GetPointer(py_mat, "matrix._matrix_C_API");
  DTYPE q = (DTYPE)PyFloat_AS_DOUBLE(py_q);
  
  PyObject* py_int_seq;
  py_int_seq = PySequence_Fast(py_idx, NULL);
  
  ITYPE* idx; ITYPE len_idx;
  parse_py_int_seq(py_int_seq, &idx, &len_idx);
  set_square(mat, idx, len_idx, q);
  free(idx);
  
  Py_DECREF(py_int_seq);
  Py_RETURN_NONE;
}

static PyObject* py_add(PyObject* self, PyObject* args)
{
  PyObject* py_mat_1;
  PyObject* py_mat_2;
  PyArg_UnpackTuple(args, NULL, 2, 2, &py_mat_1, &py_mat_2);
  
  matrix* mat_1 = (matrix *)PyCapsule_GetPointer(py_mat_1, "matrix._matrix_C_API");
  matrix* mat_2 = (matrix *)PyCapsule_GetPointer(py_mat_2, "matrix._matrix_C_API");
  
  add(mat_1, mat_2);
  
  Py_RETURN_NONE;
}

static PyObject* py_set_zeros(PyObject* self, PyObject* args)
{
  PyObject* py_mat;
  PyArg_UnpackTuple(args, NULL, 1, 1, &py_mat);
  
  matrix* mat = (matrix *)PyCapsule_GetPointer(py_mat, "matrix._matrix_C_API");
  set_zeros(mat);
  
  Py_RETURN_NONE;
}

static PyObject* py_set_ones(PyObject* self, PyObject* args)
{
  PyObject* py_mat;
  PyArg_UnpackTuple(args, NULL, 1, 1, &py_mat);

  matrix* mat = (matrix *)PyCapsule_GetPointer(py_mat, "matrix._matrix_C_API");
  set_ones(mat);

  Py_RETURN_NONE;
}
static PyObject* py_set_values(PyObject* self, PyObject* args)
{
  PyObject* py_mat;
  PyObject* py_q;
  PyArg_UnpackTuple(args, NULL, 2, 2, &py_mat, &py_q);

  matrix* mat = (matrix *)PyCapsule_GetPointer(py_mat, "matrix._matrix_C_API");
  DTYPE q = (DTYPE)PyFloat_AS_DOUBLE(py_q);

  set_values(mat, q);

  Py_RETURN_NONE;
}
static PyObject* py_set_zeros_by_idx(PyObject* self, PyObject* args)
{
  PyObject* py_mat;
  PyObject* py_idx;
  PyArg_UnpackTuple(args, NULL, 2, 2, &py_mat, &py_idx);

  matrix* mat = (matrix *)PyCapsule_GetPointer(py_mat, "matrix._matrix_C_API");

  PyObject* py_int_seq;
  py_int_seq = PySequence_Fast(py_idx, NULL);

  ITYPE* idx; ITYPE len_idx;
  parse_py_int_seq(py_int_seq, &idx, &len_idx);
  set_zeros_by_idx(mat, idx, len_idx);
  free(idx);

  Py_DECREF(py_int_seq);
  Py_RETURN_NONE;
}

static PyObject* py_print_matrix(PyObject* self, PyObject* args)
{
  PyObject* py_mat; 
  PyArg_UnpackTuple(args, NULL, 1, 1, &py_mat);
  matrix* mat = (matrix *)PyCapsule_GetPointer(py_mat, "matrix._matrix_C_API");
  print_matrix(mat);
  
  Py_RETURN_NONE;
}

static PyObject* py_export_list(PyObject* self, PyObject* args) // export to python list, elements are copied. To release memory, you have to destroy the mat_C.
{
  PyObject* py_mat; 
  PyArg_UnpackTuple(args, NULL, 1, 1, &py_mat);
  matrix* mat = (matrix *)PyCapsule_GetPointer(py_mat, "matrix._matrix_C_API");
  DTYPE* data = mat->data; ITYPE n = mat->n;
  
  PyObject *py_list = Py_BuildValue("[]");
  ITYPE i = 0;
  for (; i < n*n; i++)
  {
    PyList_Append(py_list, Py_BuildValue("d", data[i]));
  }
  
  return py_list;
}

static PyObject* py_export_ndarray(PyObject* self, PyObject* args) // export to numpy ndarray, elements are used in place. Do not destroy the mat_C.
{
  PyObject* py_mat; 
  PyArg_UnpackTuple(args, NULL, 1, 1, &py_mat);
  matrix* mat = (matrix *)PyCapsule_GetPointer(py_mat, "matrix._matrix_C_API");
  DTYPE* data = mat->data; ITYPE n = mat->n;
  
  import_array();
  ITYPE dims[2];
  dims[0] = dims[1] = n;
  PyObject *py_ndarray = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, data);
  
  return py_ndarray;
}

static PyMethodDef myMethods[] = 
{
  {"new_matrix", py_new_matrix, METH_VARARGS, "new matrix"},
  {"print_matrix", py_print_matrix, METH_VARARGS, "print matrix"},
  {"add_square", py_add_square, METH_VARARGS, "add_square"},
  {"tri_to_sq", py_tri_to_sq, METH_VARARGS, "triangle to square"},
  {"sub_mean_by_idx", py_sub_mean_by_idx, METH_VARARGS, "sub_mean_by_idx"},
  {"divide_by_idx", py_divide_by_idx, METH_VARARGS, "divide_by_idx"},
  {"set_square", py_set_square, METH_VARARGS, "set_square"},
  {"add", py_add, METH_VARARGS, "add"},
  {"set_zeros", py_set_zeros, METH_VARARGS, "set_zeros"},
  {"set_ones", py_set_ones, METH_VARARGS, "set all the elements to 1"},
  {"set_values", py_set_values, METH_VARARGS, "set all the elements to the same value"},
  {"set_zeros_by_idx", py_set_zeros_by_idx, METH_VARARGS, "set_zeros_by_idx"},
  {"destroy_matrix", py_destroy_matrix, METH_VARARGS, "destroy matrix"},
  {"export_list", py_export_list, METH_VARARGS, "export as list"},
  {"export_ndarray", py_export_ndarray, METH_VARARGS, "export as ndarray"},
  {NULL, NULL, 0, NULL},
};

static struct PyModuleDef matrixModule =
{
  PyModuleDef_HEAD_INIT,
  "matrixModule",
  "matrix Module",
  -1,
  myMethods
};

//PyMODINIT_FUNC PyInit_matrix_as(void)
PyMODINIT_FUNC PyInit_asegrm_matrix(void)
{
  return PyModule_Create(&matrixModule);
}

/*
int main()
{
  //matrix* mat = newMatrix(10);
  //int idx[3] = {2,5,7}; 
  //addSquare(mat, idx, 3, 9);
  printf("%d", fib_(10));
  return 0;
}
*/
