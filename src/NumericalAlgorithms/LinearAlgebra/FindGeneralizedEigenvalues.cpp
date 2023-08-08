// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "NumericalAlgorithms/LinearAlgebra/FindGeneralizedEigenvalues.hpp"

#include <cstddef>
#include <ostream>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Matrix.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"

// LAPACK routine to do the generalized eigenvalue problem
extern "C" {
// The final two arguments are the "hidden" lengths of the first two.
// https://gcc.gnu.org/onlinedocs/gfortran/Argument-passing-conventions.html
extern void dggev_(char*, char*, int*, double*, int*, double*, int*, double*,
                   double*, double*, double*, int*, double*, int*, double*,
                   int*, int*, size_t, size_t);
}

void find_generalized_eigenvalues(
    const gsl::not_null<DataVector*> eigenvalues_real_part,
    const gsl::not_null<DataVector*> eigenvalues_imaginary_part,
    const gsl::not_null<Matrix*> eigenvectors, Matrix matrix_a,
    Matrix matrix_b) {
  // Sanity checks on the sizes of the vectors and matrices
  const size_t number_of_rows = matrix_a.rows();
  ASSERT(number_of_rows == matrix_a.columns(),
         "Matrix A should be square, but A has "
             << matrix_a.rows() << " rows and " << matrix_a.columns()
             << " columns.");
  ASSERT(number_of_rows == matrix_b.rows() and
             number_of_rows == matrix_b.columns(),
         "Matrix A and matrix B should be the same size, but A has "
             << matrix_a.rows() << " rows and " << matrix_a.columns()
             << " columns, while B has " << matrix_b.rows() << " rows and "
             << matrix_b.columns() << " columns.");
  ASSERT(number_of_rows == eigenvectors->rows() and
             number_of_rows == eigenvectors->columns(),
         "Matrix A and matrix eigenvectors should have the same size, "
             "but A has " << matrix_a.rows() << " rows and "
             << matrix_a.columns() << " columns, while the eigenvectors matrix "
             << "has " << eigenvectors->rows() << " rows and "
             << eigenvectors->columns() << " columns.");
  ASSERT(number_of_rows == eigenvalues_real_part->size() and
             number_of_rows == eigenvalues_imaginary_part->size(),
         "eigenvalues DataVector sizes should equal number of columns "
         "in Matrix A, but A has "
             << matrix_a.columns()
             << " columns, while the real eigenvalues DataVector size is "
             << eigenvalues_real_part->size()
             << " and the imaginary eigenvalues DataVector size is "
             << eigenvalues_imaginary_part->size() << ".");

  // Set up parameters for the lapack call
  // Lapack uses chars to decide whether to compute the left eigenvectors,
  // the right eigenvectors, both, or neither. 'N' means do not compute,
  // 'V' means do compute. Note: not const because lapack does not want this
  // option const.
  char compute_left_eigenvectors = 'N';
  char compute_right_eigenvectors = 'V';

  // Lapack expects the sizes to be ints, not size_t.
  // NOTE: not const because lapack function dggev_() arguments
  // are not const.
  auto matrix_and_vector_size = static_cast<int>(number_of_rows);

  // Lapack splits the eigenvalues into unnormalized real and imaginary
  // parts, which it calls alphar and alphai, and a normalization,
  // which it calls beta. The real and imaginary parts of the eigenvalues are
  // found by dividing the unnormalized results by the normalization.
  DataVector eigenvalue_normalization(number_of_rows, 0.0);

  // Lapack uses a work vector, that should have a size 8N
  // for doing eigenvalue problems with NxN matrices
  // Note: a non-const int, not size_t, because lapack wants a non-const int
  int work_size = number_of_rows * 8;
  std::vector<double> lapack_work(static_cast<size_t>(work_size), 0.0);

  //  Lapack uses an integer called info to return its status
  //  info = 0 : success
  //  info = -i: ith argument had bad value
  //  info > 0: some other failure
  int info = 0;

  int matrix_a_spacing = matrix_a.spacing();
  int matrix_b_spacing = matrix_b.spacing();
  int eigenvectors_spacing = eigenvectors->spacing();

  dggev_(&compute_left_eigenvectors, &compute_right_eigenvectors,
         &matrix_and_vector_size, matrix_a.data(), &matrix_a_spacing,
         matrix_b.data(), &matrix_b_spacing, eigenvalues_real_part->data(),
         eigenvalues_imaginary_part->data(), eigenvalue_normalization.data(),
         eigenvectors->data(), &eigenvectors_spacing, eigenvectors->data(),
         &eigenvectors_spacing, lapack_work.data(), &work_size, &info, 1, 1);

  if (UNLIKELY(info != 0)) {
    ERROR(
        "Lapack failed to compute generalized eigenvectors. Lapack's dggev "
        "INFO = "
        << info);
  }

  *eigenvalues_real_part /= eigenvalue_normalization;
  *eigenvalues_imaginary_part /= eigenvalue_normalization;
}
