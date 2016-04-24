#ifndef GUARD_TMatrix_h
#define GUARD_TMatrix_h

//-------------------------------------------------------------------------------------------
//
// TMatrix.h
//
// Author: Struan Robertson
// Date:   30/Mar/2003
//
// This header file contains the declaration of the TMatrix class. This class inherits from
// Matrix and wraps calls to EISPACK functions.
//
//-------------------------------------------------------------------------------------------
// #include "Matrix.h"
#include <string>
#include <cstring>
#include <cmath>
#include <climits>
#include <stdio.h>
#include <vector>
#include <complex>
#include <stdexcept>
#include <stdlib.h>
#include "Persistence.h"
#include "qd/dd_real.h"
#include "qd/qd_real.h"

#include <Eigen/Dense>

namespace mesmer
{
	template<class T>
	class TMatrix : public Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> {

	public:

		//
		// Constructors
		//
		explicit TMatrix(size_t n, const T& init = T(0.0)) : m_msize(n), Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(n, n) {
			for (size_t i = 0; i < m_msize; i++)
				for (size_t j = 0; j < m_msize; j++)
					(*this)(i, j) = init;
		};

		//
		// Copy constructor
		//
		TMatrix(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& rhs) : Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(rhs) { };

		// Destructor

		virtual ~TMatrix() {};

		// Operators

		TMatrix& operator=(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& rhs) {

			if (this != &rhs) {

				// If necessary resize underlying array.

				if (m_msize != rhs.size()) {
					m_msize = rhs.size() ;
					this->resize(m_msize, m_msize);
				}
				(*this) = rhs;
			}
			return *this;

		};

		inline T* operator[](const size_t i) { return &(*this)(i); }

		inline const T* operator[](const size_t i) const { return &(*this)(i); }

		// Accessors

		size_t size() const { return m_msize; }

		// Modifiers
		void reset(const size_t n) {

			if (n < 1) {
				cerr << "Matrix must be of size one or greater";
			}

			this->resize(n, n);

			for (size_t i = 0; i < n; i++)
				for (size_t j = 0; j < n; j++)
					(*this)(i, j) = T(0.0);

		}

		//
		// Wrapped calls to diagonalise matrix.
		//
		void diagonalize(T *rr) {

			Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > eigensolver((*this));

			Eigen::Matrix<T, Eigen::Dynamic, 1> egv = eigensolver.eigenvalues();

			for (size_t i = 0; i < m_msize; i++)
				rr[i] = egv[i];

			(*this) = eigensolver.eigenvectors();

		}

		//
		// Solve a set of linear equations with a single right hand side.
		//
		void solveLinearEquationSet(T *rr) {

			Eigen::Matrix<T, Eigen::Dynamic, 1> rhs(m_msize);
			for (size_t i = 0; i < m_msize; i++)
				rhs[i] = rr[i];

			Eigen::Matrix<T, Eigen::Dynamic, 1> x = this->colPivHouseholderQr().solve(rhs);

			for (size_t i = 0; i < m_msize; i++)
				rr[i] = x[i];

		};

		// Determinant of Matrix.
		T Determinant() { return this->determinant(); };

		// Code adapted from C acording to the algorithm given at rosettacode.org/wiki/Cholesky_decomposition.
		void cholesky() {

		};

		// Matrix inversion method by LU decomposition
		bool invertLUdecomposition() {
			try {
				this->inverse();
			}
			catch (...) 
			   { return true; }
			return false;
		} ;

		void normalizeProbabilityMatrix();

		// Print out the contents of the matrix.
		void print(std::string& title, std::ostream& output_stream, int nr = -1, int nc = -1, int fr = -1, int fc = -1) const;

		// Apply the Gram-Schmidt procedure to orthogonalize the current matrix.
		void GramSchimdt(size_t root_vector);

		// Transpose matrix.
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>&  Transpose() {
			this->transpose();
			return (*this);
		};

		// Write matrix to an XML stream.
		void WriteToXML(PersistPtr pp);

		void showFinalBits(const size_t n, bool isTabbed = false);

	private:

		//
		// The following two methods permute a matrix, in order to reduce the effects
		// of numerical rounding in the Househlider method (see numerical recipes) and
		// unpermute the associated eignvector matrix. SHR 7/Mar/2011: this experiment
		// appeared to have only a minor effect on numerical values, hence the statements
		// in diagonalize() have been commented out as the benefits do not appear to out
		// way the cost in CPU time at present.
		//
		//void permuteMatrix(T **a, vector<size_t>& index);

		//void unPermuteEigenEigenvectors(T **a, vector<size_t>& index);

		size_t m_msize;

	};


	// Matrix vector mutiplication operator.
	template<class T>
	void operator*=(std::vector<T>& lhs, const TMatrix<T>& rhs) {

		size_t msize = rhs.size();
		if (lhs.size() != msize) {
			// Throw error.
		}

		Eigen::Matrix<T, Eigen::Dynamic, 1> result(msize);
		for (size_t i = 0; i < msize; i++)
			result[i] = lhs[i];

		result = rhs*result;

		for (size_t i = 0; i < msize; i++)
			lhs[i] = result[i];

	}

	//
	// This method permutes the matrix, in order to reduce the effects
	// of numerical rounding in the Householder method (see numerical recipes).
	//
	//template<class T>
	//void TMatrix<T>::permuteMatrix(T **a, vector<size_t>& index) {

	//	size_t size = this->size();
	//	vector<T> diagonalElements(size, 0.0);

	//	for (size_t i(0); i < size; i++) {
	//		diagonalElements[i] = fabs(a[i][i]);
	//		index[i] = i;
	//	}

	//	// Order matrix columns based on the magnitude of the diagonal elements.

	//	for (size_t ii(1); ii < size; ++ii) {
	//		size_t i = ii - 1;
	//		size_t k = i;
	//		T p = diagonalElements[i];
	//		for (size_t j(ii); j < size; ++j) {
	//			if (diagonalElements[j] < p) {
	//				k = j;
	//				p = diagonalElements[j];
	//			}
	//		}
	//		if (k != i) {
	//			diagonalElements[k] = diagonalElements[i];
	//			diagonalElements[i] = p;
	//			swap(index[i], index[k]);

	//			// Swap columns.
	//			for (size_t j(0); j < size; ++j) {
	//				swap(a[j][i], a[j][k]);
	//			}

	//			// Swap rows.
	//			for (size_t j(0); j < size; ++j) {
	//				swap(a[i][j], a[k][j]);
	//			}
	//		}
	//	}

	//}

	//
	// This method unpermutes the eignvector matrix rows following diagonalization.
	//
	//template<class T>
	//void TMatrix<T>::unPermuteEigenEigenvectors(T **a, vector<size_t>& index) {

	//	size_t size = this->size();

	//	vector<size_t> invIndex(size, 0);
	//	for (size_t i(0); i < size; i++) {
	//		invIndex[index[i]] = i;
	//	}

	//	for (size_t i(0); i < size; i++) {
	//		size_t k = invIndex[i];
	//		if (k != i) {
	//			invIndex[index[i]] = k;
	//			invIndex[i] = i;
	//			swap(index[i], index[k]);
	//			for (size_t j(0); j < size; j++) {
	//				swap(a[i][j], a[k][j]);
	//			}
	//		}
	//	}

	//}

	//
	// Normalization of Probability matrix.
	// Normalising coefficients are found by using the fact that column sums
	// are unity. The procedure leads to a matrix that is of upper triangular
	// form and the normalisation constants are found by back substitution.
	//
	template<class T>
	void TMatrix<T>::normalizeProbabilityMatrix() {

		int i, j; //int makes sure the comparison to negative numbers meaningful (i >=0)

		int optrsize(int(this->size()));
		vector<T> work(optrsize);// Work space.

		T scaledRemain(0.0);
		for (i = optrsize - 1; i >= 0; --i) {

			T upperSum(0.0);
			for (j = 0; j <= i; ++j)
				upperSum += (*this)[j][i];

			if (upperSum > 0.0) {
				if (i < optrsize - 1) {
					scaledRemain = 0.0;
					for (j = i + 1; j < optrsize; ++j) {
						T scale = work[j];
						scaledRemain += (*this)[j][i] * scale;
					}
				}
				work[i] = (1.0 - scaledRemain) / upperSum;
			}
		}

		//
		// Apply normalization coefficients
		//
		for (i = 0; i < optrsize; ++i) {
			(*this)[i][i] *= work[i];
			//T value = (*this)[i][i];
			for (j = i + 1; j < optrsize; ++j) {
				(*this)[j][i] *= work[j];
				(*this)[i][j] *= work[j];
			}
		}

	}

	//
	// Print out the contents of the matrix.
	//
	template<class T>
	void TMatrix<T>::print(std::string& title, std::ostream& output_stream, int nr, int nc, int fr, int fc) const {

		size_t msize = this->size();
		size_t nrows = (nr < 0) ? msize : min(msize, size_t(nr));
		size_t nclms = (nc < 0) ? msize : min(msize, size_t(nc));
		size_t frow = (fr < 0) ? 0 : min(msize, size_t(fr));
		size_t fclm = (fc < 0) ? 0 : min(msize, size_t(fc));

		output_stream << endl << title << endl << "{" << endl;
		for (size_t i(frow); i < nrows; ++i) {
			for (size_t j(fclm); j < nclms; ++j) {
				formatFloat(output_stream, (*this)(i, j), 6, 15);
				output_stream << ",";
			}
			output_stream << endl;
		}
		output_stream << "}" << endl;

	}

	//
	// Apply the Gram-Schmidt procedure to orthogonalize the current matrix.
	//
	template<class T>
	void TMatrix<T>::GramSchimdt(size_t root_vector) {

		size_t size = this->size();

		for (int i(size - 1); i > -1; i--) { // Need to use int here as size_t is unsigned.

			size_t j;
			T sum(T(0.0));
			//
			// Orthogonalize vector (Remove projections).
			//
			for (j = (size - 1); j > size_t(i); j--) {
				sum = T(0.0);
				for (size_t k = 0; k < size; k++) {
					sum += (*this)[k][j] * (*this)[k][i];
				}

				for (size_t k = 0; k < size; k++) {
					(*this)[k][i] -= sum * (*this)[k][j];
				}
			}

			//
			// Normalize vector.
			//
			sum = T(0.0);
			size_t l;
			for (l = 0; l < size; l++) {
				sum += (*this)[l][i] * (*this)[l][i];
			}
			// sum = sqrt(sum);
			for (l = 0; l < size; l++) {
				(*this)[l][i] /= sum;
			}

		}

	}

	template<class T>
	void TMatrix<T>::WriteToXML(PersistPtr pp)
	{
		//Write a CML element <matrix> as child of pp
		stringstream ss;
		for (size_t i(0); i < this->size(); i++) {
			for (size_t j(0); j <= i; j++)
				ss << to_double((*this)[i][j]) << ' ';
			ss << '\n';
		}
		PersistPtr ppmatrix = pp->XmlWriteValueElement("matrix", ss.str(), true);
		// The "true" parameter puts the matrix values in a CDATA wrapper so that
		// the new lines are preserved. If the parameter is omitted the data
		// is all space separated. Both form are read identically.
		ppmatrix->XmlWriteAttribute("matrixType", "squareSymmetricLT");
		ppmatrix->XmlWriteAttribute("rows", toString(this->size()));
	}

	template<class T>
	void TMatrix<T>::showFinalBits(const size_t n, bool isTabbed) {

		size_t size = this->size();

		// if n == 0, print the whole matrix
		size_t fb(n);
		if (n == 0) fb = size;

		//
		// Show the final n x n square of the current matrix
		//
		ctest << "{\n";
		if (!isTabbed) {
			for (size_t i = size - fb; i < size; ++i) {
				for (size_t j = size - fb; j < size; ++j) {
					formatFloat(ctest, (*this)[i][j], 5, 13);
				}
				ctest << endl;
			}
		}
		else {
			for (size_t i = size - fb; i < size; ++i) {
				for (size_t j = size - fb; j < size; ++j) {
					ctest << (*this)[i][j] << "\t";
				}
				ctest << endl;
			}
		}
		ctest << "}\n";
	}

	//-------------------------------------------------------------------------------------------
	// Function tqlev - eigenvalue only method based on tqli.
	//
	// QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors of a real,
	// symmetric, tridiagonal matrix, or of a real, symmetric matrix previously reduced by tred2.
	//
	// On input:
	//    d[1..n] contains the diagonal elements of the tridiagonal matrix.
	//    e[1..n] contains the subdiagonal elements of the tridiagonal matrix.
	// with e[1] arbitrary.
	// On output:
	//    d[1..n] returns the eigenvalues.
	//    e[1..n] is destroyed.
	//
	//-------------------------------------------------------------------------------------------
	template<class T>
	void tqlev(T *d, T *e, size_t n)
	{
		size_t m, l, iter, i, k;
		T s, r, p, g, f, dd, c, b;

		if (n == 0) return;

		for (i = 2; i <= n; ++i) e[i - 2] = e[i - 1];
		e[n - 1] = 0.0;
		for (l = 1; l <= n; ++l) {
			iter = 0;
			do {
				for (m = l; m <= n - 1; ++m) {
					dd = fabs(d[m - 1]) + fabs(d[m]);
					if (fabs(e[m - 1]) + dd == dd) break;
				}
				if (m != l) {
					// if (iter++ == 30) fprintf(stderr, "Too many iterations in tqlev");
					if (iter++ == 60) {
						fprintf(stderr, "Too many iterations in tqlev");
						throw(std::runtime_error("Too many iterations in tqlev.")); ;
					}
					/* CHL
					Source: http://www.nr.com/forum/showthread.php?t=592
					I hope that bellow words will be useful for you.
					See thread under the title: Possible convergence problems in svdcmp, jacobi, tqli, hqr by Saul Teukolsky
					in Forum: Official Bug Reports with known bugs. May be this is a reason of slow convergency.
					It is good check, that matrix is symmetric and to work with double accuracy. I have known also versions
					with increased number of iterations (200 for example). But I think that this experimental number is right
					in any case: if you have not convergency for 30 iterations, there is no convergency at all.
					SVD method used in book is an intrinsic iterative procedure, 30 iterations is a good number to
					convergency up to numerical accuracy. Evgeny
					*/
					g = (d[l] - d[l - 1]) / (2.0*e[l - 1]);
					r = sqrt((g*g) + 1.0);
					g = d[m - 1] - d[l - 1] + e[l - 1] / (g + (g < 0.0 ? -fabs(r) : fabs(r)));
					s = c = 1.0;
					p = 0.0;
					for (i = m - 1; i >= l; --i) {
						f = s*e[i - 1];
						b = c*e[i - 1];
						if (fabs(f) >= fabs(g)) {
							c = g / f;
							r = sqrt((c*c) + 1.0);
							e[i] = f*r;
							c *= (s = 1.0 / r);
						}
						else {
							s = f / g;
							r = sqrt((s*s) + 1.0);
							e[i] = g*r;
							s *= (c = 1.0 / r);
						}
						g = d[i] - p;
						r = (d[i - 1] - g)*s + 2.0*c*b;
						p = s*r;
						d[i] = g + p;
						g = c*r - b;
					}
					d[l - 1] = d[l - 1] - p;
					e[l - 1] = g;
					e[m - 1] = 0.0;
				}
			} while (m != l);
		}

		// Order eigenvalues.

		for (size_t ii = 1; ii < n; ++ii) {
			i = ii - 1;
			k = i;
			p = d[i];
			for (size_t j = ii; j < n; ++j) {
				if (d[j] < p) {
					k = j;
					p = d[j];
				}
			}
			if (k != i) {
				d[k] = d[i];
				d[i] = p;
			}
		}

	}

	// Read a square CML <matrix>:
	// The data type of the values is T and they are separated by whitespace.
	// ppmatrix points to the <matrix> element, which must have an attribute
	// rows="n" where n is a non-zero integer.
	// If it has an attribute matrixType="squareSymmetricUT" the values are the upper triangle;
	// if it is "squareSymmetricLT" they are the lower triangle; and 
	// if it is omitted or is anything else all elements are assumed present.
	// The returned matrix is fully populated.
	template<class T>
	TMatrix<T>* ReadMatrix(PersistPtr ppmatrix) {
		size_t nrows(0);
		if (ppmatrix && (nrows = ppmatrix->XmlReadInteger("rows", false)) != 0)
		{
			bool upper = false, lower = false;
			const char* ptype = ppmatrix->XmlReadValue("matrixType", false);
			if (strcmp(ptype, "squareSymmetricLT") == 0)
				lower = true;
			else if (strcmp(ptype, "squareSymmetricUT") == 0)
				upper = true;

			TMatrix<T>* m = new TMatrix<T>(nrows);
			std::stringstream ss;
			ss.str(ppmatrix->XmlRead());
			for (size_t nr = 0; nr < nrows; ++nr)
			{
				size_t a = upper ? nr : 0;
				size_t b = lower ? nr : nrows - 1;
				for (size_t nc = a; nc <= b; ++nc)
				{
					double val;
					ss >> val;
					(*m)[nr][nc] = (*m)[nc][nr] = val;
				}
			}
			return m;
		}
		return NULL; //empty
	}

	// Read a symmetrical square matrix from a CML <property> element, given
	// a pointer to the molecule or propertyList. See ReadMatrix() for details.
	template<class T>
	TMatrix<T>* ReadPropertyMatrix(const std::string& name, PersistPtr ppparent) {
		PersistPtr ppmatrix = ppparent->XmlMoveToProperty(name); //to <matrix>
		return ReadMatrix<T>(ppmatrix);
	}

}//namespace mesmer


#endif // GUARD_TMatrix_h
