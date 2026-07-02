#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "Core/Block.impl.hpp"

template <typename SBM>
void print(const SBM& sbm) {
  for (int row = 0; row < sbm.rowsOfBlocks(); ++row) {
    for (typename SBM::InnerIterator it(sbm.innerIterator(row)); it; ++it) {
      std::cout << "Block at (" << row << ", " << it.inner() << ") is \n" << sbm.block(it.ptr()) << "\n";
    }
  }
}

// template<typename SBM, typename Rhs>
// Rhs solve_LU(const SBM& sbm, const Rhs& rhs) {
//     Rhs result = rhs;
//     // TODO: Implement a block LU solver here that recursively solves small
//     // Schur complement blocks or even better a class that stores the local
//     // block decompositions such that the linear system can be solved multiple
//     // times with the same decomposition.
//     for( int row = 0 ; row < sbm.rowsOfBlocks() ; ++ row ) {
//         for( typename SBM::InnerIterator it ( sbm.innerIterator( row ) ) ; it ; ++it ) {
//             // Do something here...
//         }
//     }
//     return rhs;
// }

// template<typename Block>
// class BlockLU
// {
// public:

//     using SparseBlockMatrix = bogus::SparseBlockMatrix<Block>;
//     using DenseLU = Eigen::PartialPivLU<Block>;

//     void compute(const SparseBlockMatrix& A) {
//         int nBlocks = A.size();
//         int k;
//         int pivot = 0;
//         decltype(A.innerIterator(0)) pivot_ptr;
//         for( int row = 0 ; row < A.rowsOfBlocks() ; ++ row ) {
//             // find pivot block
//             pivot = row;
//             typename SparseBlockMatrix::InnerIterator pivot_it(A.innerIterator(row));
//             typename SparseBlockMatrix::InnerIterator it(A.innerIterator(row));
//             for(; it; ++it) {
//                 const Block& block = A.block( it.ptr() );
//                 DenseLU lu(block);
//                 if (lu.info() == Eigen::Success) {
//                     // m_pivot_rows.push_back(row);
//                     // m_pivot_cols.push_back(it.col());
//                     m_pivot_idx.push_back(it.ptr());
//                     m_diagLU.push_back(lu);
//                     break;
//                 }
//                 ++pivot;
//             }

//             if(pivot == nBlocks)
//                 throw std::runtime_error("Singular block matrix");

//             // if(pivot != k)
//             // if (pivot_it != it)
//             std::cout << "pivot_it.ptr(): " << pivot_it.ptr() << std::endl;
//             std::cout << "it.ptr(): " << it.ptr() << std::endl;
//             if (pivot_it.ptr() != it.ptr())
//             {
//                 // swap
//                 std::cout << "swap\n";
//                 // std::swap(A[k], A[pivot]);
//                 // std::swap(m_perm[k], m_perm[pivot]);
//             }
//         }

//         // for(k=0; k<nBlocks; ++k) {
//         //         // Find pivot block
//         //         int pivot = k;

//         //         while(pivot < nBlocks)
//         //         {
//         //             Block Akk = A[pivot][k];

//         //             Eigen::PartialPivLU<Block> lu(Akk);

//         //             if(lu.isInvertible()) {
//         //                 break;
//         //             }

//         //             ++pivot;
//         //         }

//         //         if(pivot == nBlocks)
//         //             throw std::runtime_error("Singular block matrix");

//         //         if(pivot != k)
//         //         {
//         //             std::swap(A[k], A[pivot]);
//         //             std::swap(m_perm[k], m_perm[pivot]);
//         //         }

//         //         m_diagLU[k].compute(A[k][k]);

//         //         for(i=k+1; i<nBlocks; ++i)
//         //         {
//         //             if(A[i][k].size()==0)
//         //                 continue;

//         //             L[i][k] =
//         //                 m_diagLU[k]
//         //                     .solve(A[i][k].transpose())
//         //                     .transpose();

//         //             for(j=k+1; j<nBlocks; ++j)
//         //             {
//         //                 if(A[k][j].size()==0)
//         //                     continue;

//         //                 A[i][j] -= L[i][k] * A[k][j];
//         //             }
//         //         }
//         //     }
//     }

//     template<typename Rhs>
//     Eigen::VectorXd solve(const Rhs& rhs) const;

// // private:

//     std::vector<int> m_pivot_rows;
//     std::vector<int> m_pivot_cols;
//     std::vector<int> m_pivot_idx;
//     std::vector<DenseLU> m_diagLU;

//     // std::vector<unsigned> m_offsets;

//     // Eigen::MatrixXd m_denseLU;
//     // Eigen::PartialPivLU<Eigen::MatrixXd> m_schurLU;

//     // std::vector<Eigen::PartialPivLU<Block>> m_diagLU;
//     std::vector<int> m_perm;
//     // std::vector<std::vector<Block>> m_L;
//     // std::vector<std::vector<Block>> m_U;
// };

template <typename Block>
class BlockLU {
 public:
  using Scalar = typename Block::Scalar;
  using DenseLU = Eigen::PartialPivLU<Block>;

  template <typename SBM>
  void compute(const SBM& sbm) {
    const int n = sbm.rowsOfBlocks();

    m_n = n;

    m_rowSizes.resize(n);
    m_colSizes.resize(n);

    for (int i = 0; i < n; ++i) {
      m_rowSizes[i] = sbm.blockRows(i);
      m_colSizes[i] = sbm.blockCols(i);
    }

    // m_L.assign(n, std::vector<Block>(n));
    // m_U.assign(n, std::vector<Block>(n));

    // std::vector<std::vector<Block>> A(
    //     n,
    //     std::vector<Block>(n)
    // );

    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        // A[i][j] =
        //     Block::Zero(
        //         m_rowSizes[i],
        //         m_colSizes[j]
        //     );
        m_A.push_back(Block::Zero(m_rowSizes[i], m_colSizes[j]));
      }
    }

    for (int row = 0; row < n; ++row) {
      for (typename SBM::InnerIterator it(sbm.innerIterator(row)); it; ++it) {
        // A[row][it.inner()] =
        //     sbm.block(it.ptr());
        block(row, it.inner()) = sbm.block(it.ptr());
      }
    }

    m_diagLU.resize(n);

    for (int k = 0; k < n; ++k) {
      // m_diagLU[k].compute(A[k][k]);
      m_diagLU[k].compute(block(k, k));

      if (m_diagLU[k].matrixLU().diagonal().array().abs().minCoeff() < 1e-14) {
        throw std::runtime_error("Singular pivot block");
      }

      // TODO
      // m_U[k][k] = A[k][k];

      // TODO
      // for(int j=k+1;j<n;++j)
      //     m_U[k][j] = A[k][j];

      for (int i = k + 1; i < n; ++i) {
        // if(A[i][k].size()==0)
        if (block(i, k).size() == 0) continue;

        Block Lik = m_diagLU[k].solve(block(i, k).transpose()).transpose();

        // TODO:
        // m_L[i][k] = Lik;

        for (int j = k + 1; j < n; ++j) {
          // A[i][j] -= Lik * A[k][j];
          block(i, j) -= Lik * block(k, j);
        }
      }
    }

    for (int i = 0; i < n; ++i) {
      // TODO:
      // m_L[i][i] =
      //     Block::Identity(
      //         m_rowSizes[i],
      //         m_rowSizes[i]
      //     );
    }
  }

  Eigen::VectorXd solve(const Eigen::VectorXd& rhs) const {
    std::vector<Eigen::VectorXd> b(m_n);

    {
      int off = 0;

      for (int i = 0; i < m_n; ++i) {
        b[i] = rhs.segment(off, m_rowSizes[i]);

        off += m_rowSizes[i];
      }
    }

    std::vector<Eigen::VectorXd> y(m_n);

    for (int i = 0; i < m_n; ++i) {
      y[i] = b[i];

      for (int j = 0; j < i; ++j) {
        // TODO
        // if(m_L[i][j].size())
        //     y[i] -= m_L[i][j] * y[j];
      }
    }

    std::vector<Eigen::VectorXd> x(m_n);

    for (int i = m_n - 1; i >= 0; --i) {
      Eigen::VectorXd rhs_i = y[i];

      for (int j = i + 1; j < m_n; ++j) {
        // TODO:
        // if(m_U[i][j].size())
        //     rhs_i -= m_U[i][j] * x[j];
      }

      x[i] = m_diagLU[i].solve(rhs_i);
    }

    Eigen::VectorXd result(totalRows());

    int off = 0;

    for (int i = 0; i < m_n; ++i) {
      result.segment(off, m_rowSizes[i]) = x[i];

      off += m_rowSizes[i];
    }

    return result;
  }

 private:
  int totalRows() const {
    int n = 0;

    for (unsigned s : m_rowSizes) n += s;

    return n;
  }

  // int m_n = 0;

  // std::vector<unsigned> m_rowSizes;
  // std::vector<unsigned> m_colSizes;

  // std::vector<DenseLU> m_diagLU;

  // std::vector<
  //     std::vector<Block>
  // > m_L;

  // std::vector<
  //     std::vector<Block>
  // > m_U;

  int m_n;

  std::vector<unsigned> m_rowSizes;
  std::vector<unsigned> m_colSizes;

  std::vector<Block> m_A;

  std::vector<Eigen::PartialPivLU<Block> > m_diagLU;

  std::vector<int> m_rowPerm;
  std::vector<int> m_colPerm;

  Block& block(int i, int j) { return m_A[i * m_n + j]; }

  const Block& block(int i, int j) const { return m_A[i * m_n + j]; }
};

int main() {
  // Block structure (of some mechanical system with compliant force laws)
  //
  // [ A11(3x3) A12(3x3) A13(3x2) ]
  // [ A21(3x3) A22(3x3) A23(3x2) ]
  // [ A31(2x3) A32(2x3) A33(2x2) ]
  // =
  // [   M1(3x3)    0(3x3) W1(3x2) ]
  // [    0(3x3)   M2(3x3) W2(3x2) ]
  // [ W1.T(2x3) W2.T(2x3)  C(2x2) ]
  using Block = Eigen::MatrixXd;
  bogus::SparseBlockMatrix<Block> A;

  unsigned rowSizes[3] = {3, 3, 2};
  unsigned colSizes[3] = {3, 3, 2};

  A.setRows(3, rowSizes);
  A.setCols(3, colSizes);
  A.reserve(7);

  double m1 = 1.25;
  Eigen::Matrix3d M1 = m1 * Eigen::Matrix3d::Identity();

  double m2 = 2.5;
  Eigen::Matrix3d M2 = m2 * Eigen::Matrix3d::Identity();

  Eigen::Matrix<double, 3, 2> W1;
  W1 << 1, 0, 0, 1, 0, 0;
  Eigen::Matrix<double, 3, 2> W2 = -W1;

  double k = 1e6;
  double c = 1.0 / k;
  Eigen::Matrix2d C = c * Eigen::Matrix2d::Identity();

  A.insertBack(0, 0) = M1;
  A.insertBack(0, 2) = W1;
  A.insertBack(1, 1) = M2;
  A.insertBack(1, 2) = W2;
  A.insertBack(2, 0) = W1.transpose();
  A.insertBack(2, 1) = W2.transpose();
  A.insertBack(2, 2) = C;

  A.finalize();

  // std::cout << "A:\n";
  // print(A);

  Eigen::VectorXd b(8);
  b.setRandom();

  // // TODO: Implement a block LU/ldlt solver here!
  // Eigen::VectorXd result = solve_LU(A, b);
  // std::cout << "result: " << result.transpose() << std::endl;

  // double error = (A * result - b).norm();
  // std::cout << "error: " << error << std::endl;

  BlockLU<Block> solver;
  solver.compute(A);

  Eigen::VectorXd result = solver.solve(b);

  std::cout << "error = " << (A * result - b).norm() << std::endl;

  // std::cout << "pivot blocks:\n";
  // for (auto idx : solver.m_pivot_idx) {
  //     std::cout << "- pivot block:\n" << A.block(idx) << std::endl;
  // }

  // // std::cout << "pivot_rows:\n";
  // // for (auto row : solver.m_pivot_rows) {
  // //     std::cout << "row: " << row << std::endl;
  // // }

  // // std::cout << "pivot_cols:\n";
  // // for (auto col : solver.m_pivot_cols) {
  // //     std::cout << "col: " << col << std::endl;
  // // }

  {
    // simple block matrix for iterating over it
    // TODO: Can we do this for sparse matrices as well?
    using Block = Eigen::MatrixXd;
    using BlockMatrix = Eigen::Matrix<Block, Eigen::Dynamic, Eigen::Dynamic>;

    // [   M1(3x3)    0(3x3) W1(3x2) ]
    // [    0(3x3)   M2(3x3) W2(3x2) ]
    // [ W1.T(2x3) W2.T(2x3)  C(2x2) ]
    BlockMatrix AA(3, 3);
    AA(0, 0) = M1;
    AA(0, 2) = W1;
    AA(1, 1) = M2;
    AA(1, 2) = W2;
    AA(2, 0) = W1.transpose();
    AA(2, 1) = W2.transpose();
    AA(2, 2) = C;

    for (int row = 0; row < AA.rows(); ++row) {
      for (int col = 0; col < AA.cols(); ++col) {
        std::cout << "Block at (" << row << ", " << col << ") is \n" << AA(row, col) << "\n";
      }
    }
  }
}