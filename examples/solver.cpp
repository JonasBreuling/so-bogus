#include <Eigen/Dense>

#include "Core/Block.impl.hpp"
#include "Core/BlockSolvers/GaussSeidel.impl.hpp"
#include "Core/BlockSolvers/Krylov.impl.hpp"
#include "Core/BlockSolvers/GaussSeidel.impl.hpp"
#include "Core/BlockSolvers/ProjectedGradient.impl.hpp"
#include "Core/BlockSolvers/ProductGaussSeidel.impl.hpp"

#include "Core/BlockSolvers/NoneLaw.hpp"
// #include "Core/BlockSolvers/LCPLaw.impl.hpp"
// #include "Core/BlockSolvers/PyramidLaw.impl.hpp"

// TODO: The issue is that the block size of the off diagonal matrices does not 
// match the NoneLaw with dimension one. This should either be dynamically sized 
// or adapted.
int main()
{
    typedef Eigen::MatrixXd Block;

    bogus::SparseBlockMatrix<Block> A;

    // Block structure:
    //
    // [ A11(3x3)  A12(3x2) ]
    // [ A21(2x3)  A22(2x2) ]
    //

    unsigned rowSizes[2] = {3, 2};
    unsigned colSizes[2] = {3, 2};

    A.setRows(2, rowSizes);
    A.setCols(2, colSizes);

    A.reserve(4);

    Eigen::Matrix3d A11;
    A11 << 4, 1, 0,
           1, 4, 1,
           0, 1, 4;

    Eigen::Matrix<double, 3, 2> A12;
    A12 << 1, 0,
           0, 1,
           0, 0;

    Eigen::Matrix<double, 2, 3> A21 = A12.transpose();

    Eigen::Matrix2d A22;
    A22 << 3, 0,
           0, 3;

    A.insertBack(0,0) = A11;
    A.insertBack(0,1) = A12;
    A.insertBack(1,0) = A21;
    A.insertBack(1,1) = A22;

    A.finalize();

    // Global vector dimension = 3 + 2 = 5
    Eigen::VectorXd b(5);
    b << 1, 2, 3, 4, 5;

    Eigen::VectorXd x = Eigen::VectorXd::Zero(5);

    {
        std::cout << "Krylov solver:\n";

        bogus::Krylov<
            bogus::SparseBlockMatrix<Block>
        > solver(A);

        solver.setMaxIters(1000);
        solver.setTol(1e-12);

        solver.solve(
            b,
            x,
            bogus::krylov::CG
        );

        std::cout << "solution: " << x.transpose() << std::endl;

        // test error
        Eigen::VectorXd error = A * x - b;
        std::cout << "error: " << error.transpose() << std::endl;
    }

    {
        std::cout << "Projected gradient solver:\n";

        bogus::ProjectedGradient<
            bogus::SparseBlockMatrix<Block>
        > solver(A);

        solver.setMaxIters(1000);
        solver.setTol(1e-12);

        solver.solve(
            bogus::NoneLaw< double >(),
            b,
            x
        );

        std::cout << "solution: " << x.transpose() << std::endl;

        // test error
        Eigen::VectorXd error = A * x - b;
        std::cout << "error: " << error.transpose() << std::endl;
    }

    {
        std::cout << "Gauss-Seidel solver:\n";

        bogus::GaussSeidel<
            bogus::SparseBlockMatrix<Block>
        > solver(A);

        solver.setMaxIters(1000);
        solver.setTol(1e-12);

        solver.solve(
            bogus::NoneLaw< double >(),
            b,
            x
        );

        std::cout << "solution: " << x.transpose() << std::endl;

        // test error
        Eigen::VectorXd error = A * x - b;
        std::cout << "error: " << error.transpose() << std::endl;
    }
}