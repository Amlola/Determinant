#include <iostream>
#include "matrix.hpp"

int main() {

    try {
        size_t n;
        std::cin >> n;

        MatrixOperation::Matrix<double> matrix(n, n);
        matrix.FillMatrix(std::cin);

        std::cout << matrix.Determinant() << "\n";

    } catch (const std::exception& err) {
        std::cerr << "Error: " << err.what() << '\n';
        return EXIT_FAILURE;
    }
    
    return 0;
}