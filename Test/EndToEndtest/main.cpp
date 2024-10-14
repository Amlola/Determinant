#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include "matrix.hpp"

int main() {

    try {
        std::ofstream output_file("../out.txt");

        std::ifstream in_file1("../data1.txt");
        if (!in_file1) {
            return EXIT_FAILURE;
        }

        MatrixOperation::Matrix<double> matrix1(3, 5);      // Check BaseCtor and FillMatrix
        matrix1.FillMatrix(in_file1);

        in_file1.close();

        matrix1.Negate();                                   // Check Negate and DumpMatrix
        matrix1.DumpMatrix(output_file);

        MatrixOperation::Matrix<double> matrix2(3, 5);      // Check CopyAssighment and Transpose
        matrix2 = matrix1;
        matrix2.Transpose();

        MatrixOperation::Matrix<double> matrix3(3, 5);      // Check MoveAssighment
        matrix3 = std::move(matrix2);

        MatrixOperation::Matrix<int>* identity_matrix = MatrixOperation::Matrix<int>::IdentityMatrix(3, 3); // Check IdentityMatrix Ctor and Trace
        identity_matrix->Trace();  
        delete identity_matrix;

        MatrixOperation::Matrix<int> matrix_value(2, 2, 15);                // Check ValueCtor and MoveCtor
        MatrixOperation::Matrix<int> move_matrix = std::move(matrix_value); 

        std::ifstream in_file2("../data2.txt");
        if (!in_file2) {
            return EXIT_FAILURE;
        }

        std::vector<int> nums;
        int num;
        while (in_file2 >> num) {
            nums.push_back(num);
        }

        in_file2.close();

        MatrixOperation::Matrix<int> determinant_matrix(3, 3, nums.begin(), nums.end());    // Check It Ctor and Determinant
        determinant_matrix.Determinant();

        output_file.close();

    } catch (const std::exception& err) {
        std::cerr << "Error: " << err.what() << '\n';
        return EXIT_FAILURE;
    }

    return 0;
}