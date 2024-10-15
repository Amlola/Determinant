#include <cmath>
#include <cstddef>
#include <gtest/gtest.h>
#include <string>
#include <vector>
#include <fstream>
#include "matrix.hpp"

namespace {
    static constexpr double DEVIATION  = 1e-3;

    double GetDifference(double num1, double num2) {

        return std::fabs(num1 - num2);
    }
}

TEST(MatrixTest, DetTraceTestFile) {

    const std::string files[3] = {"../data/test1.det_trace", "../data/test2.det_trace", "../data/test3.det_trace"};

    std::ifstream results("../expected_data/results.det_trace");

    size_t num_rows = 0;
    
    for (size_t cur_file = 0; cur_file < 3; cur_file++) {

        std::ifstream cur_matrix(files[cur_file]);
        cur_matrix >> num_rows;

        MatrixOperation::Matrix<double> matrix(num_rows, num_rows);
        matrix.FillMatrix(cur_matrix);

        cur_matrix.close();

        double det_result = 0;
        results >> det_result;

        double trace_result = 0;
        results >> trace_result;
        
        EXPECT_LT(GetDifference(matrix.Determinant(), det_result), DEVIATION) << "Wrong Determinant in test_file: " << cur_file + 1;
        EXPECT_LT(GetDifference(matrix.Trace(), trace_result),     DEVIATION) << "Wrong Trace in test_file: "       << cur_file + 1;
    }

    results.close();
}

TEST(MatrixTest, IdentityMatrixTest) {

    MatrixOperation::Matrix<int>* matrix = MatrixOperation::Matrix<int>::IdentityMatrix(100, 100);

    EXPECT_LT(GetDifference(matrix->Determinant(), 1), DEVIATION);
    EXPECT_LT(GetDifference(matrix->Trace(), 100),     DEVIATION);
    
    delete matrix;
}

TEST(MatrixTest, ItCtorTest) {

    std::vector<int> elems{1, 2, 3,
                           4, 5, 6,
                           7, 8, 9};

    MatrixOperation::Matrix<int> matrix{3, 3, elems.begin(), elems.end()};

    ASSERT_LT(GetDifference(matrix.Determinant(), 0), DEVIATION);
}

TEST(MatrixTest, MoveAssighmentTest) {

    MatrixOperation::Matrix<double> matrix1(3, 5, 2);

    MatrixOperation::Matrix<double> matrix2(3, 5);
    matrix2 = matrix1;

    ASSERT_EQ(matrix1.Equal(matrix2), 1);

    MatrixOperation::Matrix<double> matrix3(3, 5);
    matrix3 = std::move(matrix1);

    ASSERT_EQ(matrix3.Equal(matrix2), 1);
}

TEST(MatrixTest, CopyCtorTest) {

    MatrixOperation::Matrix<double> matrix1(3, 5, 2);

    MatrixOperation::Matrix<double> matrix2 = matrix1;

    ASSERT_EQ(matrix1.Equal(matrix2), 1);

    MatrixOperation::Matrix<double> matrix3 = std::move(matrix1);

    ASSERT_EQ(matrix3.Equal(matrix2), 1);
}

TEST(MatrixTest, NegTestFile) {

    const std::string files[2] = {"../data/test1.neg_transp", "../data/test2.neg_transp"};

    const std::string expected_results_neg[2] = {"../expected_data/results1.neg"   , "../expected_data/results2.neg"};

    size_t num_rows = 0;
    size_t num_cols = 0;
    
    for (size_t cur_file = 0; cur_file < 2; cur_file++) {

        std::ifstream cur_matrix(files[cur_file]);
        cur_matrix >> num_rows;
        cur_matrix >> num_cols;

        MatrixOperation::Matrix<double> matrix(num_rows, num_cols);
        matrix.FillMatrix(cur_matrix);

        matrix.Negate();

        cur_matrix.close();

        std::ifstream cur_neg_result(expected_results_neg[cur_file]);
        MatrixOperation::Matrix<double> result_matrix(num_rows, num_cols);

        result_matrix.FillMatrix(cur_neg_result);

        cur_neg_result.close();
        
        ASSERT_EQ(matrix.Equal(result_matrix), 1) << "Wrong Negate in test_file: " << cur_file + 1;
    }
}

TEST(MatrixTest, TranspTestFile) {

    const std::string files[2] = {"../data/test1.neg_transp", "../data/test2.neg_transp"};

    const std::string expected_results_transp[2] = {"../expected_data/results1.transp", "../expected_data/results2.transp"};

    size_t num_rows = 0;
    size_t num_cols = 0;
    
    for (size_t cur_file = 0; cur_file < 2; cur_file++) {

        std::ifstream cur_matrix(files[cur_file]);
        cur_matrix >> num_rows;
        cur_matrix >> num_cols;

        MatrixOperation::Matrix<double> matrix(num_rows, num_cols);
        matrix.FillMatrix(cur_matrix);

        matrix.Transpose();

        cur_matrix.close();

        std::ifstream cur_transp_result(expected_results_transp[cur_file]);
        MatrixOperation::Matrix<double> result_matrix(num_cols, num_rows);

        result_matrix.FillMatrix(cur_transp_result);

        cur_transp_result.close();
        
        ASSERT_EQ(matrix.Equal(result_matrix), 1) << "Wrong Transposed in test_file: " << cur_file + 1;
    }
}

TEST(MatrixTest, DitryTest) {

    MatrixOperation::Matrix<int> matrix1(3, 3, 0);

    ASSERT_LT(GetDifference(matrix1.Determinant(), 0), DEVIATION);
    ASSERT_LT(GetDifference(matrix1.Trace(),       0), DEVIATION);

    MatrixOperation::Matrix<int> transpose_matrix = matrix1;
    transpose_matrix.Transpose();

    ASSERT_EQ(matrix1.Equal(transpose_matrix), 1);

    MatrixOperation::Matrix<int> negate_matrix = matrix1;
    negate_matrix.Negate();

    ASSERT_EQ(matrix1.Equal(negate_matrix), 1);
}

int main(int argc, char **argv) {

    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
