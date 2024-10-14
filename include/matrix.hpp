#pragma once

#include <cstddef>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <type_traits>
#include <cmath>
#include <iomanip>

namespace MatrixOperation {

    template <typename T>
    concept Arithmetic = std::is_arithmetic_v<T>;

    template<Arithmetic T>
    class Matrix_t {

    private:
        static constexpr size_t max_size_of_matrix = 1e10;

    protected:
        size_t rows = 0;
        size_t cols = 0;
        T* matrix;

        Matrix_t(size_t init_rows, size_t init_cols) {
            
            if (init_rows > max_size_of_matrix) {
                throw std::invalid_argument("Invalid number of rows");
            }

            if (init_cols > max_size_of_matrix) {
                throw std::invalid_argument("Invalid number of columns");
            }

            rows = init_rows;
            cols = init_cols;

            matrix = new T[rows * cols]();
        }

        Matrix_t(const Matrix_t<T>& other) : rows(other.rows), cols(other.cols), matrix(new T[other.cols * other.rows]) {

            std::copy(other.matrix, other.matrix + rows * cols, matrix);
        }

        Matrix_t(Matrix_t<T>&& other) noexcept : rows(other.rows), cols(other.cols), matrix(other.matrix) {

            other.matrix = nullptr;
            other.rows = 0;
            other.cols = 0;
        }

        Matrix_t& operator=(const Matrix_t<T>& other) {

            if (this != &other) {
                delete[] matrix;

                rows = other.rows;
                cols = other.cols;
                matrix = new T[rows * cols];
                std::copy(other.matrix, other.matrix + rows * cols, matrix);
            }

            return *this;
        }

        Matrix_t& operator=(Matrix_t<T>&& other) noexcept {

            if (this != &other) {
                std::swap(rows, other.rows);
                std::swap(cols, other.cols);
                std::swap(matrix, other.matrix);
            }

            return *this;
        }

        ~Matrix_t() {

            delete[] matrix;
        }
    };

    template<Arithmetic T>
    class Matrix final : private Matrix_t<T> {

        using Matrix_t<T>::rows;
        using Matrix_t<T>::cols;
        using Matrix_t<T>::matrix;

         struct ProxyRow final {

            T* row;

            T& operator[](size_t n) {
                return row[n];
            }

            T operator[](size_t n) const {
                return row[n];
            }
        };

        bool CompareDouble(double first_number, double second_number) const {

            static constexpr double EPS = 1e-9;

            return std::fabs(first_number - second_number) < EPS;
        }

        const T& GetElement(size_t cur_row, size_t cur_col) const {

            return matrix[cur_row * cols + cur_col];
        }

        void SetElement(size_t cur_row, size_t cur_col, T value) {

            matrix[cur_row * cols + cur_col] = value;
        }

        void SetElement(size_t ind, T value) {

            matrix[ind] = value;
        }
        
    public:
        Matrix(size_t init_rows, size_t init_cols) : Matrix_t<T>(init_rows, init_cols) {};

        Matrix(size_t init_rows, size_t init_cols, const T& value) : Matrix_t<T>(init_rows, init_cols) {

            const size_t& size_of_matrix = init_cols * init_rows;

            for (size_t cur_elem_ind = 0; cur_elem_ind < size_of_matrix; cur_elem_ind++) {
                SetElement(cur_elem_ind, value);
            }
        }

        template<typename It>
        Matrix(size_t init_rows, size_t init_cols, It start, It fin) : Matrix_t<T>(init_rows, init_cols) {

            size_t cur_ind = 0;
            const size_t& size_of_matrix = init_rows * init_cols;
            It& it = start;

            while (cur_ind < size_of_matrix && it < fin) {
                SetElement(cur_ind, static_cast<T>(*it));
                it++;
                cur_ind++;
            }
        }

        static Matrix<T>* IdentityMatrix(size_t init_rows, size_t init_cols) {

            if (init_cols != init_rows) {
                throw std::invalid_argument("The parameters must correspond to a square matrix");
            }

            Matrix<T>* identity_matrix = new Matrix<T>(init_rows, init_cols, 0);

            for (size_t cur_row = 0; cur_row < init_rows; cur_row++) {
                identity_matrix->SetElement(cur_row, cur_row, 1);
            }

            return identity_matrix;
        }

        template <Arithmetic T2>
        Matrix(const Matrix<T2>& other) : Matrix_t<T>(other.GetNumberRows(), other.GetNumberCols()) {

            for(size_t cur_row = 0; cur_row < rows; cur_row++) {
                for (size_t cur_col = 0; cur_col < cols; cur_col++) {
                    SetElement(cur_row, cur_col, other[cur_row][cur_col]);
                }
            }
        }

        size_t GetNumberRows() const {

            return rows;
        }

        size_t GetNumberCols() const {

            return cols;
        }

        size_t GetSizeOfMatrix() const {

            return rows * cols;
        }

        ProxyRow operator[](size_t n) const {

            return ProxyRow{.row = matrix + n * cols};
        }

        Matrix<T>& Transpose() & {

            Matrix<T> transpose_matrix(cols, rows);

            for (size_t cur_row = 0; cur_row < rows; cur_row++) {
                for (size_t cur_col = 0; cur_col < cols; cur_col++) {
                    transpose_matrix.SetElement(cur_col, cur_row, GetElement(cur_row, cur_col));
                }
            }

            *this = std::move(transpose_matrix);

            return *this;
        }

        Matrix<T>& Negate() & {

            for (size_t cur_row = 0; cur_row < rows; cur_row++) {
                for (size_t cur_col = 0; cur_col < cols; cur_col++) {
                    SetElement(cur_row, cur_col, -GetElement(cur_row, cur_col));
                }
            }

            return *this;
        }

        T Trace() const {

            if (rows != cols) {
                throw std::invalid_argument("The parameters must correspond to a square matrix");
            }

            T trace = 0;
            for (size_t cur_row = 0; cur_row < rows; cur_row++) {
                trace += GetElement(cur_row, cur_row);
            }

            return trace;
        }

        bool Equal(const Matrix& other) const {

            if (other.GetNumberRows() != rows) {
                return false;
            }

            if (other.GetNumberCols() != cols) {
                return false;
            }

            for (size_t cur_row = 0; cur_row < rows; cur_row++) {
                for (size_t cur_col = 0; cur_col < cols; cur_col++) {
                    if (!CompareDouble((*this)[cur_row][cur_col], other[cur_row][cur_col])) {
                        return false;
                    }
                }
            }

            return true;
        }

    private:
        double GetDiagonalMult(const Matrix<double>& tmp_matrix) const {

            double determinant = 1.0;

            for (size_t cur_row = 0; cur_row < rows; cur_row++) {
                determinant *= tmp_matrix[cur_row][cur_row];
            }

            return determinant;
        }

        void MakeZeroBelowDiagonal(const Matrix<double>& tmp_matrix, size_t cur_row) const {

            for (size_t next_row = cur_row + 1; next_row < rows; next_row++) {
                const double& coeff = tmp_matrix[next_row][cur_row] / tmp_matrix[cur_row][cur_row];

                for (size_t elem_ind = cur_row; elem_ind < cols; elem_ind++) {
                    tmp_matrix[next_row][elem_ind] -= coeff * tmp_matrix[cur_row][elem_ind];
                }
            }
        }

    public:
        double Determinant() {

            if (rows != cols) {
                throw std::invalid_argument("The parameters must correspond to a square matrix");
            }

            unsigned int num_swaps = 0;   

            Matrix<double> tmp_matrix(*this);

            for (size_t cur_row = 0; cur_row < rows; cur_row++) {
                if (CompareDouble(tmp_matrix[cur_row][cur_row], 0)) {
                    bool is_swapped = false;

                    for (size_t next_row = cur_row + 1; next_row < rows; next_row++) {
                        if (!CompareDouble(tmp_matrix[next_row][cur_row], 0)) {
                            for (size_t elem_ind = 0; elem_ind < cols; elem_ind++) {
                                std::swap(tmp_matrix[cur_row][elem_ind], tmp_matrix[next_row][elem_ind]);
                            }

                            num_swaps++;
                            is_swapped = true;
                            break;
                        }
                    }

                    if (!is_swapped) {
                        return 0;
                    }
                }

                MakeZeroBelowDiagonal(tmp_matrix, cur_row);
            }

            double determinant = GetDiagonalMult(tmp_matrix);

            if (num_swaps % 2 != 0) {
                return -determinant;
            }

            return determinant;
        }

        void DumpMatrix(std::ostream& os) const {

            for (size_t cur_row = 0; cur_row < rows; cur_row++) {
                for (size_t cur_col = 0; cur_col < cols; cur_col++) {
                    os << GetElement(cur_row, cur_col) << " ";
                }

                os << "\n";
            }

            os << "\n";
        }

        void FillMatrix(std::istream& is) {

            for (size_t cur_row = 0; cur_row < rows; cur_row++) {
                for (size_t cur_col = 0; cur_col < cols; cur_col++) {
                    is >> (*this)[cur_row][cur_col];
                }
            }
        }
    };
}   