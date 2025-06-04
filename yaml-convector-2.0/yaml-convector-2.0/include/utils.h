#pragma once
#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <cmath>
#include <cstdarg>

namespace YamlConvector2 {

inline std::string formatString(const std::string& fmt) {
    return fmt;
}

template <typename T, typename... Args>
std::string formatString(const std::string& fmt, T value, Args... args) {
    std::ostringstream out;
    std::size_t pos = fmt.find("{}");
    if (pos == std::string::npos) {
        return fmt; // nothing to format
    }
    out << fmt.substr(0, pos);
    out << value;
    out << formatString(fmt.substr(pos + 2), args...);
    return out.str();
}

class CanteraError : public std::runtime_error {
public:
    CanteraError(const std::string& origin, const std::string& msg)
        : std::runtime_error(origin + ": " + msg), m_msg(origin + ": " + msg) {}
    const std::string& getMessage() const { return m_msg; }
private:
    std::string m_msg;
};

template <typename... Args>
inline void warn_user(const std::string& origin, const std::string& fmt, Args... args) {
    std::cerr << "Warning [" << origin << "]: " << formatString(fmt, args...) << std::endl;
}

inline void writelog(const std::string& msg) {
    std::cout << msg;
}

template <typename... Args>
inline void writelogf(const char* fmt, Args... args) {
    std::printf(fmt, args...);
}

class DenseMatrix {
public:
    DenseMatrix() : m_rows(0), m_cols(0) {}
    DenseMatrix(std::size_t r, std::size_t c, double val=0.0) { resize(r,c,val); }
    void resize(std::size_t r, std::size_t c, double val=0.0) {
        m_rows = r; m_cols = c; m_data.assign(r*c, val);
    }
    double& operator()(std::size_t i, std::size_t j) { return m_data[i*m_cols + j]; }
    double operator()(std::size_t i, std::size_t j) const { return m_data[i*m_cols + j]; }
    std::size_t nRows() const { return m_rows; }
    std::size_t nCols() const { return m_cols; }
private:
    std::size_t m_rows, m_cols;
    std::vector<double> m_data;
};

inline int solve(DenseMatrix& A, double* b) {
    std::size_t n = A.nRows();
    for (std::size_t i=0;i<n;i++) {
        // pivot
        std::size_t piv = i;
        double maxVal = std::fabs(A(i,i));
        for (std::size_t k=i+1;k<n;k++) {
            if (std::fabs(A(k,i)) > maxVal) { piv = k; maxVal = std::fabs(A(k,i)); }
        }
        if (maxVal < 1e-300) return -1;
        if (piv != i) {
            for (std::size_t j=0;j<n;j++) std::swap(A(i,j), A(piv,j));
            std::swap(b[i], b[piv]);
        }
        double diag = A(i,i);
        for (std::size_t j=i;j<n;j++) A(i,j) /= diag;
        b[i] /= diag;
        for (std::size_t k=i+1;k<n;k++) {
            double factor = A(k,i);
            for (std::size_t j=i;j<n;j++) A(k,j) -= factor*A(i,j);
            b[k] -= factor*b[i];
        }
    }
    for (int i=n-1; i>=0; --i) {
        for (int k=0; k<i; ++k) {
            double factor = A(k,i);
            A(k,i) = 0.0;
            b[k] -= factor*b[i];
        }
    }
    return 0;
}

template<typename Iter1, typename Iter2>
inline double dot(Iter1 a_begin, Iter1 a_end, Iter2 b_begin) {
    double sum = 0.0;
    while (a_begin != a_end) {
        sum += (*a_begin) * (*b_begin);
        ++a_begin; ++b_begin;
    }
    return sum;
}

template<typename InIt, typename OutIt>
inline void scale(InIt first, InIt last, OutIt result, double a) {
    for (; first != last; ++first, ++result) {
        *result = (*first) * a;
    }
}

} // namespace YamlConvector2
