#include <bits/stdc++.h>
#include "polynomial.h"
#include "rational.h"
#include "permutation.h"

using namespace std;

namespace linal {
    template <typename T>
    class Matrix {
      private:
        int n, m;
        vector <vector <T>> a;

      public:
        Matrix(int _n = 0, int _m = -1) : n(_n), m(_m) {
            if (m == -1)
                m = n;     // if n = m
            a.resize(n, vector <T> (m));
        }

        Matrix(vector <vector <T>> _a) : a(_a), n(_a.size()), m(_a[0].size()) {}


        // {lines, columns}
        pair<int, int> size() const {
            return {n, m};
        }

        // operators

        vector <T>& operator[](int i) {
            return a[i];
        }

        vector <T> operator[](int i) const {
            return a[i];
        }

        Matrix<T> operator*(const Matrix<T>& other) const {
            Matrix<T> ans(n, other.size().second);
            for (size_t i = 0; i != n; ++i)
                for (size_t j = 0; j != other.size().second; ++j)
                    for (size_t k = 0; k != m; ++k)
                        ans[i][j] += a[i][k] * other[k][j];
            return ans;
        }

        Matrix<T>& operator=(const Matrix& other) {
            if (this == &other)
                return *this;
            n = other.n, m = other.m;
            a = other.a;
            return *this;
        }

        // works if pow >= 0
        Matrix<T> operator^(size_t pow) const {
            Matrix<T> ans(*this);
            if (pow == 0) {
                for (size_t i = 0; i != n; ++i)
                    for (size_t j = 0; j != n; ++j)
                        ans[i][j] = i == j ? 1 : 0;
            } else  {
                for (size_t i = 1; i != pow; ++i)
                    ans = ans * (*this);
            }
            return ans;
        }

        Matrix<T> operator-(const Matrix& other) const {
            Matrix<T> ans(n);
            for (size_t i = 0; i != n; ++i)
                for (size_t j = 0; j != n; ++j)
                    ans[i][j] = a[i][j] - other[i][j];
            return ans;
        }

        void append_down(Matrix other) {
            n += other.size().first;
            for (size_t i = 0; i != other.size().first; ++i)
                a.push_back(other[i]);
        }

        void append_right(Matrix other) {
            m += other.size().second;
            for (size_t i = 0; i != other.size().first; ++i)
                for (size_t j = 0; j != other.size().second; ++j)
                    a[i].push_back(other[i][j]);
        }

        Matrix<T> transpose() const {
            Matrix<T> ans(m, n);
            for (size_t i = 0; i != n; ++i)
                for (size_t j = 0; j != m; ++j)
                    ans[j][i] = (*this)[i][j];
            return ans;
        }

        // By definition
        // Works in O(n!)
        T det() const {
            T ans = 0;
            Permutation p(n);
            do {
                T sum = 1;
                for (size_t i = 0; i != n; ++i) {
                    sum *= a[i][p[i]];
                }
                ans += p.sign() * sum;
            } while (p.next_perm());
            return ans;
        }

        // ^{-1}
        Matrix<T> inverse() const {
            T current_det = (*this).det();
            Matrix<T> ans(n);
            for (size_t i = 0; i != n; ++i) {
                for (size_t j = 0; j != n; ++j) {
                    Matrix<T> b(n - 1);
                    for (size_t k = 0; k != i; ++k)
                        for (size_t l = 0; l != j; ++l)
                            b[k][l] = a[k][l];
                    for (size_t k = 0; k != i; ++k)
                        for (size_t l = j + 1; l != n; ++l)
                            b[k][l - 1] = a[k][l];
                    for (size_t k = i + 1; k != n; ++k)
                        for (size_t l = 0; l != j; ++l)
                            b[k - 1][l] = a[k][l];
                    for (size_t k = i + 1; k != n; ++k)
                        for (size_t l = j + 1; l != n; ++l)
                            b[k - 1][l - 1] = a[k][l];
                    ans[j][i] = b.det();
                    if ((i + j) % 2)
                        ans[j][i] *= -1;
                    ans[j][i] /= current_det;
                }
            }
            return ans;
        }

        // by definition
        // O(n!) - same complexity as det()
        Polynomial<T> characteristic_polynomial() const {
            vector <T> temp = {0, -1};
            Matrix <Polynomial <T>> b(n);
            for (size_t i = 0; i != n; ++i) {
                for (size_t j = 0; j != n; ++j) {
                    if (i != j)
                        b[i][j] = Polynomial<T>(a[i][j]);
                    else
                        b[i][j] = Polynomial<T>(vector<T>{a[i][j], -1});
                }
            }
            return b.det();
        }

        // gaussian elimination, works in O(n^3)
        Matrix<T> gauss() const {
            Matrix<T> ans(size().first, size().second);
            for (size_t i = 0; i != size().first; ++i)
                for (size_t j = 0; j != size().second; ++j)
                    ans[i][j] = Rational((*this)[i][j]);
            for (size_t i = 0, main_colomn = 0; i != ans.size().first && main_colomn != ans.size().second; ++i, ++main_colomn) {
                for (size_t j = i; j != ans.size().first; ++j) {
                    if (ans[j][main_colomn] != static_cast<T>(0)) {
                        swap(ans[i], ans[j]);
                        break;
                    }
                }
                if (ans[i][main_colomn] == static_cast<T>(0)) {
                    --i;
                    continue;
                }
                Rational k = ans[i][main_colomn];
                for (size_t j = main_colomn; j != ans.size().second; ++j)
                    ans[i][j] /= k;
                for (size_t j = 0; j != ans.size().first; ++j) {
                    Rational k = ans[j][main_colomn];
                    for (size_t l = 0; l != ans.size().second; ++l) {
                        if (j != i)
                            ans[j][l] -= ans[i][l] * k;
                    }
                }
            }
            return ans;
        }

        // checks if A is a solution to (this * x = 0)
        bool is_solution(const Matrix& x) const {
            Matrix<T> res;
            res = *this * x;
            for (size_t i = 0; i != res.size().first; ++i)
                for (size_t j = 0; j != res.size().second; ++j)
                    if (res[i][j] != 0)
                        return false;
            return true;
        }

        // fundemental system of solutions
        Matrix<T> ker() const {
            Matrix<T> gaussed = this -> gauss();
            vector <int> main(gaussed.size().second, -1);   // if the variable is main, gives its line, otherwise -1
            for (size_t i = 0; i != gaussed.size().first; ++i) {
                for (size_t j = 0; j != gaussed.size().second; ++j) {
                    if (gaussed[i][j] == 1) {
                        main[j] = i;
                        break;
                    }
                }
            }
            Matrix<T> ans(0, gaussed.size().second);
            for (size_t j = 0; j != gaussed.size().second; ++j) {
                if (main[j] != -1)
                    continue;
                Matrix<T> cur(1, gaussed.size().second);
                cur[0][j] = 1;
                for (size_t i = 0; i != j; ++i)
                    if (main[i] != -1)
                        cur[0][i] = -gaussed[main[i]][j];
                ans.append_down(cur);
            }
            return ans.transpose();
        }

        // applies gaussian elimination and erases zero lines
        Matrix<T> im() const {
            Matrix<T> gaussed = (this -> transpose()).gauss();
            Matrix<T> ans(0, gaussed.size().second);
            for (size_t i = 0; i != gaussed.size().first; ++i) {
                Matrix<T> cur_line(1, gaussed.size().second);
                bool to_append = false;
                for (size_t j = 0; j != gaussed.size().second; ++j) {
                    cur_line[0][j] = gaussed[i][j];
                    if (gaussed[i][j] != 0)
                        to_append = true;
                }
                if (to_append)
                    ans.append_down(cur_line);
            }
            return ans.transpose();
        }

        size_t rk() const {
            Matrix<T> gaussed = this -> gauss();
            size_t ans = 0;
            for (size_t i = 0; i < gaussed.size().first; ++i) {
                bool non_empty = false;
                for (size_t j = 0; j < gaussed.size().second; ++j)
                    if (gaussed[i][j] != static_cast<T>(0))
                        non_empty = true;
                ans += non_empty;
            }
            return ans;
        }

        // from position pos append J(x, sz)
        void append_jordan_matrix(T x, size_t sz, size_t pos) {
            for (size_t i = pos; i != pos + sz; ++i) {
                a[i][i] = x;
                if (i != pos)
                    a[i - 1][i] = 1; 
            }
        }

        // jordan normal form
        Matrix<T> jnf(vector<T> eigenvalues) const {
            size_t last_free = 0;
            Matrix<T> ans(n);
            for (T& x : eigenvalues) {
                vector<size_t> ranks(n + 2);
                Matrix<T> cur = *this;
                for (size_t i = 0; i != n; ++i)
                    cur[i][i] -= x;
                for (size_t i = 0; i != n + 2; ++i)
                    ranks[i] = (cur ^ i).rk();
                for (size_t i = 1; i != n + 1; ++i) {
                    size_t sz = ranks[i - 1] - 2 * ranks[i] + ranks[i + 1];
                    for (size_t j = 0; j != sz; ++j) {
                        ans.append_jordan_matrix(x, i, last_free);
                        last_free += i;
                    }
                }
            }
            return ans;
        }

        Matrix<T> eigenvectors(vector <T> eigenvalues) const {
            Matrix<T> ans(n, 0);
            for (T& x : eigenvalues) {
                Matrix<T> cur = *this;
                for (size_t i = 0; i != n; ++i)
                    cur[i][i] -= x;
                ans.append_right(cur.ker());
            }
            return ans;
        }
    };

    template<typename T>
    Matrix<T> sum(Matrix<T> a, Matrix<T> b) {
        auto c = a;
        c.append_right(b);
        return c.im();
    }

    template <typename T>
    ostream& operator<<(ostream& out, Matrix<T> a) {
        for (size_t i = 0; i != a.size().first; ++i)
            for (size_t j = 0; j != a.size().second; ++j)
                out << a[i][j] << "\t\n"[j == a.size().second - 1];
        return out;
    }
}