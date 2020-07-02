#include <algorithm>
#include <vector>
#include <ostream>

class Permutation {
  private:
    std::vector <int> p;
    size_t n;

  public:

    // constructors

    // default is id
    Permutation(size_t _n) : n(_n) {
        p.resize(n);
        for (size_t i = 0; i != n; ++i)
            p[i] = i;
    }

    // assuming p[i] \in [1, n]
    Permutation(const std::vector <int>& _p) : p(_p), n(_p.size()) {
        // change to [0, n)
        for (int& i : p)
            i--;
    }

    size_t size() const {
        return n;
    }

    // operators

    int& operator[](int i) {
        return p[i];
    }

    int operator[](int i) const {
        return p[i];
    }

    // assuming a * b is a(b)
    Permutation operator*(const Permutation& other) const {
        Permutation ans(n);
        for (size_t i = 0; i != n; ++i)
            ans[i] = p[other[i]];
        return ans;
    }

    Permutation operator^(int pow) const {
        if (pow == -1) {
            // inverse
            Permutation ans(n);
            for (size_t i = 0; i != n; ++i)
                ans[p[i]] = i;
            return ans;
        }
        if (pow == 0)
            return Permutation(n);
        if (pow % 2)
            return (*this ^ (pow - 1)) * (*this);
        else
            return (*this ^ (pow / 2)) * (*this ^ (pow / 2));
    }

    bool operator==(const Permutation& other) const {
        for (size_t i = 0; i != n; ++i)
            if ((*this)[i] != other[i])
                return false;
        return true;
    }

    bool next_perm() {    // O(n)
        for (size_t i = n - 1; i != 0; --i) {
            if (p[i - 1] < p[i]) {
                int mn = i;
                for (size_t j = i; j != n; ++j)
                    if (p[i - 1] < p[j] && p[j] < p[mn])
                        mn = j;
                std::swap(p[mn], p[i - 1]);
                for (size_t j = i, k = n - 1; j < k; ++j, --k)
                    std::swap(p[j], p[k]);
                return true;
            }
        }
        return false;
    }

    // calculates by parity of inverses in O(n^2)
    // TODO: change to O(n) by finding out cycles
    int sign() const {
        int sign = 1;
        for (size_t i = 0; i != n; ++i)
            for (size_t j = i + 1; j != n; ++j)
                if (p[i] > p[j])
                    sign *= -1;
        return sign;
    }
};

std::ostream& operator<<(std::ostream& out, const Permutation& a) {
    // printed elements are [1, n]
    for (size_t i = 0; i != a.size(); ++i)
        out << a[i] + 1 << " ";
    return out;
}