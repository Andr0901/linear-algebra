// Class for polynomials for characteristic polynomial

template<typename T>
class Polynomial {
  private:
    std::vector<T> a;

  public:
    // constructors

    Polynomial(const std::vector<T>& _a) : a(_a) {
        ReDegree();
    }

    Polynomial(const T& _a = T()) : a({_a}) {
        ReDegree();
    }

    template<typename It>
    Polynomial(const It& l, const It& r) {
        for (It i = l; i != r; ++i)
            a.push_back(*i);
        ReDegree();
    }

    std::vector<T> get() const {
        return a;
    }

    size_t size() const {
        return a.size();
    }

    // operators

    T& operator[](size_t i) {
        if (i >= a.size())
            a.resize(i + 1);
        return a[i];
    }

    const T operator[](size_t i) const {
        if (i >= a.size())
            return static_cast<T>(0);
        return a[i];
    }

    T operator()(const T& x) const {
        T ans = static_cast<T>(0), pow = static_cast<T>(1);
        for (size_t i = 0; i != a.size(); ++i) {
            ans += (*this)[i] * pow;
            pow *= x;
        }
        return ans;
    }

    int Degree() const {
        if (a.size() == 0)
            return -1;
        return a.size() - 1;
    }

    // iterators

    typename std::vector<T>::const_iterator begin() const {
        return a.cbegin();
    }

    typename std::vector<T>::const_iterator end() const {
        return a.cend();
    }

    // removes leading zeroes.
    void ReDegree() {
        while (!a.empty() && (a.back() == static_cast<T>(0)))
            a.pop_back();
    }

    // makes leading coefficient 1.
    void Normalize() {
        if (Degree() >= 0 && (*this)[Degree()] != static_cast<T> (0)) {
            T k = (*this)[Degree()];
            for (int i = 0; i <= Degree(); ++i)
               (*this)[i] /= k;
        }
    }
};

// operators

template<typename T>
bool operator==(const Polynomial<T>& a, const Polynomial<T>& b) {
    return a.get() == b.get();
}

template<typename T>
bool operator==(const Polynomial<T>& a, const T& b) {
    return a == Polynomial<T>(b);
}

template<typename T>
bool operator==(const T& a, const Polynomial<T>& b) {
    return Polynomial<T>(a) == b;
}

template<typename T>
bool operator!=(const Polynomial<T>& a, const Polynomial<T>& b) {
    return a.get() != b.get();
}

template<typename T>
bool operator!=(const Polynomial<T>& a, const T& b) {
    return a != Polynomial<T>(b);
}

template<typename T>
bool operator!=(const T& a, const Polynomial<T>& b) {
    return Polynomial<T>(a) != b;
}

template<typename T>
Polynomial<T> operator+(const Polynomial<T>& a, const Polynomial<T>& b) {
    std::vector<T> ans(std::max(a.Degree(), b.Degree()) + 1);
    for (size_t i = 0; i != ans.size(); ++i) {
        if (a.size() <= i)
            ans[i] = b[i];
        else if (b.size() <= i)
            ans[i] = a[i];
        else
            ans[i] = a[i] + b[i];
    }
    Polynomial<T> res(ans);
    res.ReDegree();
    return res;
}

template<typename T>
Polynomial<T> operator+(const Polynomial<T>& a, const T& b) {
    return a + Polynomial<T>(b);
}

template<typename T>
Polynomial<T> operator+(const T& a, const Polynomial<T>& b) {
    return Polynomial<T>(a) + b;
}

template<typename T>
Polynomial<T>& operator+=(Polynomial<T>& a, const Polynomial<T>& b) {
    a = a + b;
    return a;
}

template<typename T>
Polynomial<T>& operator+=(Polynomial<T>& a, const T& b) {
    a = a + b;
    return a;
}

template<typename T>
Polynomial<T> operator-(const Polynomial<T>& a, const Polynomial<T>& b) {
    std::vector<T> ans(std::max(a.Degree(), b.Degree()) + 1);
    for (size_t i = 0; i != ans.size(); ++i) {
        if (a.size() <= i)
            ans[i] = -b[i];
        else if (b.size() <= i)
            ans[i] = a[i];
        else
            ans[i] = a[i] - b[i];
    }
    Polynomial<T> res(ans);
    res.ReDegree();
    return res;
}


template<typename T>
Polynomial<T> operator-(const Polynomial<T>& a, const T& b) {
    return a - Polynomial<T>(b);
}

template<typename T>
Polynomial<T> operator-(const T& a, const Polynomial<T>& b) {
    return Polynomial<T>(a) - b;
}

template<typename T>
Polynomial<T>& operator-=(Polynomial<T>& a, const Polynomial<T>& b) {
    a = a - b;
    return a;
}

template<typename T>
Polynomial<T>& operator-=(Polynomial<T>& a, const T& b) {
    a = a - b;
    return a;
}

template<typename T>
Polynomial<T> operator*(const Polynomial<T>& a, const Polynomial<T>& b) {
    std::vector<T> ans(a.Degree() + b.Degree() + 1, static_cast<T> (0));
    for (size_t i = 0; i != a.size(); ++i)
        for (size_t j = 0; j != b.size(); ++j)
            ans[i + j] += a[i] * b[j];
    Polynomial<T> res(ans);
    res.ReDegree();
    return res;
}

template<typename T>
Polynomial<T> operator*(const Polynomial<T>& a, const T& b) {
    return a * Polynomial<T>(b);
}

template<typename T>
Polynomial<T> operator*(const T& a, const Polynomial<T>& b) {
    return Polynomial<T>(a) * b;
}

template<typename T>
Polynomial<T>& operator*=(Polynomial<T>& a, const Polynomial<T>& b) {
    a = a * b;
    return a;
}

template<typename T>
Polynomial<T>& operator*=(Polynomial<T>& a, const T& b) {
    a = a * b;
    return a;
}

template<typename T>
std::ostream& operator<<(std::ostream& out, Polynomial<T>& a) {
    if (a.Degree() == -1) {
        out << "0";
    } else if (a.Degree() == 0) {
        out << a[0];
    } else {
        for (int i = a.Degree(); i >= 0; --i) {
            if (a[i] == static_cast<T>(0))
                continue;
            if (i == 0) {
                if (a[i] > static_cast<T>(0))
                    out << "+";
                out << a[i];
            } else if (i == a.Degree()) {
                if (a[i] == static_cast<T>(-1))
                    out << "-";
                else if (a[i] != static_cast<T>(1))
                    out << a[i] << "*";
            }  else {
                if (a[i] > static_cast<T>(0))
                    out << "+";
                if (a[i] == static_cast<T>(-1))
                    out << "-";
                else if (a[i] != static_cast<T>(1))
                    out << a[i] << "*";
            }
            if (i > 0)
                out << "x";
            if (i > 1)
                out << "^" << i;
        }
    }
    return out;
}
