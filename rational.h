#include <numeric>

int gcd(int a, int b) {
    if (!a)
        return b;
    return gcd(b % a, a);
}

// Class that implemets fractions for gaussian elimination
class Rational {
  private:
    int a, b;

  public:
    int numerator() const {
        return a;
    }

    int denominator() const {
        return b;
    }

    // makes a fraction irreducible.
    void normalize() {
        int g = gcd(a, b);
        a /= g, b /= g;
        if (b < 0) {
            a *= -1;
            b *= -1;
        }
    }

    // constructors

    Rational(int _a = 0, int _b = 1) : a(_a), b(_b) {
        normalize();
    }

    Rational(const Rational& other) : a(other.numerator()),
                                      b(other.denominator()) {}
};

// operators

Rational operator+(const Rational& a, const Rational& b) {
    return {a.numerator() * b.denominator() + b.numerator() * a.denominator(),
            a.denominator() * b.denominator()};
}

template<typename T, typename U>
Rational operator+(const T& a, const U& b) {
    return Rational(a) + Rational(b);
}

Rational operator-(const Rational& a, const Rational& b) {
    return {a.numerator() * b.denominator() - b.numerator() * a.denominator(),
            a.denominator() * b.denominator()};
}

template<typename T, typename U>
Rational operator-(const T& a, const U& b) {
    return Rational(a) - Rational(b);
}

Rational operator*(const Rational& a, const Rational& b) {
    return {a.numerator() * b.numerator(),
            a.denominator() * b.denominator()};
}

template<typename T, typename U>
Rational operator*(const T& a, const U& b) {
    return Rational(a) * Rational(b);
}

Rational operator/(const Rational& a, const Rational& b) {
    return {a.numerator() * b.denominator(),
            a.denominator() * b.numerator()};
}

template<typename T, typename U>
Rational operator/(const T& a, const U& b) {
    return Rational(a) / Rational(b);
}

Rational operator+(const Rational& a) {
    return a;
}

Rational operator-(const Rational& a) {
    return {-a.numerator(), a.denominator()};
}

template<typename T>
Rational& operator+=(Rational& a, const T& b) {
    a = a + b;
    return a;
}

template<typename T>
Rational& operator-=(Rational& a, const T& b) {
    a = a - b;
    return a;
}

template<typename T>
Rational& operator*=(Rational& a, const T& b) {
    a = a * b;
    return a;
}

template<typename T>
Rational& operator/=(Rational& a, const T& b) {
    a = a / b;
    return a;
}

bool operator==(const Rational& a, const Rational& b) {
    return a.numerator() == b.numerator() &&
           a.denominator() == b.denominator();
}

template<typename T, typename U>
bool operator==(const T& a, const U& b) {
    return Rational(a) == Rational(b);
}

template<typename T, typename U>
bool operator!=(const T& a, const U& b) {
    return !(a == b);
}

Rational& operator++(Rational& a) {
    a += 1;
    return a;
}

Rational operator++(Rational& a, int) {
    Rational old_a = a;
    a += 1;
    return old_a;
}

Rational& operator--(Rational& a) {
    a -= 1;
    return a;
}

Rational operator--(Rational& a, int) {
    Rational old_a = a;
    a -= 1;
    return old_a;
}

#include <iostream>

std::ostream& operator<<(std::ostream& out, Rational& a) {
    a.normalize();
    out << a.numerator();
    if (a.denominator() != 1)
       out << "/" << a.denominator();
    return out;
}
