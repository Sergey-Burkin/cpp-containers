#include <cmath>
#include <iostream>
#include <complex>
#include <vector>

namespace namespaceForFft {

    const int BASE = 100;
    const size_t LOG10BASE = 2;

    void fft(std::vector<std::complex<double>>& a, bool invert = false) {
        size_t n = a.size();

        for (size_t i = 1, j = 0; i < n; ++i) {
            size_t bit = n >> 1;
            for (; j >= bit; bit >>= 1) {
                j -= bit;
            }
            j += bit;
            if (i < j) {
                swap(a[i], a[j]);
            }
        }

        for (size_t len = 2; len <= n; len <<= 1) {
            double ang = 2 * M_PI / static_cast<double>(len) * (invert ? -1 : 1);
            std::complex<double> wLen(cos(ang), sin(ang));
            for (size_t i = 0; i < n; i += len) {
                std::complex<double> w(1);
                for (size_t j = 0; j < len / 2; ++j) {
                    std::complex<double> u = a[i + j], v = a[i + j + len / 2] * w;
                    a[i + j] = u + v;
                    a[i + j + len / 2] = u - v;
                    w *= wLen;
                }
            }
        }
    }

    std::string toStringBase(int x, bool first = false) {
        std::string result = std::to_string(x);
        if (!first) {
            while (result.size() < LOG10BASE) {
                result.insert(result.begin(), '0');
            }
        }
        return result;
    }
}

class BigInteger {
private:
    std::vector<int> digits;
    bool negative = false;

    bool isNull() const {
        return (digits.empty() || (digits.size() == 1 && digits[0] == 0));
    }

    void normalize(size_t numberOfDigits = -1) {
        for (size_t i = 0; i < digits.size(); ++i) {
            if (i > numberOfDigits && 0 <= digits[i] && digits[i] < namespaceForFft::BASE) {
                break;
            }
            if (digits[i] >= namespaceForFft::BASE && i + 1 == digits.size()) {
                digits.push_back(0);
            }
            while (digits[i] < 0) {
                --digits[i + 1];
                digits[i] += namespaceForFft::BASE;
            }
            digits[i + 1] += digits[i] / namespaceForFft::BASE;
            digits[i] %= namespaceForFft::BASE;
        }
        while (digits.size() > 1 && digits.back() == 0) {
            digits.pop_back();
        }
        if (isNull()) {
            negative = false;
            if (digits.empty()) {
                digits.push_back(0);
            }
        }
    }

    void addAbs(const BigInteger& q) {
        while (digits.size() < q.size()) {
            digits.push_back(0);
        }
        for (size_t i = 0; i < std::min(digits.size(), q.size()); ++i) {
            digits[i] += q.digits[i];
        }
        normalize(q.size());
    }

    bool equalAbs(const BigInteger& b) const {
        if (size() != b.size()) {
            return false;
        }
        for (size_t i = 0; i < size(); ++i) {
            if (digits[i] != b.digits[i]) {
                return false;
            }
        }
        return true;
    }

    bool lessAbs(const BigInteger& q) const {
        if (digits.size() > q.size()) {
            return false;
        }
        if (digits.size() < q.size()) {
            return true;
        }
        for (size_t i = digits.size(); i-- > 0;) {
            if (digits[i] < q.digits[i]) {
                return true;
            }
            if (digits[i] > q.digits[i]) {
                return false;
            }
        }
        return false;
    }

    void subtractAbs(const BigInteger& q) {
        bool isNegative = lessAbs(q);
        while (digits.size() < q.size()) {
            digits.push_back(0);
        }
        for (size_t i = 0; i < std::min(digits.size(), q.size()); ++i) {
            digits[i] -= q.digits[i];
        }
        if (isNegative) {
            for (size_t i = 0; i < digits.size(); ++i) {
                digits[i] *= -1;
            }
        }
        normalize(q.size());
    }


    bool canSub(size_t shift, const BigInteger& a) {
        for (size_t j = 0; j < a.size(); ++j) {
            if (digits[shift - j] < a.digits[a.size() - j - 1]) {
                return false;
            }
            if (digits[shift - j] > a.digits[a.size() - j - 1]) {
                return true;
            }
        }
        return true;
    }

    void doSub(size_t i, const BigInteger& a) {
        for (size_t j = 0; j < a.size(); ++j) {
            digits[i - j] -= a.digits[a.size() - j - 1];
            size_t newI = i;
            while (digits[newI - j] < 0) {
                digits[newI - j] += namespaceForFft::BASE;
                --digits[newI - j + 1];
                ++newI;
            }
        }
    }

    BigInteger divideBy(const BigInteger& a) {
        BigInteger result;
        if (this == &a) {
            result = *this;
            *this = 0;
            return result;
        }
        if (a.isNull()) {
            return *this;
        }
        for (size_t i = size(); i-- + 1 > a.size();) {
            while (canSub(i, a)) {
                doSub(i, a);
                ++result.digits.back();
            }
            result.digits.push_back(0);
            if (i > 0) {
                digits[i - 1] += (namespaceForFft::BASE * digits[i]);
                digits[i] = 0;
            }
        }
        result.digits.pop_back();
        for (size_t i = 0; 2 * i < result.size(); ++i) {
            std::swap(result.digits[i], result.digits[result.size() - i - 1]);
        }
        result.normalize(0);
        if (!result.isNull()) {
            result.negative ^= negative ^ a.negative;
        }
        return result;
    }


public:

    bool isNegative() {
        return negative;
    }

    void setNegative(bool x) {
        negative = x;
    }


    BigInteger() {
        digits.push_back(0);
    };

    BigInteger(int x) {
        if (x < 0) {
            negative = true;
            x *= -1;
        }
        do {
            digits.push_back(x % namespaceForFft::BASE);
            x /= namespaceForFft::BASE;
        } while (x != 0);
    }

    BigInteger& operator=(const BigInteger& b) = default;

    size_t size() const { return digits.size(); }

    friend std::istream& operator>>(std::istream& in, BigInteger& a) {
        std::string s;
        in >> s;
        a.negative = false;
        if (s[0] == '+') {
            a.negative = false;
            s.erase(s.begin());
        }
        if (s[0] == '-') {
            a.negative = true;
            s.erase(s.begin());
        }
        a.digits.resize((s.length() + namespaceForFft::LOG10BASE - 1) / namespaceForFft::LOG10BASE);
        for (size_t i = 0; i < a.digits.size(); ++i) {
            a.digits[i] = 0;
        }
        for (size_t i = s.size(); i-- > 0;) {
            a.digits[i / namespaceForFft::LOG10BASE] *= 10;
            a.digits[i / namespaceForFft::LOG10BASE] += (s[s.size() - 1 - i] - '0');
        }
        a.normalize(0);
        return in;
    }

    std::string toString() const {
        std::string res;
        if (isNull()) {
            res += '0';
            return res;
        }
        if (negative) {
            res += '-';
        }
        for (size_t i = 0; i < size(); ++i) {
            res += namespaceForFft::toStringBase(digits[size() - 1 - i], i == 0);
        }
        return res;
    }

    friend std::ostream& operator<<(std::ostream& out, const BigInteger& a) {
        out << a.toString();
        return out;
    }

    BigInteger& operator+=(const BigInteger& b) {
        if (negative == b.negative) {
            addAbs(b);
            return *this;
        }
        negative = (lessAbs(b) ? b.negative : negative);
        subtractAbs(b);
        return *this;
    }

    BigInteger& operator-=(const BigInteger& b) {
        if (this == &b) {
            *this = 0;
            return *this;
        }
        negative ^= true;
        *this += b;
        negative ^= true;
        if (isNull()) {negative = false;}
        return *this;
    }

    BigInteger& operator++() {
        *this += 1;
        return *this;
    }

    BigInteger operator++(int) {
        BigInteger copy = *this;
        ++*this;
        return copy;
    }

    BigInteger& operator--() {
        *this -= 1;
        return *this;
    };

    BigInteger operator--(int) {
        BigInteger copy = *this;
        --*this;
        return copy;
    }

    BigInteger& operator*=(const BigInteger& b) {
        size_t n = 1;
        while (n < digits.size() + b.size()) {
            n *= 2;
        }
        std::vector<std::complex<double>> x(n), y(n);
        for (size_t i = 0; i < digits.size(); ++i) {
            x[i] = std::complex<double>(digits[i], 0);
        }
        for (size_t i = 0; i < b.size(); ++i) {
            y[i] = std::complex<double>(b.digits[i], 0);
        }
        namespaceForFft::fft(x);
        namespaceForFft::fft(y);
        for (size_t i = 0; i < n; ++i) {
            x[i] *= y[i];
        }
        namespaceForFft::fft(x, true);
        digits.resize(n);
        for (size_t i = 0; i < n; ++i) {
            digits[i] = round(x[i].real() / static_cast<double>(n));
        }
        negative ^= b.negative;
        normalize();
        return *this;
    }

    BigInteger operator-() const {
        BigInteger copy = *this;
        if (!copy.isNull()) {
            copy.negative ^= true;
        }
        return copy;
    }

    BigInteger& operator/=(const BigInteger& a) { return *this = divideBy(a); }

    BigInteger& operator%=(const BigInteger& a) {
        divideBy(a);
        normalize();
        return *this;
    }

    explicit operator bool() const { return !isNull(); }

    friend bool operator==(const BigInteger& a, const BigInteger& b) {
        if (a.isNull() && b.isNull()) {
            return true;
        }
        if (a.isNull() ^ b.isNull()) {
            return false;
        }
        return (a.negative == b.negative && a.equalAbs(b));
    }

    friend bool operator<(const BigInteger& a, const BigInteger& b) {
        if (a.isNull() && b.isNull()) {
            return false;
        }
        if (a.isNull()) {
            return !b.negative;
        }
        if (b.isNull()) {
            return a.negative;
        }
        if (a.negative != b.negative) {
            return a.negative ^ !b.negative;
        }
        if (!a.negative) {
            return a.lessAbs(b);
        }
        return b.lessAbs(a);
    }
};

BigInteger operator+(const BigInteger& a, const BigInteger& b) {
    BigInteger copy = a;
    return copy += b;
}

BigInteger operator-(const BigInteger& a, const BigInteger& b) {
    BigInteger copy = a;
    return copy -= b;
}

BigInteger operator*(const BigInteger& a, const BigInteger& b) {
    BigInteger copy = a;
    return copy *= b;
}

BigInteger operator/(const BigInteger& a, const BigInteger& b) {
    BigInteger copy = a;
    return copy /= b;
}

BigInteger operator%(const BigInteger& a, const BigInteger& b) {
    BigInteger copy = a;
    return copy %= b;
}


bool operator!=(const BigInteger& a, const BigInteger& b) { return !(a == b); }

bool operator>(const BigInteger& a, const BigInteger& b) { return b < a; }

bool operator<=(const BigInteger& a, const BigInteger& b) { return !(a > b); }

bool operator>=(const BigInteger& a, const BigInteger& b) { return !(a < b); }

BigInteger gcd(BigInteger a, BigInteger b) {
    while (b) {
        a %= b;
        std::swap(a, b);
    }
    return a;
}

class Rational {
private:
    BigInteger numerator = 0;
    BigInteger denominator = 1;
    bool negative = false;

    void normalize() {
        negative ^= numerator.isNegative() ^ denominator.isNegative();
        numerator.setNegative(false);
        denominator.setNegative(false);
        BigInteger d = gcd(numerator, denominator);
        numerator /= d;
        denominator /= d;
    }

    bool isNull() const { return numerator == 0; }

public:

    Rational() = default;

    Rational(int x) : numerator(x) {};

    Rational(const BigInteger& x) : numerator(x) {};

    Rational& operator*=(const Rational& x) {
        numerator *= x.numerator;
        denominator *= x.denominator;
        negative ^= x.negative;
        normalize();
        return *this;
    }

    Rational& operator/=(const Rational& x) {
        numerator *= x.denominator;
        denominator *= x.numerator;
        negative ^= x.negative;
        normalize();
        return *this;
    }

    Rational& operator-=(const Rational& x) {
        if (this == &x) {
            numerator = 0;
            denominator = 1;
            negative = false;
            return *this;
        }
        negative ^= true;
        *this += x;
        negative ^= true;
        if (isNull()) {
            negative = false;
        }
        return *this;
    }

    Rational& operator+=(const Rational& x) {
        if (this == &x) {
            numerator += numerator;
            normalize();
            return *this;
        }
        numerator *= x.denominator;
        if (negative != x.negative) {
            numerator -= denominator * x.numerator;
        } else {
            numerator += denominator * x.numerator;
        }
        denominator *= x.denominator;
        normalize();
        return *this;
    }

    std::string toString() const {
        if (numerator == 0) {
            return "0";
        }
        std::string result;
        if (negative) {
            result += '-';
        }
        result += numerator.toString();
        if (denominator != 1) {
            result += '/';
            result += denominator.toString();
        }
        return result;
    }

    friend std::ostream& operator<<(std::ostream& out, const Rational& rational) {
        out << rational.toString();
        return out;
    }

    Rational operator-() const {
        Rational copy = *this;
        if (!copy.isNull()) {
            copy.negative ^= true;
        }
        return copy;
    }

    std::string asDecimal(size_t precision = 0) const {
        if (isNull()) {
            return "0";
        }
        BigInteger q = numerator;
        std::string ans;
        if (negative) {
            ans += '-';
        }
        ans += (q / denominator).toString();
        size_t sizeInt = ans.size();
        q %= denominator;
        if (precision > 0) {
            ans += '.';
            for (size_t i = 0; i < (precision + 1 + namespaceForFft::LOG10BASE - 1) / namespaceForFft::LOG10BASE;
                 ++i) {
                q *= namespaceForFft::BASE;
                std::string current = (q / denominator).toString();
                while (current.size() < namespaceForFft::LOG10BASE) {
                    current.insert(current.begin(), '0');
                }
                ans += current;
                q %= denominator;
            }
            while (ans.size() > sizeInt + 1 + precision) {
                ans.pop_back();
            }
        }
        return ans;
    }

    explicit operator double() const {
        double ans = stod(asDecimal(50));
        return ans;
    }

    friend bool operator==(const Rational& q, const Rational& r) {
        if (q.isNull() && r.isNull()) {
            return true;
        }
        if (q.isNull() ^ r.isNull()) {
            return false;
        }
        return (q.negative == r.negative && q.numerator == r.numerator && q.denominator == r.denominator);
    }

    friend bool operator<(const Rational& q, const Rational& r) {
        if (q.isNull() && r.isNull()) {
            return false;
        }
        if (q.isNull()) {
            return !r.negative;
        }
        if (r.isNull()) {
            return q.negative;
        }
        if (q.negative && !r.negative) {
            return true;
        }
        if (!q.negative && r.negative) {
            return false;
        }
        return (q.numerator * r.denominator < q.denominator * r.numerator);
    }
};

bool operator!=(const Rational& q, const Rational& r) { return !(q == r); }

bool operator>(const Rational& x, const Rational& y) { return y < x; }

bool operator<=(const Rational& x, const Rational& y) { return !(y < x); }

bool operator>=(const Rational& x, const Rational& y) { return !(x < y); }

Rational operator+(const Rational& a, const Rational& b) {
    Rational copy = a;
    copy += b;
    return copy;
}

Rational operator-(const Rational& a, const Rational& b) {
    Rational copy = a;
    copy -= b;
    return copy;
}

Rational operator*(const Rational& a, const Rational& b) {
    Rational copy = a;
    copy *= b;
    return copy;
}

Rational operator/(const Rational& a, const Rational& b) {
    Rational copy = a;
    copy /= b;
    return copy;
}
