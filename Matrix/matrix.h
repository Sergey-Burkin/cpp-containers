#include <iostream>
#include <complex>
#include <vector>
#include <assert.h>

namespace BigNumber {
    namespace FastFourierTransform {
        const int BASE = 100;
        const size_t LOG10BASE = 2;

        void fastFourierTransform(std::vector<std::complex<double>>& a, bool invert = false) {
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
        static const int BASE = 100;
        bool isNull() const {
            return (digits.empty() || (digits.size() == 1 && digits[0] == 0));
        }

        void normalize(size_t numberOfDigits = -1) {
            for (size_t i = 0; i < digits.size(); ++i) {
                if (i > numberOfDigits && 0 <= digits[i] && digits[i] < FastFourierTransform::BASE) {
                    break;
                }
                if (digits[i] >= FastFourierTransform::BASE && i + 1 == digits.size()) {
                    digits.push_back(0);
                }
                while (digits[i] < 0) {
                    --digits[i + 1];
                    digits[i] += FastFourierTransform::BASE;
                }
                digits[i + 1] += digits[i] / FastFourierTransform::BASE;
                digits[i] %= FastFourierTransform::BASE;
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

        bool canSubtract(size_t shift, const BigInteger& a) {
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

        void doSubtract(size_t i, const BigInteger& a) {
            for (size_t j = 0; j < a.size(); ++j) {
                digits[i - j] -= a.digits[a.size() - j - 1];
                size_t newI = i;
                while (digits[newI - j] < 0) {
                    digits[newI - j] += FastFourierTransform::BASE;
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
                while (canSubtract(i, a)) {
                    doSubtract(i, a);
                    ++result.digits.back();
                }
                result.digits.push_back(0);
                if (i > 0) {
                    digits[i - 1] += (FastFourierTransform::BASE * digits[i]);
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

        BigInteger() {
            digits.push_back(0);
        };

        BigInteger(int x) {
            if (x < 0) {
                negative = true;
                x *= -1;
            }
            do {
                digits.push_back(x % FastFourierTransform::BASE);
                x /= FastFourierTransform::BASE;
            } while (x != 0);
        }

        BigInteger& operator=(const BigInteger& b) = default;

        void setNegative(bool x) {
            negative = x;
        }

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
            a.digits.resize((s.length() + FastFourierTransform::LOG10BASE - 1) / FastFourierTransform::LOG10BASE);
            for (size_t i = 0; i < a.digits.size(); ++i) {
                a.digits[i] = 0;
            }
            for (size_t i = s.size(); i-- > 0;) {
                a.digits[i / FastFourierTransform::LOG10BASE] *= 10;
                a.digits[i / FastFourierTransform::LOG10BASE] += (s[s.size() - 1 - i] - '0');
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
                res += FastFourierTransform::toStringBase(digits[size() - 1 - i], i == 0);
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
            if (isNull()) { negative = false; }
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
            FastFourierTransform::fastFourierTransform(x);
            FastFourierTransform::fastFourierTransform(y);
            for (size_t i = 0; i < n; ++i) {
                x[i] *= y[i];
            }
            FastFourierTransform::fastFourierTransform(x, true);
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
                for (size_t i = 0;
                     i < (precision + 1 + FastFourierTransform::LOG10BASE - 1) / FastFourierTransform::LOG10BASE;
                     ++i) {
                    q *= FastFourierTransform::BASE;
                    std::string current = (q / denominator).toString();
                    while (current.size() < FastFourierTransform::LOG10BASE) {
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

        Rational inverted() {
            Rational result;
            result.numerator = denominator;
            result.denominator = numerator;
            result.negative = negative;
            return result;
        }

        friend std::istream& operator>>(std::istream& in, Rational& r) {
            in >> r.numerator;
            return in;
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
}

namespace residueHelper {
    template<size_t N, size_t i, bool stop>
    const bool isPrimeHelper = N % i != 0 && isPrimeHelper<N, i + 1, ((i + i) * (i + 1) > N)>;
    template<size_t N, size_t i>
    const bool isPrimeHelper<N, i, true> = true;
    template<size_t N>
    const bool isPrime = isPrimeHelper<N, 2, (2 * 2 > N)>;

    size_t pow(size_t x, size_t a, size_t mod) {
        size_t res = 1ll;
        size_t q = x;
        while (a) {
            if (a & 1ll) {
                res *= q;
                res %= mod;
            }
            q *= q;
            q %= mod;
            a >>= 1ll;
        }
        return res;
    }
}

template<size_t N>
class Residue {
private:
    size_t val = 0;

    void normalize() {
        val %= N;
    }

public:
    Residue() = default;

    Residue(int x) : val((std::div(x, N).rem + N) % N) {};

    Residue& operator+=(const Residue& a) {
        val += a.val;
        normalize();
        return *this;
    }

    Residue& operator-=(const Residue& a) {
        val += N - a.val % N;
        normalize();
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, const Residue& residue) {
        os << residue.val;
        return os;
    }

    Residue& operator*=(const Residue& a) {
        val *= a.val % N;
        normalize();
        return *this;
    }

    Residue inverted() const {
        static_assert(residueHelper::isPrime<N>, "Division is not defined for non Prime n");
        return static_cast<Residue>(residueHelper::pow(val, N - 2, N));
    }

    Residue& operator/=(const Residue& a) {
        val *= a.inverted().val;
        normalize();
        return *this;
    }

    Residue operator-() const {
        Residue ans;
        ans.val = N - val;
        return ans;
    }

    explicit operator int() const { return val; }

    bool operator==(const Residue& a) const {
        return (a.val % N + N) % N == (val % N + N) % N;
    }

    bool operator!=(const Residue& a) const {
        return (a.val % N + N) % N != (val % N + N) % N;
    }
};

template<size_t N>
Residue<N> operator+(const Residue<N>& a, const Residue<N>& b) {
    Residue<N> copy = a;
    copy += b;
    return copy;
}

template<size_t N>
Residue<N> operator-(const Residue<N>& a, const Residue<N>& b) {
    Residue<N> copy = a;
    copy -= b;
    return copy;
}

template<size_t N>
Residue<N> operator*(const Residue<N>& a, const Residue<N>& b) {
    Residue<N> copy = a;
    copy *= b;
    return copy;
}

template<size_t N>
Residue<N> operator/(const Residue<N>& a, const Residue<N>& b) {
    Residue<N> copy = a;
    copy /= b;
    return copy;
}

namespace Row {
    template<typename T, size_t M>
    class Row {
    private:
        std::vector<T> array = std::vector<T>(M);
    public:
        T& operator[](size_t index) {
            return array[index];
        }

        const T& operator[](size_t index) const {
            return array[index];
        }

        const std::vector<T>& getVector() const {
            return array;
        }

        void swap(Row& anotherRow) {
            array.swap(anotherRow.array);
        }
    };
}

template<size_t N, size_t M, typename Field = BigNumber::Rational>
class Matrix {
private:
    std::vector<Row::Row<Field, M>> a;

    void multiplyRow(size_t idx, Field alpha = static_cast<Field>(1)) {
        for (size_t j = 0; j < M; ++j) {
            a[idx][j] *= alpha;
        }
    }

    std::pair<size_t, size_t> simplified() {
        std::pair<size_t, size_t> ans = {0, 0};
        const auto ZERO = static_cast<Field>(0);
        for (size_t j = 0; j < M; ++j) {
            if (ans.first >= N) {
                break;
            }
            for (size_t i = ans.first; i < N; ++i) {
                if (a[i][j] != ZERO) {
                    ans.second += i != ans.first;
                    swapRows(ans.first, i);
                    break;
                }
            }
            if (a[ans.first][j] == ZERO) {
                continue;
            }
            for (size_t i = 0; i < N; ++i) {
                if (i == ans.first) {
                    continue;
                }
                addRow(ans.first, i, -a[i][j] / a[ans.first][j]);
            }
            ++ans.first;
        }
        return ans;
    }

public:
    friend std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                os << matrix[i][j] << ' ';
            }
            os << '\n';
        }
        return os;
    }

    Matrix() {
        a.resize(N);
        if (N == M) {
            for (size_t i = 0; i < N; ++i) {
                a[i][i] = static_cast<Field>(1);
            }
        }
    }

    template<typename T>
    Matrix(const std::vector<std::vector<T>>& vector) : Matrix() {
        for (size_t i = 0; i < std::min(N, vector.size()); ++i) {
            for (size_t j = 0; j < std::min(M, vector[i].size()); ++j) {
                a[i][j] = static_cast<Field>(vector[i][j]);
            }
        }
    }

    template<typename T>
    Matrix(const std::initializer_list<std::initializer_list<T>>& list) : Matrix() {
        for (auto x = list.begin(); x != list.end(); ++x) {
            for (auto y = x->begin(); y != x->end(); ++y) {
                a[x - list.begin()][y - x->begin()] = static_cast<Field>(*y);
            }
        }
    }

    Matrix& operator+=(const Matrix& anotherMatrix) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                a[i][j] += anotherMatrix[i][j];
            }
        }
        return *this;
    }

    Matrix& operator-=(const Matrix& anotherMatrix) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                a[i][j] -= anotherMatrix[i][j];
            }
        }
        return *this;
    }

    template<typename T>
    Matrix& operator*=(const T& x) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                a[i][j] *= x;
            }
        }
        return *this;
    }

    Field trace() {
        static_assert(N == M, "size error");
        Field ans = 0;
        for (size_t i = 0; i < N; ++i) {
            ans += a[i][i];
        }
        return ans;
    }

    Matrix<M, N, Field> transposed() const {
        Matrix<M, N, Field> ans;
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                ans[j][i] = a[i][j];
            }
        }
        return ans;
    }

    void swapRows(size_t i, size_t j) { a[i].swap(a[j]); }

    void addRow(size_t from, size_t to, Field alpha = static_cast<Field>(1)) {
        for (size_t j = 0; j < M; ++j) {
            a[to][j] += a[from][j] * alpha;
        }
    }

    size_t rank() const {
        Matrix copy = *this;
        size_t ans = copy.simplified().first;
        return ans;
    }

    Field det() const {
        static_assert(N == M, "Det is for square matrix");
        Matrix copy = *this;
        auto tmp = copy.simplified();
        auto ans = static_cast<Field>(1);
        for (size_t i = 0; i < N; ++i) {
            ans *= copy[i][i];
        }
        return (tmp.second % 2 == 0) ? (ans) : (-ans);
    }

    void invert() {
        static_assert(N == M, "size error");
        Matrix<N, 2 * N, Field> tmp;
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                tmp[i][j] = a[i][j];
            }
        }
        for (size_t i = 0; i < N; ++i) {
            tmp[i][i + N] = static_cast<Field>(1);
        }
        tmp.simplified();
        for (size_t i = 0; i < N; ++i) {
            assert(tmp[i][i] != static_cast<Field>(0));
            tmp.multiplyRow(i, tmp[i][i].inverted());
        }
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                a[i][j] = tmp[i][j + N];
            }
        }
    }

    Matrix inverted() const {
        Matrix copy = *this;
        copy.invert();
        return copy;
    }

    std::vector<Field> getRow(size_t idx) const { return a[idx].getVector(); }

    std::vector<Field> getColumn(size_t idx) const {
        std::vector<Field> ans;
        for (size_t i = 0; i < N; ++i) {
            ans.push_back(a[i][idx]);
        }
        return ans;
    }

    Row::Row<Field, M>& operator[](size_t idx) { return a[idx]; }

    const Row::Row<Field, M>& operator[](size_t idx) const { return a[idx]; }

    bool operator==(const Matrix& anotherMatrix) const {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                if (a[i][j] != anotherMatrix[i][j]) {
                    return false;
                }
            }
        }
        return true;
    }

    bool operator!=(const Matrix& anotherMatrix) const {
        return !this->operator==(anotherMatrix);
    }
};

template<size_t N, size_t M, typename Field = BigNumber::Rational>
Matrix<N, M, Field> operator+(const Matrix<N, M, Field>& x, const Matrix<N, M, Field>& y) {
    auto copy = x;
    copy += y;
    return copy;
}

template<size_t N, size_t M, typename Field = BigNumber::Rational>
Matrix<N, M, Field> operator-(const Matrix<N, M, Field>& x, const Matrix<N, M, Field>& y) {
    auto copy = x;
    copy -= y;
    return copy;
}

template<size_t N, size_t M, typename Field = BigNumber::Rational, typename T>
Matrix<N, M, Field> operator*(const Matrix<N, M, Field>& matrix, T x) {
    auto copy = matrix;
    copy *= x;
    return copy;
}

template<size_t N, size_t M, typename T, typename Field = BigNumber::Rational>
Matrix<N, M, Field> operator*(T x, const Matrix<N, M, Field>& matrix) {
    auto copy = matrix;
    copy *= x;
    return copy;
}

template<size_t N, size_t M, size_t K, typename Field = BigNumber::Rational>
Matrix<N, K, Field> operator*(const Matrix<N, M, Field>& a, const Matrix<M, K, Field>& b) {
    Matrix<N, K, Field> ans;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < K; ++j) {
            ans[i][j] = static_cast<Field>(0);
            for (size_t l = 0; l < M; ++l) {
                ans[i][j] += a[i][l] * b[l][j];
            }
        }
    }
    return ans;
}

template<size_t N, typename Field = BigNumber::Rational>
Matrix<N, N, Field>& operator*=(Matrix<N, N, Field>& a, const Matrix<N, N, Field> b) {
    a = a * b;
    return a;
}

template<size_t N, typename Field = BigNumber::Rational>
using SquareMatrix = Matrix<N, N, Field>;