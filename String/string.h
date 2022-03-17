#include <cstring>
#include <iostream>

class String {
private:
    size_t size = 0;
    char* str = nullptr;
    size_t realSize = 0;

    void doubleSize() {
        realSize *= 2;
        char* newString = new char[realSize];
        memcpy(newString, str, size);
        delete[] str;
        str = newString;
    }

    size_t findLR(const String& substring, bool right = false) const {
        size_t ans = length();
        size_t sLen = substring.length();
        size_t start = 0;
        size_t stop = length() - sLen;
        size_t step = 1;
        if (right) {
            std::swap(start, stop);
            step *= -1;
        }

        for (size_t i = start; i != stop + step; i += step) {
            bool goodSubstr = true;
            for (size_t j = i; j < i + sLen; ++j) {
                if (str[j] != substring[j - i]) {
                    goodSubstr = false;
                    break;
                }
            }
            if (goodSubstr) {
                ans = i;
                break;
            }
        }
        return ans;
    }

public:
    String() = default;

    explicit String(size_t n) : size(n), realSize(n) { str = new char[realSize]; }

    String(const char* s) : String(strlen(s)) {
        memcpy(str, s, strlen(s));
    }

    String(char c) : String(static_cast<size_t>(1)) {
        str[0] = c;
    }

    String(size_t n, char c) : String(n) {
        memset(str, c, size);
    }

    String(const String& s) {
        realSize = size = s.length();
        str = new char[realSize];
        memcpy(str, s.str, size);
    }

    void swap(String& s) {
        std::swap(str, s.str);
        std::swap(size, s.size);
        std::swap(realSize, s.realSize);
    }

    String& operator=(String s) {
        swap(s);
        return *this;
    }

    ~String() { delete[] str; }

    void push_back(char c) {
        if (realSize == 0) {
            size = realSize = 1;
            str = new char[1];
            str[0] = c;
            return;
        }
        if (size == realSize) {
            doubleSize();
        }
        str[size] = c;
        ++size;
    }

    void pop_back() {
        if (size > 0) {
            --size;
        }
    }

    char& front() { return str[0]; }

    const char& front() const { return str[0]; }

    char& back() { return str[size - 1]; }

    const char& back() const { return str[size - 1]; }

    char& operator[](size_t i) { return str[i]; }

    const char& operator[](size_t i) const { return str[i]; }

    String& operator+=(char c) {
        push_back(c);
        return *this;
    }

    String& operator+=(const String& s) {
        if (realSize < size + s.size) {
            realSize = 2 * (size + s.size);
            char* newString = new char[realSize];
            memcpy(newString, str, size);
            delete[] str;
            str = newString;
        }
        for (size_t i = size; i < size + s.size; ++i) {
            str[i] = s[i - size];
        }
        size += s.size;
        return *this;
    }

    size_t length() const { return size; }

    size_t find(const String& substring) const {
        return findLR(substring, false);
    }

    size_t rfind(const String& substring) const {
        return findLR(substring, true);
    }

    String substr(size_t start, size_t count) const {
        String ans(count);
        memcpy(ans.str, str + start, count);
        return ans;
    }

    bool empty() const { return (size == 0); }

    void clear() {
        size = 0;
    }

    friend std::istream& operator>>(std::istream& in, String& string) {
        char buff = ' ';
        string.clear();
        while (std::isspace(buff)) {
            buff = static_cast<char>(in.get());
        }
        while (!std::isspace(buff)) {
            if (buff == -1) {
                return in;
            }
            string.push_back(buff);
            buff = static_cast<char>(in.get());
        }
        in.putback(buff);
        return in;
    }
};

bool operator==(const String& a, const String& b) {
    if (a.length() != b.length()) {
        return false;
    }
    for (size_t i = 0; i < a.length(); ++i) {
        if (a[i] != b[i]) {
            return false;
        }
    }
    return true;
}

std::ostream& operator<<(std::ostream& out, const String& string) {
    for (size_t i = 0; i < string.length(); ++i) {
        out << string[i];
    }
    return out;
}

String operator+(const String& s1, const String& s2) {
    String res = s1;
    res += s2;
    return res;
}

