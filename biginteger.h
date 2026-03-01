#include <iomanip>
#include <iostream>
#include <charconv>
#include <vector>

class Rational;

class BigInteger {
    friend class Rational;

private:
    static constexpr int kDigitSize = 7;
    static constexpr int kBase = 10'000'000;

    bool is_negative_ = false;

    std::vector<int> digits_;

    static void karatsuba(unsigned long long [], unsigned long long [], unsigned long long [], size_t);

    void remove_leading_zeros();

    void shift_right();

public:
    BigInteger() = default;

    BigInteger(long long);

    BigInteger(const std::string_view&);

    BigInteger(const BigInteger&) = default;

    BigInteger& operator=(const BigInteger&) = default;

    BigInteger(BigInteger&&) noexcept;

    BigInteger& operator=(BigInteger&&) noexcept;

    explicit operator std::string() const;

    [[nodiscard]] std::string toString() const;

    BigInteger operator+() const;

    BigInteger operator-() const;

    friend BigInteger operator+(const BigInteger&, const BigInteger&);

    friend BigInteger operator-(const BigInteger&, const BigInteger&);

    friend BigInteger operator*(const BigInteger&, const BigInteger&);

    friend BigInteger operator/(const BigInteger&, const BigInteger&);

    friend BigInteger operator%(const BigInteger&, const BigInteger&);

    BigInteger& operator+=(const BigInteger&);

    BigInteger& operator-=(const BigInteger&);

    BigInteger& operator*=(const BigInteger&);

    BigInteger& operator/=(const BigInteger&);

    BigInteger& operator%=(const BigInteger&);

    BigInteger& operator++();

    BigInteger operator++(int);

    BigInteger& operator--();

    BigInteger operator--(int);

    bool operator==(const BigInteger&) const = default;

    friend std::strong_ordering operator<=>(const BigInteger&, const BigInteger&);

    explicit operator bool() const;

    friend std::ostream& operator<<(std::ostream&, const BigInteger&);

    static BigInteger module(const BigInteger&);

    static BigInteger gcd(const BigInteger& x, const BigInteger& y);
};

void BigInteger::karatsuba(unsigned long long first[], unsigned long long second[], unsigned long long result[],
                           size_t size) {
    constexpr int kSizeBruteForce = 4;

    if (size <= kSizeBruteForce) {
        for (size_t i = 0; i < size; ++i) {
            for (size_t j = 0; j < size; ++j) {
                result[i + j] += first[i] * second[j];
            }
        }

    } else if (size % 2 == 0) {
        size_t half = size / 2;
        unsigned long long left[half];
        unsigned long long right[half];
        unsigned long long temp[size];
        std::fill(temp, temp + size, 0);
        for (size_t i = 0; i < half; ++i) {
            left[i] = first[i] + first[half + i];
            right[i] = second[i] + second[half + i];
        }
        karatsuba(left, right, temp, half);

        karatsuba(first, second, result, half);
        karatsuba(first + half, second + half, result + size, half);

        unsigned long long* temp_1 = temp;
        unsigned long long* temp_2 = temp + half;
        unsigned long long* result_1 = result;
        unsigned long long* result_2 = result + half;
        unsigned long long* result_3 = result + size;
        unsigned long long* result_4 = result + 3 * half;

        for (size_t i = 0; i < half; ++i) {
            unsigned long long new_result_2 = result_2[i] + temp_1[i] - result_1[i] - result_3[i];
            unsigned long long new_result_3 = result_3[i] + temp_2[i] - result_2[i] - result_4[i];
            result[half + i] = new_result_2;
            result[size + i] = new_result_3;
        }

    } else {
        karatsuba(first, second, result, size - 1);
        for (size_t i = 0; i < size - 1; ++i) {
            result[i + size - 1] += first[i] * second[size - 1];
            result[i + size - 1] += second[i] * first[size - 1];
        }
        result[2 * size - 2] += first[size - 1] * second[size - 1];
    }
}

void BigInteger::remove_leading_zeros() {
    while (digits_.size() > 1 && digits_.back() == 0) {
        digits_.pop_back();
    }
    if (digits_.size() == 1 && digits_[0] == 0) {
        is_negative_ = false;
    }
}

void BigInteger::shift_right() {
    digits_.emplace_back(digits_.back());
    for (size_t i = digits_.size() - 2; i > 0; --i) {
        digits_[i] = digits_[i - 1];
    }
    digits_[0] = 0;
}

BigInteger::BigInteger(long long x) : is_negative_(x < 0) {
    if (x == 0) {
        digits_.emplace_back(0);
    }
    x = std::abs(x);
    while (x > 0) {
        digits_.emplace_back(x % kBase);
        x /= kBase;
    }
}

BigInteger::BigInteger(const std::string_view& str) {
    if (str.empty()) {
        is_negative_ = false;
        digits_.emplace_back(0);
        return;
    }

    if (str[0] == '-') {
        is_negative_ = true;
    }
    int check = str[0] == '-' || str[0] == '+';
    for (ssize_t i = str.length(); i > check; i -= kDigitSize) {
        digits_.emplace_back(0);
        std::string_view sub;
        if (i < kDigitSize) {
            sub = str.substr(check, i - check);
        } else {
            sub = str.substr(i - kDigitSize + (i == kDigitSize && check), kDigitSize - (i == kDigitSize && check));
        }
        std::from_chars(sub.data(), sub.data() + sub.size(), digits_.back());
    }
    remove_leading_zeros();
}

BigInteger::BigInteger(BigInteger&& bi) noexcept: is_negative_(bi.is_negative_), digits_(std::move(bi.digits_)) {
    bi.digits_ = {};
}

BigInteger& BigInteger::operator=(BigInteger&& bi) noexcept {
    is_negative_ = bi.is_negative_;
    digits_ = std::move(bi.digits_);
    bi.digits_ = {};
    return *this;
}

BigInteger::operator std::string() const {
    std::string res = (is_negative_ ? "-" : "");
    res += std::to_string(digits_.back());
    for (size_t i = digits_.size() - 1; i > 0; --i) {
        std::string num = std::to_string(digits_[i - 1]);
        res += std::string(kDigitSize - num.size(), '0');
        res += num;
    }
    return res;
}

std::string BigInteger::toString() const {
    return std::string(*this);
}

BigInteger BigInteger::operator+() const {
    return *this;
}

BigInteger BigInteger::operator-() const {
    BigInteger temp = *this;
    if (temp != 0) {
        temp.is_negative_ ^= true;
    }
    return temp;
}

BigInteger operator+(const BigInteger& lhs, const BigInteger& rhs) {
    BigInteger res = lhs;
    return res += rhs;
}

BigInteger operator-(const BigInteger& lhs, const BigInteger& rhs) {
    BigInteger res = lhs;
    return res -= rhs;
}

BigInteger operator*(const BigInteger& lhs, const BigInteger& rhs) {
    BigInteger res = lhs;
    return res *= rhs;
}

BigInteger operator/(const BigInteger& dividend, const BigInteger& divisor) {
    BigInteger res = dividend;
    return res /= divisor;
}

BigInteger operator%(const BigInteger& dividend, const BigInteger& divisor) {
    BigInteger res = dividend;
    return res %= divisor;
}

BigInteger& BigInteger::operator+=(const BigInteger& bi) {
    if (bi == 0) {
        return *this;
    }
    if (is_negative_ ^ bi.is_negative_) {
        return *this -= -bi;
    }
    int carry = 0;
    for (size_t i = 0; i < std::max(digits_.size(), bi.digits_.size()) || carry != 0; ++i) {
        if (i == digits_.size()) {
            digits_.emplace_back(0);
        }
        digits_[i] += carry + (i < bi.digits_.size() ? bi.digits_[i] : 0);
        carry = digits_[i] >= BigInteger::kBase;
        if (carry != 0) {
            digits_[i] -= BigInteger::kBase;
        }
    }
    return *this;
}

BigInteger& BigInteger::operator-=(const BigInteger& bi) {
    if (is_negative_ ^ bi.is_negative_) {
        return *this += -bi;
    }
    const bool swappable = is_negative_ ^ (*this < bi);
    const BigInteger& reduced = (swappable ? bi : *this);
    const BigInteger& deducted = (swappable ? *this : bi);
    if (swappable) {
        is_negative_ ^= true;
    }
    int carry = 0;
    for (size_t i = 0; i < deducted.digits_.size() || carry != 0; ++i) {
        if (i == digits_.size()) {
            digits_.emplace_back(0);
        }
        digits_[i] = reduced.digits_[i] - carry - (i < deducted.digits_.size() ? deducted.digits_[i] : 0);
        carry = digits_[i] < 0;
        if (carry != 0) {
            digits_[i] += BigInteger::kBase;
        }
    }
    for (size_t i = digits_.size(); i < reduced.digits_.size(); ++i) {
        digits_.emplace_back(reduced.digits_[i]);
    }
    remove_leading_zeros();
    return *this;
}

BigInteger& BigInteger::operator*=(const BigInteger& bi) {
    if (*this == 0 || bi == 0) {
        return *this = 0;
    }

    size_t size = std::max(digits_.size(), bi.digits_.size());
    unsigned long long first[size], second[size], res[2 * size];
    std::fill(first, first + size, 0);
    std::fill(second, second + size, 0);
    std::fill(res, res + 2 * size, 0);
    for (size_t i = size - digits_.size(); i < size; ++i) {
        first[i] = digits_[size - i - 1];
    }
    for (size_t i = size - bi.digits_.size(); i < size; ++i) {
        second[i] = bi.digits_[size - i - 1];
    }
    BigInteger::karatsuba(first, second, res, size);

    is_negative_ ^= bi.is_negative_;
    for (size_t i = 2 * size; i > 0; --i) {
        res[i] = res[i - 1];
    }
    res[0] = 0;

    for (size_t i = 2 * size - 1; i > 0; --i) {
        if (res[i] >= BigInteger::kBase) {
            res[i - 1] += res[i] / BigInteger::kBase;
            res[i] %= BigInteger::kBase;
        }
        if (2 * size - 1 - i == digits_.size()) {
            digits_.emplace_back(0);
        }
        digits_[2 * size - 1 - i] = res[i];
    }
    if (2 * size - 1 == digits_.size()) {
        digits_.emplace_back(0);
    }
    digits_.back() = res[0];

    remove_leading_zeros();
    return *this;
}

BigInteger& BigInteger::operator/=(const BigInteger& bi) {
    if (*this == 0) {
        return *this;
    }

    BigInteger cur = 0;
    is_negative_ ^= bi.is_negative_;

    std::vector<int> res(digits_.size(), 0);

    for (size_t i = digits_.size(); i > 0; --i) {
        cur.shift_right();
        cur.digits_[0] = digits_[i - 1];
        cur.remove_leading_zeros();

        int ith = 0;
        int left = 0, right = BigInteger::kBase;
        while (right - left > 1) {
            int mid = right - (right - left) / 2;
            BigInteger temp = BigInteger::module(bi * mid);
            if (temp <= cur) {
                ith = mid;
                left = mid;
            } else {
                right = mid;
            }
        }
        res[i - 1] = ith;
        cur -= BigInteger::module(bi * ith);
    }

    digits_ = std::move(res);

    remove_leading_zeros();
    return *this;
}

BigInteger& BigInteger::operator%=(const BigInteger& bi) {
    return *this -= (*this / bi) * bi;
}

BigInteger& BigInteger::operator++() {
    return *this += 1;
}

BigInteger BigInteger::operator++(int) {
    BigInteger temp = *this;
    ++*this;
    return temp;
}

BigInteger& BigInteger::operator--() {
    return *this -= 1;
}

BigInteger BigInteger::operator--(int) {
    BigInteger temp = *this;
    --*this;
    return temp;
}

std::strong_ordering operator<=>(const BigInteger& lhs, const BigInteger& rhs) {
    bool negation = false;
    if (lhs.is_negative_) {
        if (rhs.is_negative_) {
            negation = true;
        } else {
            return std::strong_ordering::less;
        }
    } else if (rhs.is_negative_) {
        return std::strong_ordering::greater;
    }

    if (lhs.digits_.size() != rhs.digits_.size()) {
        bool compare = lhs.digits_.size() > rhs.digits_.size();
        return compare ^ negation ? std::strong_ordering::greater : std::strong_ordering::less;
    }

    for (ssize_t i = lhs.digits_.size() - 1; i >= 0; --i) {
        if (lhs.digits_[i] != rhs.digits_[i]) {
            bool compare = lhs.digits_[i] > rhs.digits_[i];
            return compare ^ negation ? std::strong_ordering::greater : std::strong_ordering::less;
        }
    }

    return std::strong_ordering::equal;
}

BigInteger::operator bool() const {
    return *this != 0;
}

std::istream& operator>>(std::istream& is, BigInteger& bi) {
    while (std::isspace(is.peek()) || is.peek() == '\n') {
        is.get();
    }
    std::string num;
    if (is.peek() == '+' || is.peek() == '-') {
        num += static_cast<char>(is.get());
    }
    while (std::isdigit(is.peek())) {
        num += static_cast<char>(is.get());
    }
    bi = BigInteger(num == "+" || num == "-" ? "0" : num);
    return is;
}

std::ostream& operator<<(std::ostream& os, const BigInteger& bi) {
    if (bi.is_negative_) os << '-';
    os << bi.digits_.back();
    char old_fill = os.fill('0');
    for (size_t i = bi.digits_.size() - 1; i > 0; --i) {
        os << std::setw(BigInteger::kDigitSize) << bi.digits_[i - 1];
    }
    os.fill(old_fill);
    return os;
}

BigInteger operator ""_bi(unsigned long long x) {
    return x;
}

BigInteger operator ""_bi(const char* c_str, size_t size) {
    return std::string_view(c_str, size);
}

BigInteger BigInteger::module(const BigInteger& bi) {
    return bi >= 0 ? bi : -bi;
}

BigInteger BigInteger::gcd(const BigInteger& x, const BigInteger& y) {
    if (x == 0) {
        return y;
    }
    return gcd(y % x, x);
}

class Rational {
private:
    BigInteger numerator_ = 0, denominator_ = 1;

    void simplify();

public:
    Rational() = default;

    Rational(int);

    Rational(BigInteger);

    Rational operator+() const;

    Rational operator-() const;

    friend Rational operator+(const Rational&, const Rational&);

    friend Rational operator-(const Rational&, const Rational&);

    friend Rational operator*(const Rational&, const Rational&);

    friend Rational operator/(const Rational&, const Rational&);

    Rational& operator+=(const Rational&);

    Rational& operator-=(const Rational&);

    Rational& operator*=(const Rational&);

    Rational& operator/=(const Rational&);

    bool operator==(const Rational&) const = default;

    friend std::strong_ordering operator<=>(const Rational&, const Rational&);

    explicit operator double() const;

    [[nodiscard]] std::string toString() const;

    [[nodiscard]] std::string asDecimal(size_t) const;
};

void Rational::simplify() {
    BigInteger reduce = BigInteger::gcd(BigInteger::module(numerator_), denominator_);
    numerator_ /= reduce;
    denominator_ /= reduce;
}

Rational::Rational(int x) : numerator_(x) {
}

Rational::Rational(BigInteger bi) : numerator_(std::move(bi)) {
}

Rational Rational::operator+() const {
    return *this;
}

Rational Rational::operator-() const {
    Rational temp = *this;
    if (numerator_ != 0) {
        temp.numerator_.is_negative_ ^= true;
    }
    return temp;
}

Rational operator+(const Rational& lhs, const Rational& rhs) {
    Rational res = lhs;
    return res += rhs;
}

Rational operator-(const Rational& lhs, const Rational& rhs) {
    Rational res = lhs;
    return res -= rhs;
}

Rational operator*(const Rational& lhs, const Rational& rhs) {
    Rational res = lhs;
    return res *= rhs;
}

Rational operator/(const Rational& dividend, const Rational& divisor) {
    Rational res = dividend;
    return res /= divisor;
}

Rational& Rational::operator+=(const Rational& r) {
    numerator_ *= r.denominator_;
    numerator_ += r.numerator_ * denominator_;
    denominator_ *= r.denominator_;
    simplify();
    return *this;
}

Rational& Rational::operator-=(const Rational& r) {
    numerator_ *= r.denominator_;
    numerator_ -= r.numerator_ * denominator_;
    denominator_ *= r.denominator_;
    simplify();
    return *this;
}

Rational& Rational::operator*=(const Rational& r) {
    if (numerator_ == 0 || r.numerator_ == 0) {
        return *this = 0;
    }
    numerator_ *= r.numerator_;
    denominator_ *= r.denominator_;
    simplify();
    return *this;
}

Rational& Rational::operator/=(const Rational& r) {
    if (numerator_ == 0) {
        return *this;
    }
    numerator_ *= r.denominator_;
    denominator_ *= r.numerator_;
    if (denominator_.is_negative_) {
        if (numerator_ != 0) {
            numerator_.is_negative_ ^= true;
        }
        denominator_.is_negative_ ^= true;
    }
    simplify();
    return *this;
}

std::strong_ordering operator<=>(const Rational& lhs, const Rational& rhs) {
    return lhs.numerator_ * rhs.denominator_ <=> rhs.numerator_ * lhs.denominator_;
}

Rational::operator double() const {
    return std::stod(asDecimal(20));
}

std::string Rational::toString() const {
    std::string res;
    res += numerator_.toString();
    if (denominator_ != 1) {
        res += "/" + denominator_.toString();
    }
    return res;
}

std::string Rational::asDecimal(size_t precision = 0) const {
    BigInteger intPart = numerator_ / denominator_;
    BigInteger fracPart = BigInteger::module(numerator_) % denominator_;
    std::string res = (intPart == 0 && numerator_ < 0 ? "-" : "");
    res += intPart.toString();
    if (precision > 0) {
        res += ".";
        fracPart *= 10;
        for (size_t i = 0; i < precision; ++i) {
            BigInteger digit = fracPart / denominator_;
            res += digit.toString();
            fracPart = fracPart % denominator_ * 10;
        }
    }
    return res;
}
