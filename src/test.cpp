#include <iostream>
#include <cmath>
#include <limits>

// Define the number of fractional bits
constexpr int FRACTION_BITS = 10; // max 13

// Define the maximum and minimum values for the fixed-point number
constexpr int MAX_VALUE = (1 << (std::numeric_limits<int>::digits - FRACTION_BITS - 1)) - 1;
constexpr int MIN_VALUE = -(1 << (std::numeric_limits<int>::digits - FRACTION_BITS - 1));

class FixedPoint {
private:
    int value;

public:
    FixedPoint(int value = 0) : value(value) {}
    FixedPoint(float f) : value(static_cast<int>(f * (1 << FRACTION_BITS))) {}

    operator float() const {
        return static_cast<float>(value) / (1 << FRACTION_BITS);
    }

    FixedPoint operator+(const FixedPoint& other) const {
        int sum = value + other.value;
        if (sum > MAX_VALUE) {
            sum = MAX_VALUE;
        } else if (sum < MIN_VALUE) {
            sum = MIN_VALUE;
        }
        return FixedPoint(sum);
    }

    FixedPoint operator-(const FixedPoint& other) const {
        int diff = value - other.value;
        if (diff > MAX_VALUE) {
            diff = MAX_VALUE;
        } else if (diff < MIN_VALUE) {
            diff = MIN_VALUE;
        }
        return FixedPoint(diff);
    }

    FixedPoint operator*(const FixedPoint& other) const {
        long long product = static_cast<long long>(value) * other.value;
        product >>= FRACTION_BITS;
        if (product > MAX_VALUE) {
            product = MAX_VALUE;
        } else if (product < MIN_VALUE) {
            product = MIN_VALUE;
        }
        return FixedPoint(static_cast<int>(product));
    }

    FixedPoint operator/(const FixedPoint& other) const {
        if (other.value == 0) {
            return FixedPoint(MAX_VALUE);
        }
        long long dividend = static_cast<long long>(value) << FRACTION_BITS;
        long long quotient = dividend / other.value;
        if (quotient > MAX_VALUE) {
            quotient = MAX_VALUE;
        } else if (quotient < MIN_VALUE) {
            quotient = MIN_VALUE;
        }
        return FixedPoint(static_cast<int>(quotient));
    }
};


int main() {
    FixedPoint a (18.3874f);
    FixedPoint b (47.139438f);
    FixedPoint c = a + b;
    FixedPoint d = a - b;
    FixedPoint e = a * b;
    FixedPoint f = a / b;
    FixedPoint g = a << 2;

    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "a + b = " << c << std::endl;
    std::cout << "a - b = " << d << std::endl;
    std::cout << "a * b = " << e << std::endl;
    std::cout << "a / b = " << f << std::endl;

    return 0;
}