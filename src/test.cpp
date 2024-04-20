#include <iostream>
#include <array>
#include <cmath>
#include <cassert>
#include <limits>

using namespace std;


/// Real

// Define the number of fractional bits
constexpr int FRACTION_BITS = 10; // max 13

// Define the maximum and minimum values for the fixed-point number
constexpr int MAX_VALUE = (1 << (std::numeric_limits<int>::digits - FRACTION_BITS - 1)) - 1;
constexpr int MIN_VALUE = -(1 << (std::numeric_limits<int>::digits - FRACTION_BITS - 1));

class FixedPoint {
private:
    int value;
    FixedPoint(int value = 0) : value(value) {}

public:
    FixedPoint(float f) : value(static_cast<int>(f * (1 << FRACTION_BITS))) {}
    FixedPoint(double f) : value(static_cast<int>(f * (1 << FRACTION_BITS))) {}

    static FixedPoint from(int x) {
        return FixedPoint(x * (1 << FRACTION_BITS));
    }
    static FixedPoint from(float x) {
        return FixedPoint(x);
    }


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

    FixedPoint powi(int exp) {
        FixedPoint result = FixedPoint(1);
        for (int i = 0; i < exp; i++) {
            result = result * *this;
        }
        return result;
    }
};

FixedPoint factorial(int n) {
    int result = 1;
    for (int i = 1; i <= n; i++) {
        result = result * i;
    }
    return FixedPoint::from(result);
}

typedef FixedPoint R;

// Complex
const int PRECISION = 10;

class cpx {
public:
    FixedPoint real;
    FixedPoint imag;

    cpx () : real(0.0), imag(0.0) {}
    cpx (FixedPoint real, FixedPoint imag) : real(real), imag(imag) {}
    static cpx from(int x) {
        return cpx(R::from(x), 0.0);       
    }
    static cpx from(R x) {
        return cpx(x, 0.0);
    }

    // arithmetic operations
    cpx operator+(const cpx& other) const {
        return cpx(real + other.real, imag + other.imag);
    }
    cpx operator-(const cpx& other) const {
        return cpx(real - other.real, imag - other.imag);
    }
    cpx operator*(const cpx& other) const {
        return cpx(real * other.real - imag * other.imag, real * other.imag + imag * other.real);
    }
    cpx operator/(const cpx& other) const {
        FixedPoint denominator = other.real * other.real + other.imag * other.imag;
        if (denominator == 0) {
            return cpx(FixedPoint::from(MAX_VALUE), FixedPoint::from(MAX_VALUE));
        }
        cpx numerator = *this * cpx(other.real, -other.imag);
        return cpx(numerator.real / denominator, numerator.imag / denominator);
    }

    cpx powi(int exp) {
        cpx result = cpx::from(1);
        for (int i = 0; i < exp; i++) {
            result = result * *this;
        }
        return result;
    }
};
typedef cpx C;


R sin_taylor(R x) {
    R result (0.0);
    int sign = 1;
    for (int n=0; n<PRECISION; n++) {
        R term = R::from(sign) * x.powi(2*n+1) / factorial(2*n+1);
        result = result + term;
        sign *= -1;
    }
    return result;
}

R cos_taylor(R x) {
    R result (0.0);
    int sign = 1;
    for (int n=0; n<PRECISION; n++) {
        R term = R::from(sign) * x.powi(2*n) / factorial(2*n);
        result = result + term;
        sign *= -1;
    }
    return result;
}

#define M_PI 3.14159265358979323846
C find_mth_root_of_unity(int m) {
    R ang = R::from((float)M_PI * 2 / m);
    C zeta (cos_taylor(ang), sin_taylor(ang));
    return zeta;
}



int main() {
    FixedPoint a (1.1234f);
    FixedPoint b (1.0f);
    FixedPoint c = a + b;
    FixedPoint d = a - b;
    FixedPoint e = a * b;
    FixedPoint f = a / b;

    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "a + b = " << c << std::endl;
    std::cout << "a - b = " << d << std::endl;
    std::cout << "a * b = " << e << std::endl;
    std::cout << "a / b = " << f << std::endl;

    return 0;
}