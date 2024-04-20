#include <iostream>
#include <array>
#include <cmath>
#include <cassert>
#include <limits>

using namespace std;

const float eps = 1e-2;


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
        FixedPoint result = FixedPoint(1.0);
        for (int i = 0; i < exp; i++) {
            result = result * *this;
        }
        return result;
    }

    static FixedPoint factorial(int n) {
        int result = 1;
        for (int i = 1; i <= n; i++) {
            result = result * i;
        }
        return FixedPoint::from(result);
    }

    FixedPoint abs() {
        if (*this < 0) return -*this;
        else return *this;
    }
};

typedef FixedPoint R;

// Complex
const int PRECISION = 4;

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
        R term = R::from(sign) * x.powi(2*n+1) / R::factorial(2*n+1);
        result = result + term;
        sign *= -1;
    }
    return result;
}

R cos_taylor(R x) {
    R result (0.0);
    int sign = 1;
    for (int n=0; n<PRECISION; n++) {
        R term = R::from(sign) * x.powi(2*n) / R::factorial(2*n);
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

void FpArithmeticTest();
void FpFunctionTest();
void CpxArithmeticTest();
void CpxFunctionTest();
void sinTest();
void cosTest();
void rootTest();

int main() {
    FpArithmeticTest();
    FpFunctionTest();
    CpxArithmeticTest();
    CpxFunctionTest();
    sinTest();
    cosTest();
    rootTest();
    return 0;
}

void FpArithmeticTest() {
    cout << "FpArithmeticTest start" << "\n";
    R f1 = 1.23, f2 = 4.56, f3 = 7.89;
    R res = f1 + f2 + f3;
    R diff = res - R(13.68);

    if (diff < R(0.0)) diff = -diff;
    assert(diff < 1);

    f1 = 12.34, f2 = 56.78;
    res = f1 * f2;
    diff = res - R(700.6652);
    if (diff < R(0.0)) diff = -diff;
    assert(diff < 1);
    cout << "FpArithmeticTest end" << "\n";
}

void FpFunctionTest() { // powi, factorial, abs
    cout << "FpFunctionTest start" << "\n";
    R f1 = 1.23, f2 = 4.56, f3 = -1.23;

    assert(abs(f1) == abs(f3));

    R res = f1.powi(5);
    R diff = abs(res - R(2.28886641));
    assert(diff < 1);

    res = R::factorial(7);
    assert(R(5040.0) == res);

    cout << "FpFunctionTest end" << "\n";
}

void CpxArithmeticTest() {
    cout << "CpxArithmeticTest start" << "\n";
    cpx A (1.23, 4.56), B (1.23, -2.34);
    cpx C = A + B;
    cpx diff = C - cpx(2.46, 2.22);
    assert(diff.real.abs() < 1 && diff.imag.abs() < 1);

    C = A * B;
    diff = C - cpx(12.1833, 2.7306);
    assert(diff.real.abs() < 1 && diff.imag.abs() < 1);
    cout << "CpxArithmeticTest end" << "\n";
}

void CpxFunctionTest() {
    cout << "CpxFunctionTest start" << "\n";
    cpx A (1.23, 4.56);
    cpx diff = A.powi(3) - cpx(-74.8675, -74.1223);
    assert(diff.real.abs() < 1 && diff.imag.abs() < 1);
    cout << "CpxFunctionTest end" << "\n";
}

void sinTest() {
    
}

void cosTest() {

}

void rootTest() {

}