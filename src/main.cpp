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


bool is_power_of_two(int n) {
    return (n & (n - 1)) == 0;
}

int number_of_bits(int n) {
    int count = 0;
    while (n > 0) {
        n >>= 1;
        count++;
    }
    return count;
}

int bit_reverse_value(int x, int n) {
    int r = 0;
    for (int i = 0; i < n; i++) {
        r <<= 1;
        r |= x & 1;
        x >>= 1;
    }
    return r;
}

template<size_t N>
array<C, N> get_bit_reverse_array(const array<C, N>& a, int n) {
    assert(is_power_of_two(n));

    array<C, N> r = {};
    for (int i = 0; i < n; i++) {
        int bit_reverse_i = bit_reverse_value(i, number_of_bits(n));
        r[bit_reverse_i] = a[i];
    }

    return r;
}



template<size_t N>
array<C, N> get_psi_powers(int m) {
    // m^th primitive root of unity
    C psi = find_mth_root_of_unity(m);
    // powers of m^th primitive root of unity
    array<C, N> psi_powers = {};
    for (int i = 0; i <= m; i++) {
        psi_powers[i] = psi.powi(i);
    }

    return psi_powers;
}

template<size_t N>
array<int, N> get_rot_group(int N_half, int M) {
    int p = 1;
    array<int, N> rot_group = {};
    for (int i = 0; i < N_half; i++) {
        rot_group[i] = p;
        p *= 5;
        p %= M;
    }

    return rot_group;
}


const int N = 4;
[[circuit]] array<C, N> specialFFT(array<int, N> a_input, int n, int M) {
    assert(a_input.size() == n);
    assert(n == N);
    assert(is_power_of_two(n));

    array<C, N> a;
    for (int i = 0; i < n; i++) {
        a[i] = C::from(a_input[i]);
    }

    array<C, N> b = get_bit_reverse_array(a, n);
    array<C, N> psi_powers = get_psi_powers<N>(M);
    array<int, N> rot_group = get_rot_group<N>(M >> 2, M);

    int length_n = 2;
    while (length_n <= n) {
        for (int i = 0; i < n; i += length_n) {
            int lenh = length_n >> 1;
            int lenq = length_n << 2;
            int gap = M / lenq;
            for (int j = 0; j < lenh; j++) {
                int idx = (rot_group[j] % lenq) * gap;
                C u = b[i + j];
                C v = b[i + j + lenh];
                v = v * psi_powers[idx];
                b[i + j] = u + v;
                b[i + j + lenh] = u - v;
            }
        }
        length_n *= 2;
    }

    return b;
}

template<size_t N>
array<C, N> specialIFFT(const array<C, N>& a, int n, int M) {
    assert(a.size() == n);
    assert(is_power_of_two(n));

    array<C, N> b = a;

    int length_n = n;
    array<C, N> psi_powers = get_psi_powers<N>(M);
    array<int, N> rot_group = get_rot_group<N>(M >> 2, M);

    while (length_n >= 1) {
        for (int i = 0; i < n; i += length_n) {
            int lenh = length_n >> 1;
            int lenq = length_n << 2;
            int gap = M / lenq;
            for (int j = 0; j < lenh; j++) {
                int idx = (lenq - (rot_group[j] % lenq)) * gap;
                C u = b[i + j] + b[i + j + lenh];
                C v = b[i + j] - b[i + j + lenh];
                v = v * psi_powers[idx];
                b[i + j] = u;
                b[i + j + lenh] = v;
            }
        }
        length_n >>= 1;
    }

    b = get_bit_reverse_array(b, n);

    // multiply by 1/n and return
    for (int i = 0; i < n; i++) {
        b[i] /= n;
    }

    return b;
}
