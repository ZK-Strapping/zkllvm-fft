#include <array>
// #include <cmath>
// #include <cassert>

using namespace std;


/// Real

// Define the number of fractional bits
const int FRACTION_BITS = 10; // max 13

// Define the maximum and minimum values for the fixed-point number
const int MAX_VALUE = (1 << (32 - FRACTION_BITS - 1)) - 1;
const int MIN_VALUE = -(1 << (32 - FRACTION_BITS - 1));



class FixedPoint {
private:
public:
    int value;
    
    FixedPoint(int value = 0, bool value_alloc = 0) : value((value_alloc ? value : value * (1 << FRACTION_BITS))) {}
    static FixedPoint from(int x) {
        return FixedPoint(x);
    }

    // float not supported in zk llvm
    // FixedPoint(float f) : value(static_cast<int>(f * (1 << FRACTION_BITS))) {}
    // FixedPoint(double f) : value(static_cast<int>(f * (1 << FRACTION_BITS))) {}
    // static FixedPoint from(float x) {
    //     return FixedPoint(x);
    // }
    // operator float() const {
    //     return static_cast<float>(value) / (1 << FRACTION_BITS);
    // }

    operator int() const {
        return value / (1 << FRACTION_BITS);
    }

    FixedPoint operator+(const FixedPoint& other) const {
        int sum = value + other.value;
        if (sum > MAX_VALUE) {
            sum = MAX_VALUE;
        } else if (sum < MIN_VALUE) {
            sum = MIN_VALUE;
        }
        return FixedPoint(sum, 1);
    }

    FixedPoint operator-(const FixedPoint& other) const {
        int diff = value - other.value;
        if (diff > MAX_VALUE) {
            diff = MAX_VALUE;
        } else if (diff < MIN_VALUE) {
            diff = MIN_VALUE;
        }
        return FixedPoint(diff, 1);
    }

    FixedPoint operator*(const FixedPoint& other) const {
        long long product = static_cast<long long>(value) * other.value;
        product >>= FRACTION_BITS;
        if (product > MAX_VALUE) {
            product = MAX_VALUE;
        } else if (product < MIN_VALUE) {
            product = MIN_VALUE;
        }
        return FixedPoint(static_cast<int>(product), 1);
    }

    FixedPoint operator/(const FixedPoint& other) const {
        if (other.value == 0) {
            return FixedPoint(MAX_VALUE, 1);
        }
        int dividend = static_cast<int>(value) << FRACTION_BITS;
        if (dividend < value) {
            throw "Divide Overflow";
        }

        int quotient = dividend / other.value;
        if (quotient > MAX_VALUE) {
            quotient = MAX_VALUE;
        } else if (quotient < MIN_VALUE) {
            quotient = MIN_VALUE;
        }
        return FixedPoint(static_cast<int>(quotient), 1);
    }

    FixedPoint powi(int exp) {
        FixedPoint result = FixedPoint(1);
        for (int i = 0; i < exp; i++) {
            result = result * *this;
        }
        return result;
    }
};

FixedPoint from2fp(int x) {
    return FixedPoint(x * (1 << FRACTION_BITS));
}

int fp2int(FixedPoint x) {
    return x.value / (1 << FRACTION_BITS);
}

// FixedPoint factorial(int n) {
//     int result = 1;
//     for (int i = 1; i <= n; i++) {
//         result = result * i;
//     }
//     return FixedPoint::from(result);
// }

// typedef FixedPoint R;
#define R FixedPoint

// Complex
const int PRECISION = 10;

class cpx {
public:
    FixedPoint real;
    FixedPoint imag;

    cpx () : real(0.0), imag(0.0) {}
    cpx (FixedPoint real, FixedPoint imag) : real(real), imag(imag) {}
    /*
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
        if (fp2int(denominator) == 0) {
            throw "Division by zero";
            return cpx(MAX_VALUE, MAX_VALUE);
        }
        cpx numerator = *this * cpx(other.real, R(0)-other.imag);
        return cpx(numerator.real / denominator, numerator.imag / denominator);
    }

    cpx powi(int exp) {
        cpx result = cpx(1, 0);
        for (int i = 0; i < exp; i++) {
            result = result * *this;
        }
        return result;
    }*/
};
#define C cpx

cpx from2cpx(int x) {
    return cpx(from2fp(x), 0.0);       
}
cpx from2cpx(R x) {
    return cpx(x, 0.0);
}


/*
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

// #define M_PI 3.14159265358979323846
C find_mth_root_of_unity(int m) {
    // R ang = R::from((float)M_PI * 2 / m);
    R ang = R::from(4 * 2 / m);
    C zeta (cos_taylor(ang), sin_taylor(ang));
    return zeta;
}
*/

bool is_power_of_two(int n) {
    return (n & (n - 1)) == 0;
}

/*
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
    __builtin_assigner_exit_check(is_power_of_two(n));

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


template <size_t N>
// [[circuit]] array<C, N> specialFFT([[private_input]] array<int, N> a_input, int n, int M) {
array<C, N> specialFFT(array<C, N> a, int n, int M) {
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
    __builtin_assigner_exit_check(a.size() == n);
    __builtin_assigner_exit_check(is_power_of_two(n));

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
*/

/// Circuit
#include <iostream>
const int N = 4;

array<int, 2*N> complex2int(array<C, N>& a, int n) {
    array<int, 2*N> r = {};
    for (int i = 0; i < N; i++) {
        r[2*i] = fp2int(a[i].real);
        r[2*i+1] = fp2int(a[i].imag);
        // std::cout << r[2*i] << " " << r[2*i+1] << std::endl;
    }

    return r;
}

[[circuit]] array<int, 2*N> main_circuit(int n, int M) {
// [[circuit]] array<int, 2*N> main_circuit([[private_input]] array<int, N> a_input, int n, int M) {
    array<int, N> a_input = {1, 2, 3, 4};

    // __builtin_assigner_exit_check(a_input.size() == n);
    // __builtin_assigner_exit_check(n == N);
    // __builtin_assigner_exit_check(is_power_of_two(n));

    array<cpx, N> a = {};
    for (int i = 0; i < N; i++) {
        a[i] = from2cpx(a_input[i]);
        // std::cout << a[i].real.value << " " << a[i].imag.value << std::endl;
    }

    array<int, 2*N> output = complex2int(a, n);
    // array<C, N> res = specialFFT<N>(a, n, M);
    // array<int, 2*N> output = complex2int<N>(res);

    return output;
}
