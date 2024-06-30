#include <array>
// #include <cmath>
// #include <cassert>
#include <cstdint>
#include <iostream>

using namespace std;

// N
const int N = 4;
const int M = 8;
const int number_of_bits = 2;

/// Integer
int sign(int64_t x) {
    return (x >= 0) - (x < 0);
}

uint64_t abs_int(int64_t x) {
    return x * sign(x);
}

int64_t div_int(int64_t x, int64_t y) {
    int64_t res = abs_int(x) / abs_int(y);
    return sign(x * y) * res;
}

/// if condition is true, return a, else return b
int mux(int cond, int a, int b) {
    return (cond * a) + ((1 - cond) * b);
}

/// Real number (Fixed-point)

// Define the number of fractional bits
const int FRACTION_BITS = 10; // max 13

// Define the maximum and minimum values for the fixed-point number
const int MAX_VALUE = (1 << (32 - FRACTION_BITS - 1)) - 1;
const int MIN_VALUE = -(1 << (32 - FRACTION_BITS - 1));

// TODO: typedef int로 대체
struct FixedPoint {
    int value;
};

FixedPoint fp_from(int x) {
    FixedPoint res;
    res.value = x * (1 << FRACTION_BITS);
    return res;
}

FixedPoint fp_from_direct(int x) {
    FixedPoint res;
    res.value = x;
    return res;
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

int fp_to_int(FixedPoint x) {
    return div_int(x.value, (1 << FRACTION_BITS));
}

FixedPoint fp_add(FixedPoint self, FixedPoint other) {
    int sum = self.value + other.value;
    // if (sum > MAX_VALUE) {
    //     sum = MAX_VALUE;
    // } else if (sum < MIN_VALUE) {
    //     sum = MIN_VALUE;
    // }
    return fp_from_direct(sum);
}

FixedPoint fp_sub(FixedPoint self, FixedPoint other) {
    int diff = self.value - other.value;
    // if (diff > MAX_VALUE) {
    //     diff = MAX_VALUE;
    // } else if (diff < MIN_VALUE) {
    //     diff = MIN_VALUE;
    // }
    return fp_from_direct(diff);
}

FixedPoint fp_mul(FixedPoint self, FixedPoint other) {
    int64_t product = static_cast<int64_t>(self.value) * other.value;
    product >>= FRACTION_BITS;
    // if (product > MAX_VALUE) {
    //     product = MAX_VALUE;
    // } else if (product < MIN_VALUE) {
    //     product = MIN_VALUE;
    // }
    return fp_from_direct(static_cast<int>(product));
}

FixedPoint fp_div(FixedPoint self, FixedPoint other) {
    if (other.value == 0) {
        throw "Divide by zero";
        return fp_from_direct(MAX_VALUE);
    }
    int64_t dividend = static_cast<int64_t>(self.value) << FRACTION_BITS;
    // if (dividend < self.value) {
    //     throw "Divide Overflow";
    // }

    int quotient = div_int(dividend, other.value);
    // if (quotient > MAX_VALUE) {
    //     quotient = MAX_VALUE;
    // } else if (quotient < MIN_VALUE) {
    //     quotient = MIN_VALUE;
    // }
    return fp_from_direct(quotient);
}

FixedPoint mux_fp(int cond, FixedPoint a, FixedPoint b) {
    return fp_from_direct(mux(cond, a.value, b.value));
}

FixedPoint fp_powi(FixedPoint self, int exp) {
    FixedPoint res = fp_from(1);
    FixedPoint id = fp_from(1);

    for (int i = 0; i < 4096; i++) {
        FixedPoint term = mux_fp(i < exp, self, id);
        res = fp_mul(res, term);
    }
    return res;
}

// factorial list
array<int, 11> factorial_list = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880,
                                          3628800, };

FixedPoint factorial(size_t n) {
    return fp_from(factorial_list[n]);
}


// # of iteration for taylor series 
const int PRECISION = 10;

// Complex number
struct cpx {
    FixedPoint real;
    FixedPoint imag;
};

cpx cpx_from(int real, int imag) {
    cpx res;
    res.real = fp_from(real);
    res.imag = fp_from(imag);
    return res;   
}

cpx cpx_from_fp(FixedPoint real, FixedPoint imag) {
    cpx res;
    res.real = real;
    res.imag = imag;
    return res;
}

// arithmetic operations
cpx cpx_add(cpx self, cpx other) {
    return cpx_from_fp(fp_add(self.real, other.real), fp_add(self.imag, other.imag));
}

cpx cpx_sub(cpx self, cpx other) {
    return cpx_from_fp(fp_sub(self.real, other.real), fp_sub(self.imag, other.imag));
}

cpx cpx_mul(cpx self, cpx other) {
    return cpx_from_fp(fp_sub(fp_mul(self.real, other.real), fp_mul(self.imag, other.imag)), 
                       fp_add(fp_mul(self.real, other.imag), fp_mul(self.imag, other.real)));
}

cpx cpx_div(cpx self, cpx other) {
    FixedPoint denominator = fp_add(fp_mul(other.real, other.real), fp_mul(other.imag, other.imag));
    if (fp_to_int(denominator) == 0) {
        throw "Division by zero";
        return cpx_from(MAX_VALUE, MAX_VALUE);
    }

    // (a+bi)(c-di) = (ac+bd) + (bc-ad)i
    FixedPoint real = fp_add(fp_mul(self.real, other.real), fp_mul(self.imag, other.imag));
    FixedPoint imag = fp_sub(fp_mul(self.imag, other.real), fp_mul(self.real, other.imag));
    return cpx_from_fp(fp_div(real, denominator), fp_div(imag, denominator));
}

cpx mux_cpx(int cond, cpx a, cpx b) {
    return cpx_from_fp(mux_fp(cond, a.real, b.real), mux_fp(cond, a.imag, b.imag));
}

cpx cpx_powi(cpx self, int exp) {
    cpx res = cpx_from(1, 0);
    cpx id = cpx_from(1, 0);

    for (int i = 0; i < 4096; i++) {
        cpx term = mux_cpx(i < exp, self, id);
        res = cpx_mul(res, term);
    }
    return res;
}

// typedef FixedPoint R;
#define R FixedPoint
#define C cpx


R sin_taylor(R x) {
    R result = fp_from_direct(0);
    int sign = 1;
    for (int n=0; n<PRECISION; n++) {
        R term = fp_mul(fp_from(sign), fp_powi(x, 2*n+1));
        term = fp_div(term, factorial(2*n+1));
        result = fp_add(result, term);
        sign *= -1;
    }
    return result;
}

R cos_taylor(R x) {
    R result = fp_from_direct(0);
    int sign = 1;
    for (int n=0; n<PRECISION; n++) {
        R term = fp_mul(fp_from(sign), fp_powi(x, 2*n));
        term = fp_div(term, factorial(2*n));
        result = fp_add(result, term);
        sign *= -1;
    }
    return result;
}

// #define M_PI 3.14159265358979323846
C find_mth_root_of_unity(int m) { // m = 8
    // R ang = R::from((float)M_PI * 2 / m);
    R ang = fp_from(1);
    C zeta = cpx_from_fp(cos_taylor(ang), sin_taylor(ang));
    return zeta;
}

bool is_power_of_two(int n) {
    return (n & (n - 1)) == 0;
}

// int number_of_bits(int n) {
//     int count = 0;
//     while (n > 0) {
//         count++;
//         n >>= 1;
//     }
//     return count;
// }

int bit_reverse_value(int x) {
    int r = 0;
    for (int i = 0; i < number_of_bits; i++) {
        r = r * 2;
        r |= x & 1;
        x = div_int(x, 2);
    }
    return r;
}

template<size_t N>
array<C, N> get_bit_reverse_array(const array<C, N> a) {

    array<C, N> r = {};
    for (int i = 0; i < N; i++) {
        int bit_reverse_i = bit_reverse_value(i);
        r[bit_reverse_i] = a[i];
    }

    return r;
}


template<size_t N>
array<C, M + 1> get_psi_powers() {
    // m^th primitive root of unity
    C psi = find_mth_root_of_unity(M);
    // powers of m^th primitive root of unity
    array<C, M + 1> psi_powers = {};
    for (int i = 0; i <= M; i++) {
        psi_powers[i] = cpx_powi(psi, i);
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
array<C, N> specialFFT(array<C, N> a) {
    cout << "pass" << "\n";
    array<C, N> b = get_bit_reverse_array(a);
    cout << "pass" << "\n";
    array<C, M + 1> psi_powers = get_psi_powers<M + 1>(); // 9, 8
    cout << "pass" << "\n";
    array<int, N> rot_group = get_rot_group<N>(M >> 2, M); // 2, 8
    cout << "pass" << "\n";
    int length_n = 2;
    while (length_n <= N) {
        for (int i = 0; i < N; i += length_n) {
            int lenh = length_n >> 1;
            int lenq = length_n << 2;
            int gap = div_int(M, lenq);
            for (int j = 0; j < lenh; j++) {
                int idx = (rot_group[j] % lenq) * gap;
                C u = b[i + j];
                C v = b[i + j + lenh];
                v = cpx_mul(v, psi_powers[idx]);
                b[i + j] = cpx_add(u, v);
                b[i + j + lenh] = cpx_sub(u, v);
            }
        }
        length_n *= 2;
    }

    return b;
}

template<size_t N>
array<C, N> specialIFFT(const array<C, N> a) {
    // __builtin_assigner_exit_check(a.size() == n);
    // __builtin_assigner_exit_check(is_power_of_two(n));

    array<C, N> b = a;

    int length_n = N;
    array<C, N> psi_powers = get_psi_powers<N>(M);
    array<int, N> rot_group = get_rot_group<N>(M >> 2, M);

    int cnt = 0;
    while (length_n >= 1) {
        cnt++;
        for (int i = 0; i < N; i += length_n) {
            int lenh = length_n >> 1;
            int lenq = length_n << 2;
            int gap = div_int(M, lenq);
            for (int j = 0; j < lenh; j++) {
                int idx = (lenq - (rot_group[j] % lenq)) * gap;
                C u = cpx_add(b[i + j], b[i + j + lenh]);
                C v = cpx_sub(b[i + j], b[i + j + lenh]);
                v = cpx_mul(v, psi_powers[idx]);
                b[i + j] = u;
                b[i + j + lenh] = v;
            }
        }
        length_n >>= 1;
    }
    cout << cnt << "\n";

    b = get_bit_reverse_array(b);

    // multiply by 1/n and return
    for (int i = 0; i < N; i++) {
        b[i] = cpx_div(b[i], cpx_from_fp(fp_from(N), fp_from(0)));
    }

    return b;
}


/// Circuit
#include <iostream>

template<size_t N>
array<int, 8> complex2int(array<C, N> a) { // N * 2 == 8
    array<int, N*2> r = {};
    for (int i = 0; i < N; i++) {
        r[2*i] = fp_to_int(a[i].real);
        r[2*i+1] = fp_to_int(a[i].imag);
        // std::cout << r[2*i] << " " << r[2*i+1] << std::endl;
    }

    return r;
}

[[circuit]] array<int, 2*N> main_circuit(int n) {
// [[circuit]] array<int, 2*N> main_circuit([[private_input]] array<int, N> a_input, int n, int M) {
    array<int, N> a_input = {1, 2, 3, 4};

    // __builtin_assigner_exit_check(a_input.size() == n);
    // __builtin_assigner_exit_check(n == N);
    // __builtin_assigner_exit_check(is_power_of_two(n));

    array<cpx, N> a = {};
    for (int i = 0; i < N; i++) {
        a[i] = cpx_from(a_input[i], 0);
        // std::cout << a[i].real.value << " " << a[i].imag.value << std::endl;
    }

    // array<C, N> res = specialFFT<N>(a);
    array<int, 2*N> output = complex2int<N>(a);

    return output;
}
