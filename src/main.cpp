#include <array>
// #include <cmath>
// #include <cassert>
// #include <iostream>
#include <cstdint>
#include <limits.h>

using namespace std;

// const parameters
constexpr uint32_t N = 4;
constexpr uint32_t M = 16;
constexpr int number_of_bits = 2; // log2(N)
constexpr int N_half = M / 4;

// Define the maximum and minimum values for the fixed-point number
// instead of these, we use <limits.h>
// const int MAX_VALUE = (1 << (32 - FRACTION_BITS - 1)) - 1;
// const int MIN_VALUE = -(1 << (32 - FRACTION_BITS - 1));
// #define INT_MAX 0x7fffffff;
// #define INT_MIN 0x80000000;

typedef uint32_t uint;

/// Integer
uint sign_int(uint x) {
    return (x >= INT_MIN) ? -1 : 1;
}

uint abs_int(uint x) {
    return x * sign_int(x);
}

/// if condition is true, return a, else return b
int mux(int cond, int a, int b) {
    return (cond * a) + ((1 - cond) * b);
}

/// Real number (Fixed-point)

// Define the number of fractional bits
constexpr uint32_t FRACTION_BITS = 10; // max 13

struct FixedPoint {
    uint sign;
    uint value;
};

FixedPoint fp_from(uint x) {
    FixedPoint res;
    res.sign = sign_int(x);
    res.value = res.sign * (x << FRACTION_BITS);
    return res;
}

FixedPoint fp_from_sign(uint sign, uint x) {
    FixedPoint res;
    res.sign = sign;
    res.value = x << FRACTION_BITS;
    return res;
}

FixedPoint fp_from_direct(uint x) {
    FixedPoint res;
    res.sign = sign_int(x);
    res.value = res.sign * x;
    return res;
}

FixedPoint fp_from_direct_sign(uint sign, uint32_t x) {
    FixedPoint res;
    res.sign = sign;
    res.value = x;
    return res;
}

int fp_to_int(FixedPoint x) {
    return x.sign * (x.value >> FRACTION_BITS);
}

FixedPoint fp_add(FixedPoint self, FixedPoint other) {
    uint sum = self.sign * self.value + other.sign * other.value;
    return fp_from_direct(sum);
}

FixedPoint fp_sub(FixedPoint self, FixedPoint other) {
    uint diff = self.sign * self.value - other.sign * other.value;
    return fp_from_direct(diff);
}

FixedPoint fp_mul(FixedPoint self, FixedPoint other) {
    uint64_t product = (self.value * other.value);
    product >>= FRACTION_BITS;

    return fp_from_direct_sign(self.sign * other.sign, static_cast<uint>(product));
}

FixedPoint fp_div(FixedPoint self, FixedPoint other) {
    if (other.value == 0) {
        // throw "Divide by zero";
        return fp_from_direct_sign(1, INT_MAX);
    }
    uint64_t dividend = static_cast<uint64_t>(self.value) << FRACTION_BITS;
    uint64_t quotient = dividend / other.value;
    return fp_from_direct_sign(self.sign * other.sign, static_cast<uint>(quotient));
}

FixedPoint mux_fp(int cond, FixedPoint a, FixedPoint b) {
    // return fp_from_direct(mux(cond, a.sign, b.sign), mux(cond, a.value, b.value));
    return cond ? a : b;
}

/*
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
*/

// # of iteration for taylor series 
const int PRECISION = 10;

// Complex number
struct cpx {
    FixedPoint real;
    FixedPoint imag;
};

cpx cpx_from(uint real, uint imag) {
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

cpx cpx_mul_scalar(cpx self, FixedPoint scalar) {
    return cpx_from_fp(fp_mul(self.real, scalar), fp_mul(self.imag, scalar));
}

cpx cpx_div(cpx self, cpx other) {
    FixedPoint denominator = fp_add(fp_mul(other.real, other.real), fp_mul(other.imag, other.imag));
    if (denominator.value == 0) {
        // throw "Division by zero";
        return cpx_from(INT_MAX, INT_MAX);
    }

    // (a+bi)(c-di) = (ac+bd) + (bc-ad)i
    FixedPoint real = fp_add(fp_mul(self.real, other.real), fp_mul(self.imag, other.imag));
    FixedPoint imag = fp_sub(fp_mul(self.imag, other.real), fp_mul(self.real, other.imag));
    return cpx_from_fp(fp_div(real, denominator), fp_div(imag, denominator));
}

cpx cpx_div_scalar(cpx self, FixedPoint scalar) {
    return cpx_from_fp(fp_div(self.real, scalar), fp_div(self.imag, scalar));
}

cpx mux_cpx(int cond, cpx a, cpx b) {
    // return cpx_from_fp(mux_fp(cond, a.real, b.real), mux_fp(cond, a.imag, b.imag));
    return cond ? a : b;
}
/*
cpx cpx_powi(cpx self, int exp) {
    cpx res = cpx_from(1, 0);
    cpx id = cpx_from(1, 0);

    for (int i = 0; i < 4096; i++) {
        cpx term = mux_cpx(i < exp, self, id);
        res = cpx_mul(res, term);
    }
    return res;
}
*/

// typedef FixedPoint R;
#define R FixedPoint
#define C cpx


R sin_taylor(R x) {
    R result = x;
    R term = x;
    for (int n = 1; n <= PRECISION; n++) {
        term = fp_mul(term, fp_mul(x, x));
        term = fp_div(term, fp_from_sign(-1, (2*n)*(2*n+1)));
        
        result = fp_add(result, term);
    }
    return result;
}

R cos_taylor(R x) {
    R result = fp_from_sign(1, 1);
    R term = fp_from_sign(1, 1);
    for (int n = 1; n <= PRECISION; n++) {
        term = fp_mul(term, fp_mul(x, x));
        term = fp_div(term, fp_from_sign(-1, (2*n-1)*(2*n)));
        
        result = fp_add(result, term);
    }
    return result;
}

#define PI 3.14159265358979323846
constexpr uint32_t FP_PI = 2 * PI * (1 << FRACTION_BITS);
constexpr uint32_t M_PI = FP_PI / M;
// const R TWO_PI = fp_from_direct(1, FP_PI);

C find_mth_root_of_unity() { // m is constant
    // R ang = fp_div(TWO_PI, fp_from(m));
    R ang = fp_from_direct_sign(1, M_PI);
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

uint32_t bit_reverse_value(uint32_t x) {
    uint32_t r = 0;
    for (int i = 0; i < number_of_bits; i++) {
        r = r * 2;
        r += x % 2;
        x = x / 2;
    }
    return r;
}

template<size_t N>
array<C, N> get_bit_reverse_array(const array<C, N> a) {
    array<C, N> r = {};
    for (int i = 0; i < N; i++) {
        uint32_t bit_reverse_i = bit_reverse_value(i);
        r[bit_reverse_i] = a[i];
    }

    return r;
}


// template<size_t N>
array<C, M + 1> get_psi_powers() {
    // m^th primitive root of unity
    C psi = find_mth_root_of_unity();
    // powers of m^th primitive root of unity
    array<C, M + 1> psi_powers = {};

    C curr = cpx_from(1, 0);
    for (int i = 0; i <= M; i++) {
        psi_powers[i] = curr;
        curr = cpx_mul(curr, psi);
    }

    return psi_powers;
}

template<size_t N>
array<uint32_t, N> get_rot_group() {
    uint32_t p = 1;
    array<uint32_t, N> rot_group = {};
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
    array<C, N> b = get_bit_reverse_array(a);
    array<C, M + 1> psi_powers = get_psi_powers(); // 9, 8
    array<uint32_t, N> rot_group = get_rot_group<N>(); // 2, 8
    
    uint32_t length_n = 2;
    // while (length_n <= N) {
    for (int l = 0; l < number_of_bits; l++) {
        // for (int i = 0; i < N; i += length_n) {
        for (int i = 0; i < N; i++) {
            uint32_t lenh = length_n >> 1;
            uint32_t lenq = length_n << 2;
            uint32_t gap = M / lenq;
            // for (int j = 0; j < lenh; j++) {
            for (int j = 0; j < N / 2; j++) {
                uint32_t cond = (i % length_n == 0) && (j < lenh);
                // int index1 = cond ? i + j : 0;
                uint32_t index1 = cond * (i + j);
                uint32_t index2 = cond * (i + j + lenh);

                uint32_t idx = (rot_group[j] % lenq) * gap;
                C u = b[index1];
                C v = b[index2];
                v = cpx_mul(v, psi_powers[idx]);
                b[index1] = cond ? cpx_add(u, v) : b[index1];
                b[index2] = cond ? cpx_sub(u, v) : b[index2];
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
    array<C, M+1> psi_powers = get_psi_powers();
    array<uint32_t, N> rot_group = get_rot_group<N>();

    uint32_t length_n = N;
    // while (length_n >= 1) {
    for (int l = 0; l < number_of_bits; l++) {
        // for (int i = 0; i < N; i += length_n) {
        for (int i = 0; i < N; i++) {
            uint32_t lenh = length_n >> 1;
            uint32_t lenq = length_n << 2;
            uint32_t gap = M / lenq;
            // for (int j = 0; j < lenh; j++) {
            for (int j = 0; j < N / 2; j++) {
                uint32_t cond = (i % length_n == 0) && (j < lenh);
                uint32_t index1 = cond * (i + j);
                uint32_t index2 = cond * (i + j + lenh);

                uint32_t idx = (lenq - (rot_group[j] % lenq)) * gap;
                C u = cpx_add(b[index1], b[index2]);
                C v = cpx_sub(b[index1], b[index2]);
                v = cpx_mul(v, psi_powers[idx]);
                b[index1] = cond ? u : b[index1];
                b[index2] = cond ? v : b[index2];
            }
        }
        length_n >>= 1;
    }
    // cout << cnt << "\n";

    b = get_bit_reverse_array(b);

    // multiply by 1/n and return
    R n = fp_from_sign(1, N);
    for (int i = 0; i < N; i++) {
        b[i] = cpx_div_scalar(b[i], n);
    }

    return b;
}

template<size_t N>
array<int, 2*N> complex2int(array<C, N> a) {
    array<int, 2*N> r = {};
    for (int i = 0; i < N; i++) {
        r[2*i] = fp_to_int(a[i].real);
        r[2*i+1] = fp_to_int(a[i].imag);
        // std::cout << r[2*i] << " " << r[2*i+1] << std::endl;
    }

    return r;
}

/// Circuit
[[circuit]] auto main_circuit([[private_input]] array<int, N> a_input, int n, int m) {
// [[circuit]] array<int, 2*N> main_circuit(int n, int m) {
    // array<int, N> a_input = {1, 2, 3, 4};

    __builtin_assigner_exit_check(a_input.size() == n);
    __builtin_assigner_exit_check(n == N);
    __builtin_assigner_exit_check(m == M);
    __builtin_assigner_exit_check(is_power_of_two(n));

    array<cpx, N> a = {};
    for (int i = 0; i < N; i++) {
        a[i] = cpx_from(a_input[i], 0);
        // std::cout << a[i].real.value << " " << a[i].imag.value << std::endl;
    }

    array<C, N> b = specialFFT<N>(a);
    array<C, N> c = specialIFFT<N>(b);
    array<int, 2*N> output = complex2int<N>(c);
    // array<int, 2*N> output = complex2int<N>(a);

    return output;
}
