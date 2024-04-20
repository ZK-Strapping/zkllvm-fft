#include <nil/crypto3/algebra/curves/pallas.hpp>

using namespace nil::crypto3::algebra::curves;
typename pallas::base_field_type::value_type pow_2(typename pallas::base_field_type::value_type a) {

    typename pallas::base_field_type::value_type res = 1;
    for (int i = 0; i < 2; ++i) {
        res *= a;
    }
    return res;
}

[[circuit]] typename pallas::base_field_type::value_type
    field_arithmetic_example(typename pallas::base_field_type::value_type a,
                             [[private_input]] typename pallas::base_field_type::value_type p,
                             typename pallas::base_field_type::value_type b) {

    typename pallas::base_field_type::value_type c = (a + b) * a + b * (a + b) * (a + b);
    const typename pallas::base_field_type::value_type constant = 0x12345678901234567890_cppui255;
    return c * c * c / (b - a) + pow_2(a) + constant + p;
}





#include <complex>
#include <vector>
#include <cmath>
#include <cassert>

using namespace std;
typedef complex<double> C;

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

vector<C> get_bit_reverse_array(const vector<C>& a, int n) {
    assert(is_power_of_two(n));

    vector<C> r(n, C(0));
    for (int i = 0; i < n; i++) {
        int bit_reverse_i = bit_reverse_value(i, number_of_bits(n));
        r[bit_reverse_i] = a[i];
    }

    return r;
}

C find_mth_root_of_unity(int m) {
    C zeta = exp(C(0, 2 * M_PI / m));
    return zeta;
}

vector<C> get_psi_powers(int m) {
    // m^th primitive root of unity
    C psi = find_mth_root_of_unity(m);
    // powers of m^th primitive root of unity
    vector<C> psi_powers(m + 1);
    for (int i = 0; i <= m; i++) {
        psi_powers[i] = pow(psi, i);
    }

    return psi_powers;
}

vector<int> get_rot_group(int N_half, int M) {
    int p = 1;
    vector<int> rot_group;
    for (int i = 0; i < N_half; i++) {
        rot_group.push_back(p);
        p *= 5;
        p %= M;
    }

    return rot_group;
}

vector<C> specialFFT(const vector<C>& a, int n, int M) {
    assert(a.size() == n);
    assert(is_power_of_two(n));

    vector<C> b = get_bit_reverse_array(a, n);
    vector<C> psi_powers = get_psi_powers(M);
    vector<int> rot_group = get_rot_group(M >> 2, M);

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
                v *= psi_powers[idx];
                b[i + j] = u + v;
                b[i + j + lenh] = u - v;
            }
        }
        length_n *= 2;
    }

    return b;
}

vector<C> specialIFFT(const vector<C>& a, int n, int M) {
    assert(a.size() == n);
    assert(is_power_of_two(n));

    vector<C> b = a;

    int length_n = n;
    vector<C> psi_powers = get_psi_powers(M);
    vector<int> rot_group = get_rot_group(M >> 2, M);

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
