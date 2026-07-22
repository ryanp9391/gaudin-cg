#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>
#include <chrono>   // added: timing for the highest-weight benchmark
#include <string>   // added: benchmark driver uses std::string

// --- benchmark instrumentation (added; not part of the original algorithm) ---
// compute_highest_weight_coeffs() records its own wall-clock duration here, so
// the pure highest-weight solve can be timed in isolation from everything else.
static double g_hw_solve_seconds = 0.0;
namespace benchtimer {
    struct ScopeTimer {
        double &acc;
        std::chrono::steady_clock::time_point t0;
        ScopeTimer(double &a) : acc(a), t0(std::chrono::steady_clock::now()) {}
        ~ScopeTimer() {
            acc = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - t0).count();
        }
    };
}

// Declaration of LAPACK subroutines
// Make sure the data types match your version of LAPACK

extern "C" void dgesvd_(char const* JOBU,
                        char const* JOBVT,
                        int const* M,
                        int const* N,
                        double* A,
                        int const* LDA,
                        double* S,
                        double* U,
                        int const* LDU,
                        double* VT,
                        int const* LDVT,
                        double* WORK,
                        int const* LWORK,
                        int *INFO);

extern "C" void dgels_(char const* TRANS,
                       int const* M,
                       int const* N,
                       int const* NRHS,
                       double* A,
                       int const* LDA,
                       double* B,
                       int const* LDB,
                       double* WORK,
                       int const* LWORK,
                       int *INFO);

namespace clebsch {
    const double EPS = 1e-12;

    // binomial coefficients
    class binomial_t {
        std::vector<int> cache;
        int N;

    public:
        int operator()(int n, int k);
    } binomial;

    // Eq. (19) and (25)
    class weight {
        std::vector<int> elem;

    public:
        // the N in "SU(N)"
        const int N;

        // create a non-initialized weight
        weight(int N);

        // create irrep weight of given index
        // Eq. (C2)
        weight(int N, int index);

        // assign from another instance
        clebsch::weight &operator=(const clebsch::weight &w);

        // access elements of this weight (k = 1, ..., N)
        int &operator()(int k);
        const int &operator()(int k) const;

        // compare weights
        // Eq. (C1)
        bool operator<(const weight &w) const;
        bool operator==(const weight &w) const;

        // element-wise sum of weights
        clebsch::weight operator+(const weight &w) const;

        // returns the index of this irrep weight (index = 0, 1, ...)
        // Eq. (C2)
        int index() const;

        // returns the dimension of this irrep weight
        // Eq. (22)
        long long dimension() const;
    };

    // Eq. (20)
    class pattern {
        std::vector<int> elem;

    public:
        // the N in "SU(N)"
        const int N;

        // copy constructor
        pattern(const pattern &pat);

        // create pattern of given index from irrep weight
        // Eq. (C7)
        pattern(const weight &irrep, int index = 0);

        // access elements of this pattern (l = 1, ..., N; k = 1, ..., l)
        int &operator()(int k, int l);
        const int &operator()(int k, int l) const;

        // find succeeding/preceding pattern, return false if not possible
        // Eq. (C9)
        bool operator++();
        bool operator--();

        // returns the pattern index (index = 0, ..., dimension - 1)
        // Eq. (C7)
        int index() const;

        // returns the pattern weight
        // Eq. (25)
        clebsch::weight get_weight() const;

        // returns matrix element of lowering operator J^(l)_-
        // between this pattern minus M^(k,l) and this pattern
        // (l = 1, ..., N; k = 1, ..., l)
        // Eq. (28)
        double lowering_coeff(int k, int l) const;

        // returns matrix element of raising operator J^(l)_+
        // between this pattern plus M^(k,l) and this pattern
        // (l = 1, ..., N; k = 1, ..., l)
        // Eq. (29)
        double raising_coeff(int k, int l) const;
    };

    class decomposition {
        std::vector<clebsch::weight> weights;
        std::vector<int> multiplicities;

    public:
        // the N in "SU(N)"
        const int N;

        // save given irreps for later use
        const weight factor1, factor2;

        // construct the decomposition of factor1 times factor2 into irreps
        // Eq. (31)
        decomposition(const weight &factor1, const weight &factor2);

        // return the number of occurring irreps
        int size() const;

        // access the occurring irreps
        // j = 0, ..., size() - 1
        const clebsch::weight &operator()(int j) const;

        // return the outer multiplicity of irrep in this decomposition
        int multiplicity(const weight &irrep) const;
    };

    class index_adapter {
        std::vector<int> indices;
        std::vector<int> multiplicities;

    public:
        // the N in "SU(N)"
        const int N;

        // save given irreps for later use
        const int factor1, factor2;

        // construct this index_adapter from a given decomposition
        index_adapter(const clebsch::decomposition &decomp);

        // return the number of occurring irreps
        int size() const;

        // access the occurring irreps
        int operator()(int j) const;

        // return the outer multiplicity of irrep in this decomposition
        int multiplicity(int irrep) const;
    };

    class coefficients {
        std::map<std::vector<int>, double> clzx;

        // access Clebsch-Gordan coefficients in convenient manner
        void set(int factor1_state,
                 int factor2_state,
                 int multiplicity_index,
                 int irrep_state,
                 double value);

        // internal functions, doing most of the work
        void highest_weight_normal_form(); // Eq. (37)
        void compute_highest_weight_coeffs(); // Eq. (36)
        void compute_lower_weight_coeffs(int multip_index, int state, std::vector<char> &done); // Eq. (40)

    public:
        // the N in "SU(N)"
        const int N;

        // save irreps and their dimensions for later use
        const weight factor1, factor2, irrep;
        const int factor1_dimension, factor2_dimension, irrep_dimension;

        // outer multiplicity of irrep in this decomposition
        const int multiplicity;

        // construct all Clebsch-Gordan coefficients of this decomposition.
        // highest_weight_only = true stops after the highest-weight state
        // (skips building all descendant states) -- added for the benchmark.
        coefficients(const weight &irrep, const weight &factor1, const weight &factor2,
                     bool highest_weight_only = false);

        // access Clebsch-Gordan coefficients (read-only)
        // multiplicity_index = 0, ..., multiplicity - 1
        // factor1_state = 0, ..., factor1_dimension - 1
        // factor2_state = 0, ..., factor2_dimension - 1
        // irrep_state = 0, ..., irrep_dimension
        double operator()(int factor1_state,
                          int factor2_state,
                          int multiplicity_index,
                          int irrep_state) const;
    };
};

// implementation of "binomial_t" starts here

int clebsch::binomial_t::operator()(int n, int k) {
    if (N <= n) {
        for (cache.resize((n + 1) * (n + 2) / 2); N <= n; ++N) {
            cache[N * (N + 1) / 2] = cache[N * (N + 1) / 2 + N] = 1;
            for (int k = 1; k < N; ++k) {
                cache[N * (N + 1) / 2 + k] = cache[(N - 1) * N / 2 + k]
                                           + cache[(N - 1) * N / 2 + k - 1];
            }
        }
    }

    return cache[n * (n + 1) / 2 + k];
}

// implementation of "weight" starts here

clebsch::weight::weight(int N) : elem(N), N(N) {}

clebsch::weight::weight(int N, int index) : elem(N, 0), N(N) {
    for (int i = 0; index > 0 && i < N; ++i) {
        for (int j = 1; binomial(N - i - 1 + j, N - i - 1) <= index; j <<= 1) {
            elem[i] = j;
        }

        for (int j = elem[i] >> 1; j > 0; j >>= 1) {
            if (binomial(N - i - 1 + (elem[i] | j), N - i - 1) <= index) {
                elem[i] |= j;
            }
        }

        index -= binomial(N - i - 1 + elem[i]++, N - i - 1);
    }
}

clebsch::weight &clebsch::weight::operator=(const clebsch::weight &w) {
    int &n = const_cast<int &>(N);
    elem = w.elem;
    n = w.N;
    return *this;
}

int &clebsch::weight::operator()(int k) {
    assert(1 <= k && k <= N);
    return elem[k - 1];
}

const int &clebsch::weight::operator()(int k) const {
    assert(1 <= k && k <= N);
    return elem[k - 1];
}

bool clebsch::weight::operator<(const weight &w) const {
    assert(w.N == N);
    for (int i = 0; i < N; ++i) {
        if (elem[i] - elem[N - 1] != w.elem[i] - w.elem[N - 1]) {
            return elem[i] - elem[N - 1] < w.elem[i] - w.elem[N - 1];
        }
    }
    return false;
}

bool clebsch::weight::operator==(const weight &w) const {
    assert(w.N == N);

    for (int i = 1; i < N; ++i) {
        if (w.elem[i] - w.elem[i - 1] != elem[i] - elem[i - 1]) {
            return false;
        }
    }

    return true;
}

clebsch::weight clebsch::weight::operator+(const weight &w) const {
    weight result(N);

    transform(elem.begin(), elem.end(), w.elem.begin(), result.elem.begin(), std::plus<int>());

    return result;
}

int clebsch::weight::index() const {
    int result = 0;

    for (int i = 0; elem[i] > elem[N - 1]; ++i) {
        result += binomial(N - i - 1 + elem[i] - elem[N - 1] - 1, N - i - 1);
    }

    return result;
}

long long clebsch::weight::dimension() const {
    long long numerator = 1, denominator = 1;

    for (int i = 1; i < N; ++i) {
        for (int j = 0; i + j < N; ++j) {
            numerator *= elem[j] - elem[i + j] + i;
            denominator *= i;
        }
    }

    return numerator / denominator;
}

// implementation of "pattern" starts here

clebsch::pattern::pattern(const pattern &p) : elem(p.elem), N(p.N) {}

clebsch::pattern::pattern(const weight &irrep, int index) :
        elem((irrep.N * (irrep.N + 1)) / 2), N(irrep.N) {
    for (int i = 1; i <= N; ++i) {
        (*this)(i, N) = irrep(i);
    }

    for (int l = N - 1; l >= 1; --l) {
        for (int k = 1; k <= l; ++k) {
            (*this)(k, l) = (*this)(k + 1, l + 1);
        }
    }

    while (index-- > 0) {
        bool b = ++(*this);

        assert(b);
    }
}

int &clebsch::pattern::operator()(int k, int l) {
    return elem[(N * (N + 1) - l * (l + 1)) / 2 + k - 1];
}

const int &clebsch::pattern::operator()(int k, int l) const {
    return elem[(N * (N + 1) - l * (l + 1)) / 2 + k - 1];
}

bool clebsch::pattern::operator++() {
    int k = 1, l = 1;

    while (l < N && (*this)(k, l) == (*this)(k, l + 1)) {
        if (--k == 0) {
            k = ++l;
        }
    }

    if (l == N) {
        return false;
    }

    ++(*this)(k, l);

    while (k != 1 || l != 1) {
        if (++k > l) {
            k = 1;
            --l;
        }

        (*this)(k, l) = (*this)(k + 1, l + 1);
    }

    return true;
}

bool clebsch::pattern::operator--() {
    int k = 1, l = 1;

    while (l < N && (*this)(k, l) == (*this)(k + 1, l + 1)) {
        if (--k == 0) {
            k = ++l;
        }
    }

    if (l == N) {
        return false;
    }

    --(*this)(k, l);

    while (k != 1 || l != 1) {
        if (++k > l) {
            k = 1;
            --l;
        }

        (*this)(k, l) = (*this)(k, l + 1);
    }

    return true;
}

int clebsch::pattern::index() const {
    int result = 0;

    for (pattern p(*this); --p; ++result) {}

    return result;
}

clebsch::weight clebsch::pattern::get_weight() const {
    clebsch::weight result(N);

    for (int prev = 0, l = 1; l <= N; ++l) {
        int now = 0;

        for (int k = 1; k <= l; ++k) {
            now += (*this)(k, l);
        }

        result(l) = now - prev;
        prev = now;
    }

    return result;
}

double clebsch::pattern::lowering_coeff(int k, int l) const {
    double result = 1.0;

    for (int i = 1; i <= l + 1; ++i) {
        result *= (*this)(i, l + 1) - (*this)(k, l) + k - i + 1;
    }
    
    for (int i = 1; i <= l - 1; ++i) {
        result *= (*this)(i, l - 1) - (*this)(k, l) + k - i;
    }

    for (int i = 1; i <= l; ++i) {
        if (i == k) continue;
        result /= (*this)(i, l) - (*this)(k, l) + k - i + 1;
        result /= (*this)(i, l) - (*this)(k, l) + k - i;
    }

    return std::sqrt(-result);
}

double clebsch::pattern::raising_coeff(int k, int l) const {
    double result = 1.0;

    for (int i = 1; i <= l + 1; ++i) {
        result *= (*this)(i, l + 1) - (*this)(k, l) + k - i;
    }

    for (int i = 1; i <= l - 1; ++i) {
        result *= (*this)(i, l - 1) - (*this)(k, l) + k - i - 1;
    }

    for (int i = 1; i <= l; ++i) {
        if (i == k) continue;
        result /= (*this)(i, l) - (*this)(k, l) + k - i;
        result /= (*this)(i, l) - (*this)(k, l) + k - i  - 1;
    }

    return std::sqrt(-result);
}

// implementation of "decomposition" starts here

clebsch::decomposition::decomposition(const weight &factor1, const weight &factor2) :
        N(factor1.N), factor1(factor1), factor2(factor2) {
    assert(factor1.N == factor2.N);
    std::vector<clebsch::weight> result;
    pattern low(factor1), high(factor1);
    weight trial(factor2);
    int k = 1, l = N;

    do {
        while (k <= N) {
            --l;
            if (k <= l) {
                low(k, l) = std::max(high(k + N - l, N), high(k, l + 1) + trial(l + 1) - trial(l));
                high(k, l) = high(k, l + 1);
                if (k > 1 && high(k, l) > high(k - 1, l - 1)) {
                    high(k, l) = high(k - 1, l - 1);
                }
                if (l > 1 && k == l && high(k, l) > trial(l - 1) - trial(l)) {
                    high(k, l) = trial(l - 1) - trial(l);
                }
                if (low(k, l) > high(k, l)) {
                    break;
                }
                trial(l + 1) += high(k, l + 1) - high(k, l);
            } else {
                trial(l + 1) += high(k, l + 1);
                ++k;
                l = N;
            }
        }

        if (k > N) {
            result.push_back(trial);
            for (int i = 1; i <= N; ++i) {
                result.back()(i) -= result.back()(N);
            }
        } else {
            ++l;
        }

        while (k != 1 || l != N) {
            if (l == N) {
                l = --k - 1;
                trial(l + 1) -= high(k, l + 1);
            } else if (low(k, l) < high(k, l)) {
                --high(k, l);
                ++trial(l + 1);
                break;
            } else {
                trial(l + 1) -= high(k, l + 1) - high(k, l);
            }
            ++l;
        }
    } while (k != 1 || l != N);

    sort(result.begin(), result.end());
    for (std::vector<clebsch::weight>::iterator it = result.begin(); it != result.end(); ++it) {
        if (it != result.begin() && *it == weights.back()) {
            ++multiplicities.back();
        } else {
            weights.push_back(*it);
            multiplicities.push_back(1);
        }
    }
}

int clebsch::decomposition::size() const {
    return weights.size();
}

const clebsch::weight &clebsch::decomposition::operator()(int j) const {
    return weights[j];
}

int clebsch::decomposition::multiplicity(const weight &irrep) const {
    assert(irrep.N == N);
    std::vector<clebsch::weight>::const_iterator it
        = std::lower_bound(weights.begin(), weights.end(), irrep);

    return it != weights.end() && *it == irrep ? multiplicities[it - weights.begin()] : 0;
}

// implementation of "index_adapter" starts here

clebsch::index_adapter::index_adapter(const clebsch::decomposition &decomp) :
        N(decomp.N),
        factor1(decomp.factor1.index()),
        factor2(decomp.factor2.index()) {
    for (int i = 0, s = decomp.size(); i < s; ++i) {
        indices.push_back(decomp(i).index());
        multiplicities.push_back(decomp.multiplicity(decomp(i)));
    }
}

int clebsch::index_adapter::size() const {
    return indices.size();
}

int clebsch::index_adapter::operator()(int j) const {
    return indices[j];
}

int clebsch::index_adapter::multiplicity(int irrep) const {
    std::vector<int>::const_iterator it = std::lower_bound(indices.begin(), indices.end(), irrep);

    return it != indices.end() && *it == irrep ? multiplicities[it - indices.begin()] : 0;
}

// implementation of "clebsch" starts here

void clebsch::coefficients::set(int factor1_state,
                                int factor2_state,
                                int multiplicity_index,
                                int irrep_state,
                                double value) {
    assert(0 <= factor1_state && factor1_state < factor1_dimension);
    assert(0 <= factor2_state && factor2_state < factor2_dimension);
    assert(0 <= multiplicity_index && multiplicity_index < multiplicity);
    assert(0 <= irrep_state && irrep_state < irrep_dimension);

    int coefficient_label[] = { factor1_state,
                                factor2_state,
                                multiplicity_index,
                                irrep_state };
    clzx[std::vector<int>(coefficient_label, coefficient_label
            + sizeof coefficient_label / sizeof coefficient_label[0])] = value;
}

void clebsch::coefficients::highest_weight_normal_form() {
    int hws = irrep_dimension - 1;

    // bring CGCs into reduced row echelon form
    for (int h = 0, i = 0; h < multiplicity - 1 && i < factor1_dimension; ++i) {
        for (int j = 0; h < multiplicity - 1 && j < factor2_dimension; ++j) {
            int k0 = h;

            for (int k = h + 1; k < multiplicity; ++k) {
                if (fabs((*this)(i, j, k, hws)) > fabs((*this)(i, j, k0, hws))) {
                    k0 = k;
                }
            }

            if ((*this)(i, j, k0, hws) < -EPS) {
                for (int i2 = i; i2 < factor1_dimension; ++i2) {
                    for (int j2 = i2 == i ? j : 0; j2 < factor2_dimension; ++j2) {
                        set(i2, j2, k0, hws, -(*this)(i2, j2, k0, hws));
                    }
                }
            } else if ((*this)(i, j, k0, hws) < EPS) {
                continue;
            }

            if (k0 != h) {
                for (int i2 = i; i2 < factor1_dimension; ++i2) {
                    for (int j2 = i2 == i ? j : 0; j2 < factor2_dimension; ++j2) {
                        double x = (*this)(i2, j2, k0, hws);
                        set(i2, j2, k0, hws, (*this)(i2, j2, h, hws));
                        set(i2, j2, h, hws, x);
                    }
                }
            }

            for (int k = h + 1; k < multiplicity; ++k) {
                for (int i2 = i; i2 < factor1_dimension; ++i2) {
                    for (int j2 = i2 == i ? j : 0; j2 < factor2_dimension; ++j2) {
                        set(i2, j2, k, hws, (*this)(i2, j2, k, hws) - (*this)(i2, j2, h, hws)
                                * (*this)(i, j, k, hws) / (*this)(i, j, h, hws));
                    }
                }
            }

            // next 3 lines not strictly necessary, might improve numerical stability
            for (int k = h + 1; k < multiplicity; ++k) {
                set(i, j, k, hws, 0.0);
            }

            ++h;
        }
    }

    // Gram-Schmidt orthonormalization
    for (int h = 0; h < multiplicity; ++h) {
        for (int k = 0; k < h; ++k) {
            double overlap = 0.0;
            for (int i = 0; i < factor1_dimension; ++i) {
                for (int j = 0; j < factor2_dimension; ++j) {
                    overlap += (*this)(i, j, h, hws) * (*this)(i, j, k, hws);
                }
            }

            for (int i = 0; i < factor1_dimension; ++i) {
                for (int j = 0; j < factor2_dimension; ++j) {
                    set(i, j, h, hws, (*this)(i, j, h, hws) - overlap * (*this)(i, j, k, hws));
                }
            }
        }

        double norm = 0.0;
        for (int i = 0; i < factor1_dimension; ++i) {
            for (int j = 0; j < factor2_dimension; ++j) {
                norm += (*this)(i, j, h, hws) * (*this)(i, j, h, hws);
            }
        }
        norm = std::sqrt(norm);

        for (int i = 0; i < factor1_dimension; ++i) {
            for (int j = 0; j < factor2_dimension; ++j) {
                set(i, j, h, hws, (*this)(i, j, h, hws) / norm);
            }
        }
    }
}

void clebsch::coefficients::compute_highest_weight_coeffs() {
    benchtimer::ScopeTimer _hw_timer(g_hw_solve_seconds);  // added: time this call

    if (multiplicity == 0) {
        return;
    }

    std::vector<std::vector<int> > map_coeff(factor1_dimension,
                                             std::vector<int>(factor2_dimension, -1));
    std::vector<std::vector<int> > map_states(factor1_dimension,
                                              std::vector<int>(factor2_dimension, -1));
    int n_coeff = 0, n_states = 0;
    pattern p(factor1, 0);

    for (int i = 0; i < factor1_dimension; ++i, ++p) {
        weight pw(p.get_weight());
        pattern q(factor2, 0);
        for (int j = 0; j < factor2_dimension; ++j, ++q) {
            if (pw + q.get_weight() == irrep) {
                map_coeff[i][j] = n_coeff++;
            }
        }
    }

    if (n_coeff == 1) {
        for (int i = 0; i < factor1_dimension; ++i) {
            for (int j = 0; j < factor2_dimension; ++j) {
                if (map_coeff[i][j] >= 0) {
                    set(i, j, 0, irrep_dimension - 1, 1.0);
                    return;
                }
            }
        }
    }

    double *hw_system = new double[n_coeff * (factor1_dimension * factor2_dimension)];
    assert(hw_system != NULL);
    memset(hw_system, 0, n_coeff * (factor1_dimension * factor2_dimension) * sizeof (double));

    pattern r(factor1, 0);
    for (int i = 0; i < factor1_dimension; ++i, ++r) {
        pattern q(factor2, 0);

        for (int j = 0; j < factor2_dimension; ++j, ++q) {
            if (map_coeff[i][j] >= 0) {
                for (int l = 1; l <= N - 1; ++l) {
                    for (int k = 1; k <= l; ++k) {
                        if ((k == 1 || r(k, l) + 1 <= r(k - 1, l - 1)) && r(k, l) + 1 <= r(k, l + 1)) {
                            ++r(k, l);
                            int h = r.index();
                            --r(k, l);

                            if (map_states[h][j] < 0) {
                                map_states[h][j] = n_states++;
                            }

                            hw_system[n_coeff * map_states[h][j] + map_coeff[i][j]]
                                += r.raising_coeff(k, l);
                        }

                        if ((k == 1 || q(k, l) + 1 <= q(k - 1, l - 1)) && q(k, l) + 1 <= q(k, l + 1)) {
                            ++q(k, l);
                            int h = q.index();
                            --q(k, l);

                            if (map_states[i][h] < 0) {
                                map_states[i][h] = n_states++;
                            }


                            hw_system[n_coeff * map_states[i][h] + map_coeff[i][j]]
                                += q.raising_coeff(k, l);
                        }
                    }
                }
            }
        }
    }

    int lwork = -1, info;
    double worksize;

    double *singval = new double[std::min(n_coeff, n_states)];
    assert(singval != NULL);
    double *singvec = new double[n_coeff * n_coeff];
    assert(singvec != NULL);

    dgesvd_("A",
            "N",
            &n_coeff,
            &n_states,
            hw_system,
            &n_coeff,
            singval,
            singvec,
            &n_coeff,
            NULL,
            &n_states,
            &worksize,
            &lwork,
            &info);
    assert(info == 0);

    lwork = worksize;
    double *work = new double[lwork];
    assert(work != NULL);

    dgesvd_("A",
            "N",
            &n_coeff,
            &n_states,
            hw_system,
            &n_coeff,
            singval,
            singvec,
            &n_coeff,
            NULL,
            &n_states,
            work,
            &lwork,
            &info);
    assert(info == 0);

    for (int i = 0; i < multiplicity; ++i) {
        for (int j = 0; j < factor1_dimension; ++j) {
            for (int k = 0; k < factor2_dimension; ++k) {
                if (map_coeff[j][k] >= 0) {
                    double x = singvec[n_coeff * (n_coeff - 1 - i) + map_coeff[j][k]];

                    if (fabs(x) > EPS) {
                        set(j, k, i, irrep_dimension - 1, x);
                    }
                }
            }
        }
    }

    // uncomment next line to bring highest-weight coefficients into "normal form"
    // highest_weight_normal_form();

    delete[] work;
    delete[] singvec;
    delete[] singval;
    delete[] hw_system;
}

void clebsch::coefficients::compute_lower_weight_coeffs(int multip_index,
                                                        int state,
                                                        std::vector<char> &done) {
    weight statew(pattern(irrep, state).get_weight());
    pattern p(irrep, 0);
    std::vector<int> map_parent(irrep_dimension, -1),
                     map_multi(irrep_dimension, -1),
                     which_l(irrep_dimension, -1);
    int n_parent = 0, n_multi = 0;

    for (int i = 0; i < irrep_dimension; ++i, ++p) {
        weight v(p.get_weight());

        if (v == statew) {
            map_multi[i] = n_multi++;
        } else for (int l = 1; l < N; ++l) {
            --v(l);
            ++v(l + 1);
            if (v == statew) {
                map_parent[i] = n_parent++;
                which_l[i] = l;
                if (!done[i]) {
                    compute_lower_weight_coeffs(multip_index, i, done);
                }
                break;
            }
            --v(l + 1);
            ++v(l);
        }
    }

    double *irrep_coeffs = new double[n_parent * n_multi];
    assert(irrep_coeffs != NULL);
    memset(irrep_coeffs, 0, n_parent * n_multi * sizeof (double));

    double *prod_coeffs = new double[n_parent * factor1_dimension * factor2_dimension];
    assert(prod_coeffs != NULL);
    memset(prod_coeffs, 0, n_parent * factor1_dimension * factor2_dimension * sizeof (double));

    std::vector<std::vector<int> > map_prodstat(factor1_dimension,
                                                std::vector<int>(factor2_dimension, -1));
    int n_prodstat = 0;

    pattern r(irrep, 0);
    for (int i = 0; i < irrep_dimension; ++i, ++r) {
        if (map_parent[i] >= 0) {
            for (int k = 1, l = which_l[i]; k <= l; ++k) {
                if (r(k, l) > r(k + 1, l + 1) && (k == l || r(k, l) > r(k, l - 1))) {
                    --r(k, l);
                    int h = r.index();
                    ++r(k, l);

                    irrep_coeffs[n_parent * map_multi[h] + map_parent[i]] += r.lowering_coeff(k, l);
                }
            }

            pattern q1(factor1, 0);
            for (int j1 = 0; j1 < factor1_dimension; ++j1, ++q1) {
                pattern q2(factor2, 0);

                for (int j2 = 0; j2 < factor2_dimension; ++j2, ++q2) {
                    if (std::fabs((*this)(j1, j2, multip_index, i)) > EPS) {
                        for (int k = 1, l = which_l[i]; k <= l; ++k) {
                            if (q1(k, l) > q1(k + 1, l + 1) && (k == l || q1(k, l) > q1(k, l - 1))) {
                                --q1(k, l);
                                int h = q1.index();
                                ++q1(k, l);

                                if (map_prodstat[h][j2] < 0) {
                                    map_prodstat[h][j2] = n_prodstat++;
                                }

                                prod_coeffs[n_parent * map_prodstat[h][j2] + map_parent[i]] +=
                                        (*this)(j1, j2, multip_index, i) * q1.lowering_coeff(k, l);
                            }

                            if (q2(k, l) > q2(k + 1, l + 1) && (k == l || q2(k, l) > q2(k, l - 1))) {
                                --q2(k, l);
                                int h = q2.index();
                                ++q2(k, l);

                                if (map_prodstat[j1][h] < 0) {
                                    map_prodstat[j1][h] = n_prodstat++;
                                }

                                prod_coeffs[n_parent * map_prodstat[j1][h] + map_parent[i]] +=
                                        (*this)(j1, j2, multip_index, i) * q2.lowering_coeff(k, l);
                            }
                        }
                    }
                }
            }
        }
    }

    double worksize;
    int lwork = -1, info;

    dgels_("N",
           &n_parent,
           &n_multi,
           &n_prodstat,
           irrep_coeffs,
           &n_parent,
           prod_coeffs,
           &n_parent,
           &worksize,
           &lwork,
           &info);
    assert(info == 0);

    lwork = worksize;
    double *work = new double[lwork];
    assert(work != NULL);

    dgels_("N",
           &n_parent,
           &n_multi,
           &n_prodstat,
           irrep_coeffs,
           &n_parent,
           prod_coeffs,
           &n_parent,
           work,
           &lwork,
           &info);
    assert(info == 0);

    for (int i = 0; i < irrep_dimension; ++i) {
        if (map_multi[i] >= 0) {
            for (int j = 0; j < factor1_dimension; ++j) {
                for (int k = 0; k < factor2_dimension; ++k) {
                    if (map_prodstat[j][k] >= 0) {
                        double x = prod_coeffs[n_parent * map_prodstat[j][k] + map_multi[i]];

                        if (fabs(x) > EPS) {
                            set(j, k, multip_index, i, x);
                        }
                    }
                }
            }

            done[i] = true;
        }
    }

    delete[] work;
    delete[] prod_coeffs;
    delete[] irrep_coeffs;
}

clebsch::coefficients::coefficients(const weight &irrep, const weight &factor1, const weight &factor2,
                                    bool highest_weight_only) :
        N(irrep.N),
        factor1(factor1),
        factor2(factor2),
        irrep(irrep),
        factor1_dimension(factor1.dimension()),
        factor2_dimension(factor2.dimension()),
        irrep_dimension(irrep.dimension()),
        multiplicity(decomposition(factor1, factor2).multiplicity(irrep)) {
    assert(factor1.N == irrep.N);
    assert(factor2.N == irrep.N);

    compute_highest_weight_coeffs();

    if (highest_weight_only) {
        return;   // added: stop here -- do NOT build the descendant states
    }

    for (int i = 0; i < multiplicity; ++i) {
        std::vector<char> done(irrep_dimension, 0);
        done[irrep_dimension - 1] = true;
        for (int j = irrep_dimension - 1; j >= 0; --j) {
            if (!done[j]) {
                compute_lower_weight_coeffs(i, j, done);
            }
        }
    }
}

double clebsch::coefficients::operator()(int factor1_state,
                                         int factor2_state,
                                         int multiplicity_index,
                                         int irrep_state) const {
    assert(0 <= factor1_state && factor1_state < factor1_dimension);
    assert(0 <= factor2_state && factor2_state < factor2_dimension);
    assert(0 <= multiplicity_index && multiplicity_index < multiplicity);
    assert(0 <= irrep_state && irrep_state < irrep_dimension);

    int coefficient_label[] = { factor1_state,
                                factor2_state,
                                multiplicity_index,
                                irrep_state };
    std::map<std::vector<int>, double>::const_iterator it(
            clzx.find(std::vector<int>(coefficient_label, coefficient_label
                    + sizeof coefficient_label / sizeof coefficient_label[0])));

    return it != clzx.end() ? it->second : 0.0;
}


// ============================================================================
// Highest-weight CGC benchmark driver  (added -- replaces the interactive menu)
//
// Computes ONLY the highest-weight-state Clebsch-Gordan coefficients of a
// target irrep S'' in the product S (x) S', and reports how long that takes,
// WITHOUT building the descendant states the original program also computes.
//
// Usage:
//     cgc_hw_bench N  s1..sN   p1..pN   q1..qN   [reps] [verify]
//
//   N          the N of SU(N)
//   s1..sN     irrep S   as N integers (highest-weight / Young-row labels)
//   p1..pN     irrep S'  as N integers
//   q1..qN     irrep S'' as N integers (the target whose h.w. state we want)
//   reps       number of timed repetitions (default 11); reported headline is
//              the MINIMUM over reps (least affected by OS jitter)
//   verify     if present, also run the FULL original computation once and
//              check the highest-weight coefficients agree (slow for big reps)
//
// Example:  cgc_hw_bench 4  2 1 0 0   2 1 0 0   3 2 1 0   101
// ============================================================================

static std::string fmt_time(double s) {
    char buf[64];
    if      (s < 1e-6) std::snprintf(buf, sizeof buf, "%8.3f ns", s * 1e9);
    else if (s < 1e-3) std::snprintf(buf, sizeof buf, "%8.3f us", s * 1e6);
    else if (s < 1.0)  std::snprintf(buf, sizeof buf, "%8.3f ms", s * 1e3);
    else               std::snprintf(buf, sizeof buf, "%8.4f s ", s);
    return std::string(buf);
}

static void stats(std::vector<double> v, double &mn, double &med, double &mean) {
    std::sort(v.begin(), v.end());
    mn = v.front();
    med = v[v.size() / 2];
    double sum = 0.0;
    for (std::size_t i = 0; i < v.size(); ++i) sum += v[i];
    mean = sum / v.size();
}

static void report(const char *label, std::vector<double> &v) {
    double mn, med, mean;
    stats(v, mn, med, mean);
    std::cout << "  " << label
              << "   min " << fmt_time(mn)
              << "   median " << fmt_time(med)
              << "   mean " << fmt_time(mean) << "\n";
}

static void print_irrep(const char *name, const clebsch::weight &w, int N) {
    std::cout << name << " = (";
    for (int k = 1; k <= N; ++k) std::cout << w(k) << (k < N ? " " : "");
    std::cout << ")";
}

int main(int argc, char **argv) {
    const char *usage =
        "Usage: cgc_hw_bench N  s1..sN   p1..pN   q1..qN   [reps] [verify]\n"
        "  e.g. cgc_hw_bench 4  2 1 0 0   2 1 0 0   3 2 1 0   101\n";

    if (argc < 2) { std::cerr << usage; return 1; }

    // ---- sub-command: "decompose N s1..sN p1..pN" -> list irreps mu (mult, dim) ----
    // Fast: only runs the Littlewood-Richardson decomposition, no CGCs. One line
    // per mu:  "MU <m1> .. <mN> mult <k> dim <d>"  (machine-readable for test.sh).
    if (std::string(argv[1]) == "decompose") {
        if (argc < 3) { std::cerr << "Usage: cgc_hw_bench decompose N s1..sN p1..pN\n"; return 1; }
        int Nd = std::atoi(argv[2]);
        if (Nd < 1) { std::cerr << "N must be >= 1\n"; return 1; }
        int need_d = 3 + 2 * Nd;
        if (argc < need_d) { std::cerr << "Expected " << 2 * Nd << " labels for S and S'.\n"; return 1; }
        clebsch::weight A(Nd), B(Nd);
        int p = 3;
        for (int k = 1; k <= Nd; ++k) A(k) = std::atoi(argv[p++]);
        for (int k = 1; k <= Nd; ++k) B(k) = std::atoi(argv[p++]);
        clebsch::decomposition d(A, B);
        for (int i = 0; i < d.size(); ++i) {
            std::cout << "MU";
            for (int k = 1; k <= Nd; ++k) std::cout << ' ' << d(i)(k);
            std::cout << " mult " << d.multiplicity(d(i))
                      << " dim " << d(i).dimension() << "\n";
        }
        return 0;
    }

    int N = std::atoi(argv[1]);
    if (N < 1) { std::cerr << "N must be >= 1\n" << usage; return 1; }

    int need = 2 + 3 * N;   // argv[1]=N, then 3N integer labels
    if (argc < need) {
        std::cerr << "Expected " << 3 * N << " integer labels (N=" << N
                  << " each for S, S', S'').\n" << usage;
        return 1;
    }

    clebsch::weight S(N), Sp(N), Spp(N);
    int a = 2;
    for (int k = 1; k <= N; ++k) S(k)   = std::atoi(argv[a++]);
    for (int k = 1; k <= N; ++k) Sp(k)  = std::atoi(argv[a++]);
    for (int k = 1; k <= N; ++k) Spp(k) = std::atoi(argv[a++]);

    int reps = 11;
    bool verify = false;
    bool machine = false;
    for (int i = need; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "verify") verify = true;
        else if (arg == "machine") machine = true;
        else { int r = std::atoi(argv[i]); if (r > 0) reps = r; }
    }

    // Confirm S'' occurs in S (x) S' and get its outer multiplicity.
    clebsch::decomposition dtest(S, Sp);
    int mult = dtest.multiplicity(Spp);
    if (mult == 0) {
        std::cerr << "Target irrep S'' does not occur in S (x) S'. Irreps present:\n";
        for (int i = 0; i < dtest.size(); ++i) {
            std::cerr << "   ";
            for (int k = 1; k <= N; ++k) std::cerr << dtest(i)(k) << ' ';
            std::cerr << "(mult " << dtest.multiplicity(dtest(i)) << ")\n";
        }
        return 1;
    }

    long long dimS = S.dimension(), dimSp = Sp.dimension(), dimSpp = Spp.dimension();

    if (!machine) {
        std::cout << "SU(" << N << ")   ";
        print_irrep("S", S, N);  std::cout << " (x) ";
        print_irrep("S'", Sp, N); std::cout << " -> target ";
        print_irrep("S''", Spp, N); std::cout << "\n";
        std::cout << "dim S = " << dimS << ",  dim S' = " << dimSp
                  << ",  dim S'' = " << dimSpp
                  << ",  outer multiplicity = " << mult << "\n";
        std::cout << "timed repetitions = " << reps
                  << "  (headline = minimum)\n\n";
    }

    // Warm up caches (binomials etc.) so they don't skew the first timing.
    { clebsch::coefficients warm(Spp, S, Sp, true); (void) warm; }

    std::vector<double> t_dec, t_solve, t_total;
    double sink = 0.0;
    for (int r = 0; r < reps; ++r) {
        // (a) Littlewood-Richardson decomposition alone
        auto d0 = std::chrono::steady_clock::now();
        clebsch::decomposition dd(S, Sp);
        int m = dd.multiplicity(Spp);
        auto d1 = std::chrono::steady_clock::now();
        sink += m;
        t_dec.push_back(std::chrono::duration<double>(d1 - d0).count());

        // (b) full: labels -> highest-weight CGCs, in highest-weight-ONLY mode
        auto c0 = std::chrono::steady_clock::now();
        clebsch::coefficients C(Spp, S, Sp, true);
        auto c1 = std::chrono::steady_clock::now();
        t_total.push_back(std::chrono::duration<double>(c1 - c0).count());

        // (c) the isolated highest-weight solve, recorded inside the ctor
        t_solve.push_back(g_hw_solve_seconds);

        sink += C(0, 0, 0, (int) dimSpp - 1);   // touch a coeff (defeat DCE)
    }

    if (machine) {
        // one compact, parseable line (raw seconds = minimum over reps)
        double sMin, sMed, sMean, dMin, dMed, dMean, tMin, tMed, tMean;
        stats(t_solve, sMin, sMed, sMean);
        stats(t_dec, dMin, dMed, dMean);
        stats(t_total, tMin, tMed, tMean);
        std::printf("RESULT dim=%lld mult=%d hw_solve_s=%.9g decomp_s=%.9g total_s=%.9g\n",
                    dimSpp, mult, sMin, dMin, tMin);
        return 0;
    }

    std::cout << "TIMINGS\n";
    report("highest-weight solve only      ", t_solve);
    report("decomposition (LR) only        ", t_dec);
    report("total: labels -> h.w. CGCs     ", t_total);
    std::cout << "\n";

    // Print the actual highest-weight coefficients (from a fresh h.w.-only run)
    // and a normalization check, so the result can be trusted.
    clebsch::coefficients C(Spp, S, Sp, true);
    int hws = (int) dimSpp - 1;   // highest-weight state index Q(M'') = dim(S'')
    std::cout << "HIGHEST-WEIGHT CGCs  (Q(M''] = " << dimSpp << ")\n";
    std::cout.setf(std::ios::fixed);
    std::cout.precision(12);
    for (int al = 0; al < mult; ++al) {
        std::cout << "  alpha = " << al + 1 << ":\n";
        std::cout << "    Q(M)  Q(M')      CGC\n";
        double norm = 0.0;
        for (int j = 0; j < dimS; ++j) {
            for (int k = 0; k < dimSp; ++k) {
                double x = C(j, k, al, hws);
                if (std::fabs(x) > clebsch::EPS) {
                    std::cout << "    " << (j + 1) << "     " << (k + 1)
                              << "     " << x << "\n";
                    norm += x * x;
                }
            }
        }
        std::cout << "    (sum of CGC^2 = " << norm << ")\n";
    }

    // Orthogonality between multiplicity copies (should be ~0).
    if (mult > 1) {
        std::cout << "  overlaps between alpha copies (should be ~0):\n";
        for (int al = 0; al < mult; ++al) {
            for (int be = al + 1; be < mult; ++be) {
                double ov = 0.0;
                for (int j = 0; j < dimS; ++j)
                    for (int k = 0; k < dimSp; ++k)
                        ov += C(j, k, al, hws) * C(j, k, be, hws);
                std::cout << "    <" << al + 1 << "|" << be + 1 << "> = " << ov << "\n";
            }
        }
    }

    // Optional: cross-check against the FULL original computation.
    if (verify) {
        std::cout << "\nVERIFY: running full (all-descendants) computation for cross-check...\n";
        auto f0 = std::chrono::steady_clock::now();
        clebsch::coefficients Cfull(Spp, S, Sp, false);
        auto f1 = std::chrono::steady_clock::now();
        double full_secs = std::chrono::duration<double>(f1 - f0).count();
        double hw_min, hw_med, hw_mean;
        stats(t_total, hw_min, hw_med, hw_mean);
        std::cout << "  full computation (all " << dimSpp << " states) took "
                  << fmt_time(full_secs) << "\n";
        std::cout << "  highest-weight only (min)          "
                  << fmt_time(hw_min) << "\n";
        if (hw_min > 0)
            std::cout << "  speedup (full / hw-only)           "
                      << (full_secs / hw_min) << "x\n";
        double maxdiff = 0.0;
        for (int al = 0; al < mult; ++al)
            for (int j = 0; j < dimS; ++j)
                for (int k = 0; k < dimSp; ++k) {
                    double d = std::fabs(C(j, k, al, hws) - Cfull(j, k, al, hws));
                    if (d > maxdiff) maxdiff = d;
                }
        std::cout << "  max |hw_only - full| over highest-weight coeffs = "
                  << maxdiff
                  << (maxdiff < 1e-9 ? "   OK\n" : "   *** MISMATCH ***\n");
    }

    if (sink == 1234567.89) std::cerr << "";  // keep 'sink' alive
    return 0;
}
