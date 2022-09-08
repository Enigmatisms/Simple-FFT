#include "../include/simple_fft/fft.h"
#include <iostream>
#include <complex>
#include <fstream>
#include <algorithm>
#include <chrono>

constexpr float PI = 3.141592653589793;

template <typename T>
std::complex<T> freqButterworth4(int k, T cutoff_f, T base_f) {
    const std::complex<T> s(0., base_f * static_cast<T>(k) / cutoff_f);
    const std::complex<T> s2 = s * s;   
    return 1. / (s2 + 0.765367 * s + 1.) / (s2 + 1.847759 * s + 1.);
}

template <typename T>
void modulateThenTransform(const std::vector<T>& tx, const std::vector<T>& rx, std::vector<T>& spect, std::vector<T>& unfilter, T cutoff_f = 0.001) {
    using complex_t = std::complex<T>;
    const char* error_description = nullptr;
    std::vector<T> modulated(tx.size());
    std::vector<complex_t> outputs(tx.size());
    std::transform(tx.begin(), tx.end(), rx.begin(), modulated.begin(), [](T v1, T v2){return v1 * v2;});
    simple_fft::FFT(modulated, outputs, modulated.size(), error_description);
    auto get_amp_func = [](const complex_t& c) {return sqrt(pow(c.real(), 2) + pow(c.imag(), 2));};

    if (cutoff_f > 1e-9) {
        const T base_f = 2 * PI / static_cast<T>(tx.size());
        std::vector<complex_t> butterworth_lp(modulated.size());
        std::generate(butterworth_lp.begin(), butterworth_lp.end(), 
            [k = 0, cutoff_f, base_f]() mutable {
                return freqButterworth4(k, cutoff_f, base_f);
            }
        );
        std::transform(outputs.begin(), outputs.end(), unfilter.begin(), get_amp_func);
        std::transform(outputs.begin(), outputs.end(), butterworth_lp.begin(), outputs.begin(), [](const auto& v1, const auto& v2){return v1 * v2;});
    }
    std::transform(outputs.begin(), outputs.end(), spect.begin(), get_amp_func);
}

int main() {
    std::ofstream of2("../unfilter.txt");
    std::ofstream of1("../filter.txt");
    std::vector<double> vs1(1000), vs2(1000);
    for (int i = 0; i < 24; i++) {
        vs1.push_back(0.);
        vs2.push_back(0.);
    }
    std::vector<double> unfilter_out(1024), filter_out(1024);
    std::generate(vs1.begin(), vs1.end(), [n = 0]() mutable {return n++;});
    std::generate(vs2.begin(), vs2.end(), [n = 0]() mutable {return n++;});
    std::transform(vs1.begin(), vs1.end(), vs1.begin(), [](double val) {return sin(0.1 * val);});
    std::transform(vs2.begin(), vs2.end(), vs2.begin(), [](double val) {return sin(0.11 * val);});
    const char * error_description = 0;
    
    // auto start_t = std::chrono::system_clock::now();
    // simple_fft::FFT(vals, outputs, vals.size(), error_description);
    // auto end_t = std::chrono::system_clock::now();
    // printf("Time consumed: %f ms\n", static_cast<double>((end_t - start_t).count()) / 1e6);

    // std::vector<double> amp(1024);
    // std::transform(outputs.begin(), outputs.end(), amp.begin(), [](std::complex<double>& cv) {return sqrt(pow(cv.real(), 2) + pow(cv.imag(), 2));});

    modulateThenTransform(vs1, vs2, filter_out, unfilter_out, 0.00001);
    for (double val: filter_out) {
        of1 << val << std::endl;
    }
    for (double val: unfilter_out) {
        of2 << val << std::endl;
    }
    of1.close();
    of2.close();
    return 0;
}