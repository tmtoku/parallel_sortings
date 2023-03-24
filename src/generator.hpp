#ifndef PARALLEL_SORTINGS_SRC_GENERATOR_HPP_
#define PARALLEL_SORTINGS_SRC_GENERATOR_HPP_

#include <cassert>
#include <limits>
#include <random>

template <typename T, typename Size>
void UniformInt(T A[], Size N, int seed, T min, T max) {
    // random integers in [min, max]
    std::uniform_int_distribution<T> dist(min, max);
    std::mt19937 mt(seed);
    for (Size i = 0; i < N; ++i)
        A[i] = dist(mt);
}

template <typename T, typename Size> void UniformInt(T A[], Size N, int seed) {
    UniformInt(A, N, seed, std::numeric_limits<T>::min(),
               std::numeric_limits<T>::max());
}

template <typename T, typename Size> void UniformFlt(T A[], Size N, int seed) {
    // random numbers in [0, 1)
    std::uniform_real_distribution<T> dist(0, 1);
    std::mt19937 mt(seed);
    for (Size i = 0; i < N; ++i)
        A[i] = dist(mt);
}

template <typename Size> void AlmostSorted(Size A[], Size N, int seed) {
    // a sorted sequence
    for (Size i = 0; i < N; ++i)
        A[i] = i;

    // random integers in [0, N-1]
    std::uniform_int_distribution<Size> dist(0, N - 1);
    std::mt19937 mt(seed);

    // pick sqrt(N) unique numbers and swap them
    Size sqrtN = std::sqrt(N);
    bool *swapped = new bool[N];
    std::fill(swapped, swapped + N, false);

    Size prev = dist(mt);
    swapped[prev] = true;
    for (Size i = 1; i < sqrtN; ++i) {
        Size pos;
        do {
            pos = dist(mt);
        } while (swapped[pos]);
        swapped[pos] = true;
        std::swap(A[prev], A[pos]);
        prev = pos;
    }
    delete[] swapped;

    // check that there are sqrtN out-of-place elements
    Size nOutOfPlace = 0;
    for (Size i = 0; i < N; ++i) {
        if (A[i] != i)
            ++nOutOfPlace;
    }
    assert(nOutOfPlace == sqrtN);
}

#endif // PARALLEL_SORTINGS_SRC_GENERATOR_HPP_
