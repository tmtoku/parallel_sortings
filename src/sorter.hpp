#ifndef PARALLEL_SORTINGS_SRC_SORTER_HPP_
#define PARALLEL_SORTINGS_SRC_SORTER_HPP_

#include "../include/PSES.hpp"
#include "../include/PSRS.hpp"
#include "../include/sequential_sort.hpp"
#include "utils.hpp"

#include <parallel/algorithm>

template <typename T, typename Size>
double sort_time(T A[], Size N, int algoType) {
    switch (algoType) {
    case AlgoType::SEQUENTIAL: {
        const auto t = timer<void (*)(T *, T *)>(
            parallel_sortings::sequentialSort)(A, A + N);
        assert(isSorted(A, A + N));
        return t;
    }
    case AlgoType::GNU_PARALLEL: {
        const auto t =
            timer<void (*)(T *, T *)>(__gnu_parallel::sort)(A, A + N);
        assert(isSorted(A, A + N));
        return t;
    }
    case AlgoType::PSRS: {
        T *B = new T[N];
#pragma omp parallel for
        for (Size i = 0; i < N; ++i) {
            B[i] = A[i];
        }
        const auto t = timer<void (*)(T *, T *, T *)>(parallel_sortings::PSRS)(
            A, A + N, B);
        assert(isSorted(B, B + N));
        delete[] B;
        return t;
    }
    case AlgoType::PSES: {
        T *B = new T[N];
#pragma omp parallel for
        for (Size i = 0; i < N; ++i) {
            B[i] = A[i];
        }
        const auto t = timer<void (*)(T *, T *, T *)>(parallel_sortings::PSES)(
            A, A + N, B);
        assert(isSorted(B, B + N));
        delete[] B;
        return t;
    }
    default: {
        fprintf(stderr, "algo type %d is not supported.\n", algoType);
        exit(1);
    }
    }
}

#endif // PARALLEL_SORTINGS_SRC_SORTER_HPP_
