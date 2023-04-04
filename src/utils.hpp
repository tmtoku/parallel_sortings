#ifndef PARALLEL_SORTINGS_SRC_UTILS_HPP_
#define PARALLEL_SORTINGS_SRC_UTILS_HPP_

#include <chrono>
#include <functional>
#include <string>

// measure the elapsed time of 'func'
template <typename Func> auto timer(Func &&func) {
    return [func = std::forward<Func>(func)](auto &&...args) {
        auto start_time = std::chrono::high_resolution_clock::now();
        func(std::forward<decltype(args)>(args)...);
        auto end_time = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::milli>(end_time - start_time)
            .count();
    };
}

// ensure that the input array is sorted
template <typename RandomAccessIterator>
bool isSorted(RandomAccessIterator begin, RandomAccessIterator end) {
    for (auto it = begin; it + 1 < end; ++it)
        if (*it > *(it + 1))
            return false;
    return true;
}

inline std::string getSequentialString() {
#if defined(_SEQ_STD)
    return "std";
#elif defined(_SEQ_PDQ)
    return "boost pdq";
#else
    return "blockq";
#endif
}

enum AlgoType {
    SEQUENTIAL,
    GNU_PARALLEL,
    PSRS,
    PSES,
};

inline std::string getAlgoString(int algoType) {
    switch (algoType) {
    case AlgoType::SEQUENTIAL: {
        return getSequentialString();
    }
    case AlgoType::GNU_PARALLEL: {
        return "gnu parallel";
    }
    case AlgoType::PSRS: {
        return "PSRS " + getSequentialString();
    }
    case AlgoType::PSES: {
        return "PSES " + getSequentialString();
    }
    default: {
        fprintf(stderr, "algo type %d is not supported.\n", algoType);
        exit(1);
    }
    }
}

template <typename T, typename Size> struct Pair {
    T key;
    Size index;

    bool operator<(const Pair &b) const { return key < b.key; }
    bool operator>(const Pair &b) const { return key > b.key; }
    bool operator>=(const Pair &b) const { return key >= b.key; }
};

template <typename T> struct Particle {
    T key;
    double position[3];
    double velocity[3];
    double acceleration[3];
    double potential;
    double mass;

    bool operator<(const Particle &b) const { return key < b.key; }
    bool operator>(const Particle &b) const { return key > b.key; }
    bool operator>=(const Particle &b) const { return key >= b.key; }
};

#endif // PARALLEL_SORTINGS_SRC_UTILS_HPP_
