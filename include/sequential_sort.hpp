#ifndef PARALLEL_SORTINGS_INCLUDE_SEQUENTIAL_SORT_HPP_
#define PARALLEL_SORTINGS_INCLUDE_SEQUENTIAL_SORT_HPP_

#include <algorithm>
#include <functional>

#if !defined(_SEQ_STD) && !defined(_SEQ_PDQ)
    #include "../blockquicksort/blocked_double_pivot_check_mosqrt.h++"
#endif

#if defined(_SEQ_PDQ)
    #include <boost/sort/sort.hpp>
#endif

namespace parallel_sortings {
template <typename RandomAccessIterator, typename Compare>
void sequentialSort(RandomAccessIterator first, RandomAccessIterator last,
                    Compare comp) {
#if defined(_SEQ_STD)
    std::sort(first, last, comp);
#elif defined(_SEQ_PDQ)
    boost::sort::pdqsort(first, last, comp);
#else
    blocked_double_pivot_check_mosqrt::sort(first, last, comp);
#endif
}

template <typename RandomAccessIterator>
void sequentialSort(RandomAccessIterator first, RandomAccessIterator last) {
    sequentialSort(first, last, std::less<>());
}
} // namespace parallel_sortings

#endif // PARALLEL_SORTINGS_INCLUDE_SEQUENTIAL_SORT_HPP_
