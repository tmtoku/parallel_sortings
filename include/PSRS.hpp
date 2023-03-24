#ifndef PARALLEL_SORTINGS_INCLUDE_PSRS_HPP_
#define PARALLEL_SORTINGS_INCLUDE_PSRS_HPP_

#include "kway_merge.hpp"
#include "sequential_sort.hpp"

#include <omp.h>

namespace parallel_sortings {
namespace detail {
template <typename RandomAccessIterator1, typename RandomAccessIterator2,
          typename Size, typename Compare>
void PSRS(RandomAccessIterator1 array, Size n, RandomAccessIterator2 result,
          Compare comp) {
    // the number of blocks
    auto nBlock = omp_get_max_threads();
    // block size
    const auto w = (n + nBlock - 1) / nBlock;
    // some blocks may be unnecessary
    nBlock = (n + w - 1) / w;

    // the number of partitions
    const auto nPartition = omp_get_max_threads();

    RandomAccessIterator1 samples[nBlock * nPartition];
    RandomAccessIterator1 pivots[nPartition];
    Size subListStarts[nBlock][nPartition + 1];
    Size partitionStarts[nPartition + 1];

#pragma omp parallel num_threads(nBlock)
    {
        const auto blockID = omp_get_thread_num();
        const auto blockBegin = blockID * w;
        const auto blockEnd = (blockBegin + w > n ? n : blockBegin + w);

        // sort elements in a block
        sequentialSort(array + blockBegin, array + blockEnd, comp);

        // sampling interval
        const auto r = (blockEnd - blockBegin + nPartition - 1) / nPartition;
        // sample nPartition - 1 elements
        for (int iP = 1; iP < nPartition; ++iP) {
            if (blockBegin + iP * r < blockEnd) {
                samples[blockID * (nPartition - 1) + iP - 1] =
                    array + (blockBegin + iP * r);
            } else {
                samples[blockID * (nPartition - 1) + iP - 1] =
                    array + (blockEnd - 1);
            }
        }
#pragma omp barrier

        // sort all samples
#pragma omp single
        sequentialSort(
            samples, samples + nBlock * (nPartition - 1),
            [comp = comp](RandomAccessIterator1 a, RandomAccessIterator1 b) {
                return comp(*a, *b);
            });

        // pick nPartition - 1 pivots
#pragma omp for
        for (int iP = 0; iP < nPartition - 1; ++iP) {
            pivots[iP] = samples[iP * nBlock + (nBlock / 2)];
        }

        // divide locally sorted lists into sublists
        subListStarts[blockID][0] = blockBegin;
        subListStarts[blockID][nPartition] = blockEnd;

        int lteq = blockBegin - 1;
        for (int iP = 0; iP < nPartition - 1; ++iP) {
            const auto pivot = pivots[iP];
            int gt = blockEnd;
            while (gt - lteq > 1) {
                int mid = (lteq + gt) / 2;
                if (comp(*pivot, array[mid]))
                    gt = mid;
                else
                    lteq = mid;
            }
            subListStarts[blockID][iP + 1] = gt;
        }
    }

#pragma omp parallel num_threads(nPartition)
    {
        const auto partitionID = omp_get_thread_num();
        // the start/end pointers of each sublist
        RandomAccessIterator1 subListBegins[nBlock], subListEnds[nBlock];
        // total size of sublists
        Size sumSubListSizes = 0;
        for (int iB = 0; iB < nBlock; ++iB) {
            // the j-th sublist
            Size subListBeginIndex = subListStarts[iB][partitionID];
            Size subListEndIndex = subListStarts[iB][partitionID + 1];
            subListBegins[iB] = array + subListBeginIndex;
            subListEnds[iB] = array + subListEndIndex;
            sumSubListSizes += subListEndIndex - subListBeginIndex;
        }
        partitionStarts[partitionID + 1] = sumSubListSizes;
#pragma omp barrier

        // calculate the start index of each partition
#pragma omp single
        {
            partitionStarts[0] = 0;
            for (int iP = 0; iP < nPartition; ++iP) {
                partitionStarts[iP + 1] += partitionStarts[iP];
            }
        }

        // merge partitions
#if defined(_MERGE_WITH_HEAP)
        KwayMerge_heap(subListBegins, subListEnds, nBlock,
                       result + partitionStarts[partitionID], comp);
#elif defined(_MERGE_BY_SORT)
        RandomAccessIterator2 copyDest = result + partitionStarts[partitionID];
        for (int iB = 0; iB < nBlock; ++iB) {
            std::copy(subListBegins[iB], subListEnds[iB], copyDest);
            copyDest += subListEnds[iB] - subListBegins[iB];
        }
        sequentialSort(result + partitionStarts[partitionID],
                       result + partitionStarts[partitionID + 1], comp);
#else
        KwayMerge(subListBegins, subListEnds, nBlock,
                  result + partitionStarts[partitionID], comp);
#endif
    }
}
} // namespace detail

template <typename RandomAccessIterator1, typename RandomAccessIterator2,
          typename Compare>
void PSRS(RandomAccessIterator1 first, RandomAccessIterator1 last,
          RandomAccessIterator2 result, Compare comp) {
    detail::PSRS(first, last - first, result, comp);
}

template <typename RandomAccessIterator1, typename RandomAccessIterator2>
void PSRS(RandomAccessIterator1 first, RandomAccessIterator1 last,
          RandomAccessIterator2 result) {
    detail::PSRS(first, last - first, result, std::less<>());
}

} // namespace parallel_sortings

#endif // PARALLEL_SORTINGS_INCLUDE_PSRS_HPP_
