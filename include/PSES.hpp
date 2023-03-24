#ifndef PARALLEL_SORTINGS_INCLUDE_PSES_HPP_
#define PARALLEL_SORTINGS_INCLUDE_PSES_HPP_

#include "kway_merge.hpp"
#include "sequential_sort.hpp"

#include <omp.h>

namespace parallel_sortings {
namespace detail {
template <typename RandomAccessIterator, typename Size, typename Compare>
Size countLtEq(RandomAccessIterator sortedArray, Size n,
               const RandomAccessIterator target, Compare comp) {
    if (!comp(*target, sortedArray[0])) {
        // sortedArray[0] <= x
        Size lteq = 0, gt = n;
        while (gt - lteq > 1) {
            const auto mid = lteq + (gt - lteq) / 2;
            if (comp(*target, sortedArray[mid]))
                gt = mid;
            else
                lteq = mid;
        }
        return gt;
    } else {
        return 0;
    }
}

template <typename RandomAccessIterator, typename Size, typename Compare>
Size countLt(RandomAccessIterator sortedArray, Size n,
             const RandomAccessIterator target, Compare comp) {
    if (comp(sortedArray[0], *target)) {
        Size lt = 0, gteq = n;
        while (gteq - lt > 1) {
            const auto mid = lt + (gteq - lt) / 2;
            if (comp(sortedArray[mid], *target))
                lt = mid;
            else
                gteq = mid;
        }
        return gteq;
    } else {
        return 0;
    }
}

template <typename Size> struct Block {
    Size begin;
    Size end;
};

template <typename RandomAccessIterator, typename Size, typename Compare>
Size sumLtEq(RandomAccessIterator array, Block<Size> blocks[], int nBlock,
             const RandomAccessIterator target, Compare comp) {
    Size nLessThanEq = 0;
    for (int i = 0; i < nBlock; ++i) {
        nLessThanEq += countLtEq(array + blocks[i].begin,
                                 blocks[i].end - blocks[i].begin, target, comp);
    }
    return nLessThanEq;
}

template <typename RandomAccessIterator, typename Size, typename Compare>
Size sumLt(RandomAccessIterator array, Block<Size> blocks[], int nBlock,
           const RandomAccessIterator target, Compare comp) {
    Size nLessThan = 0;
    for (int i = 0; i < nBlock; ++i) {
        nLessThan += countLt(array + blocks[i].begin,
                             blocks[i].end - blocks[i].begin, target, comp);
    }
    return nLessThan;
}

template <typename Iterator, typename Compare>
Iterator medianOfThree(const Iterator a, const Iterator b, const Iterator c,
                       Compare comp) {
    if (comp(*b, *a) ^ comp(*c, *a))
        return a;
    else if (comp(*b, *a) ^ comp(*b, *c))
        return b;
    else
        return c;
}

template <typename Iterator, typename Compare>
Iterator calcGlobalMedian(const Iterator localMedians[], int nLocalMedian,
                          Compare comp) {
    if (nLocalMedian <= 3) {
        if (nLocalMedian <= 2) {
            return localMedians[0];
        } else {
            return medianOfThree(localMedians[0], localMedians[1],
                                 localMedians[2], comp);
        }
    } else {
        const auto n3 = nLocalMedian / 3;

        const auto a = calcGlobalMedian(localMedians, n3, comp);
        const auto b = calcGlobalMedian(localMedians + n3, n3, comp);
        const auto c = calcGlobalMedian(localMedians + 2 * n3,
                                        nLocalMedian - 2 * n3, comp);
        return medianOfThree(a, b, c, comp);
    }
}

template <typename Iterator, typename Size> struct Splitter {
    Iterator pivot;
    Size numSameKeyIncluded;
};

template <typename RandomAccessIterator, typename Size, typename Compare>
void findSplitter(RandomAccessIterator array, Block<Size> blocks[], int nBlock,
                  Size count, const RandomAccessIterator globalMin,
                  const RandomAccessIterator globalMax,
                  Splitter<RandomAccessIterator, Size> &splitter,
                  Compare comp) {
    // the number of elements <= 'globalMin'
    if (count <= sumLtEq(array, blocks, nBlock, globalMin, comp)) {
        // the number of elements < 'globalMin' is always 0,
        // so 'globalMin' can be a splitter
        splitter = {globalMin, count};
        return;
    }

    // the number of elements < 'globalMax'
    const auto nLtMax = sumLt(array, blocks, nBlock, globalMax, comp);
    if (nLtMax <= count) {
        // the number of elements <= 'globalMax' is always n (>= count)
        // so 'globalMax' can be a splitter
        splitter = {globalMax, count - nLtMax};
        return;
    }

    RandomAccessIterator localMedians[nBlock];
    RandomAccessIterator globalMedian;

    auto lb = globalMin, ub = globalMax;
    while (true) {
        /* pivot to be searched is in (lb, ub) */
        // 1. gather local medians
        int nLocalMedian = 0;
        for (int i = 0; i < nBlock; ++i) {
            const auto blockBegin = blocks[i].begin;
            const auto blockEnd = blocks[i].end;

            // array[x] is the first element that is greater than 'lb'
            const auto x =
                countLtEq(array + blockBegin, blockEnd - blockBegin, lb, comp);
            // array[y-1] is the last element that is less than 'ub'
            const auto y =
                countLt(array + blockBegin, blockEnd - blockBegin, ub, comp);
            assert(x <= y);

            // The middle element of [x,y) is a local median
            if (x < y) {
                const auto mid = blockBegin + (x + y) / 2;
                assert(mid < blockEnd);
                localMedians[nLocalMedian++] = array + mid;
            }
        }
        // 2. determine 'globalMedian'
        globalMedian = calcGlobalMedian(localMedians, nLocalMedian, comp);

        const auto nLessThan = sumLt(array, blocks, nBlock, globalMedian, comp);
        // 'globalMedian' is too large to be a splitter
        if (nLessThan > count) {
            ub = globalMedian;
            continue;
        }

        const auto nLessThanEq =
            sumLtEq(array, blocks, nBlock, globalMedian, comp);
        // 'globalMedian' is too small to be a splitter
        if (nLessThanEq < count) {
            lb = globalMedian;
            continue;
        }

        // 'globalMedian' is eligible to be a splitter
        splitter = {globalMedian, count - nLessThan};
        return;
    }
}

template <typename RandomAccessIterator1, typename RandomAccessIterator2,
          typename Size, typename Compare>
void PSES(RandomAccessIterator1 array, Size n, RandomAccessIterator2 result,
          Compare comp) {
    // the number of blocks
    auto nBlock = omp_get_max_threads();
    // block size
    const auto w = (n + nBlock - 1) / nBlock;
    // some blocks may be unnecessary
    nBlock = (n + w - 1) / w;

    // the number of partitions
    const auto nPartition = omp_get_max_threads();

    Block<Size> blocks[nBlock];
    Splitter<RandomAccessIterator1, Size> splitters[nPartition];
    Size subListStarts[nBlock][nPartition + 1];
    Size partitionStarts[nPartition + 1];

    auto globalMin = array;
    auto globalMax = globalMin;
    Size sameKeySum[nPartition][nBlock + 1];

#pragma omp parallel num_threads(nBlock)
    {
        const auto blockID = omp_get_thread_num();
        const auto blockBegin = blockID * w;
        const auto blockEnd = (blockBegin + w > n ? n : blockBegin + w);
        blocks[blockID] = {blockBegin, blockEnd};
        /* locally sort */
        sequentialSort(array + blockBegin, array + blockEnd, comp);
#pragma omp barrier

        // select multi-pivots
        // calculate global min / max
        const auto localMin = array + blockBegin;
        const auto localMax = array + (blockEnd - 1);
#pragma omp critical
        {
            if (comp(*localMin, *globalMin))
                globalMin = localMin;
            if (comp(*globalMax, *localMax))
                globalMax = localMax;
        }
#pragma omp barrier

#pragma omp for
        for (int i = 0; i < nPartition - 1; ++i) {
            findSplitter(array, blocks, nBlock, blocks[i].end, globalMin,
                         globalMax, splitters[i], comp);
        }

#pragma omp for
        for (int i = 0; i < nPartition - 1; ++i) {
            sameKeySum[i][0] = 0;
            for (int j = 0; j < nBlock; ++j) {
                sameKeySum[i][j + 1] =
                    sameKeySum[i][j] +
                    (countLtEq(array + blocks[j].begin,
                               blocks[j].end - blocks[j].begin,
                               splitters[i].pivot, comp) -
                     countLt(array + blocks[j].begin,
                             blocks[j].end - blocks[j].begin,
                             splitters[i].pivot, comp));
            }
        }

        /* divide each locally sorted list into sublists */
        subListStarts[blockID][0] = blockBegin;
        subListStarts[blockID][nPartition] = blockEnd;

        int lteq = blockBegin - 1;
        for (int iP = 0; iP < nPartition - 1; ++iP) {
            const auto pivot = splitters[iP].pivot;
            int gt = blockEnd;
            while (gt - lteq > 1) {
                int mid = (lteq + gt) / 2;
                if (comp(*pivot, array[mid]))
                    gt = mid;
                else
                    lteq = mid;
            }
            Size numSameKeyInBlock =
                sameKeySum[iP][blockID + 1] - sameKeySum[iP][blockID];
            /* all the same keys in the block are temporarily excluded */
            subListStarts[blockID][iP + 1] = gt - numSameKeyInBlock;

            /* still need to include the same keys */
            if (sameKeySum[iP][blockID] < splitters[iP].numSameKeyIncluded) {
                Size numSameKeyRequired =
                    splitters[iP].numSameKeyIncluded - sameKeySum[iP][blockID];
                if (numSameKeyInBlock < numSameKeyRequired) {
                    /* all the same keys in the block are included */
                    subListStarts[blockID][iP + 1] += numSameKeyInBlock;
                } else {
                    /* Some of the same keys in the block are included */
                    subListStarts[blockID][iP + 1] += numSameKeyRequired;
                }
            }
        }
    }

#pragma omp parallel num_threads(nPartition)
    {
        const auto partitionID = omp_get_thread_num();
        /* the start/end pointer of each sublist */
        RandomAccessIterator1 subListBegins[nBlock], subListEnds[nBlock];
        /* sum of the sublist sizes */
        Size sumSubListSizes = 0;
        for (int iB = 0; iB < nBlock; ++iB) {
            /* j-th sublist */
            Size subListBeginIndex = subListStarts[iB][partitionID];
            Size subListEndIndex = subListStarts[iB][partitionID + 1];
            subListBegins[iB] = array + subListBeginIndex;
            subListEnds[iB] = array + subListEndIndex;
            sumSubListSizes += subListEndIndex - subListBeginIndex;
        }
        /* the size of threadID-th partition */
        partitionStarts[partitionID + 1] = sumSubListSizes;
#pragma omp barrier

        /* the start index of each partition */
#pragma omp single
        {
            partitionStarts[0] = 0;
            for (int iP = 0; iP < nPartition; ++iP) {
                partitionStarts[iP + 1] += partitionStarts[iP];
            }
        }

        /* merge threadID-th partition */
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
void PSES(RandomAccessIterator1 first, RandomAccessIterator1 last,
          RandomAccessIterator2 result, Compare comp) {
    detail::PSES(first, last - first, result, comp);
}

template <typename RandomAccessIterator1, typename RandomAccessIterator2>
void PSES(RandomAccessIterator1 first, RandomAccessIterator1 last,
          RandomAccessIterator2 result) {
    PSES(first, last, result, std::less<>());
}

} // namespace parallel_sortings

#endif // PARALLEL_SORTINGS_INCLUDE_PSES_HPP_
