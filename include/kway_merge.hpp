#ifndef PARALLEL_SORTINGS_INCLUDE_KWAY_MERGE_HPP_
#define PARALLEL_SORTINGS_INCLUDE_KWAY_MERGE_HPP_

#include <queue>
#include <tuple>
#include <utility>

namespace parallel_sortings {
namespace detail {
template <typename InputIterator, typename Compare>
bool compete(int p1, int p2, const InputIterator begins[],
             const InputIterator ends[], Compare comp) {
    // Player1 beats Player2, if and only if
    // Player2 has no elements or
    // Both players have elements and Player1 has a smaller element
    return (begins[p2] == ends[p2]) ||
           ((begins[p1] != ends[p1]) && comp(*begins[p1], *begins[p2]));
}

template <typename InputIterator, typename Compare>
void buildTree(const InputIterator begins[], const InputIterator ends[], int k,
               int loserTree[], Compare comp) {
    int winnerTree[2 * k];
    // initialize 'winnerTree'
    for (int i = 0; i < k; ++i) {
        winnerTree[k + i] = i;
    }
    // build the loser tree bottom up
    for (int i = k - 1; i > 0; --i) {
        const auto lWinner = winnerTree[2 * i];
        const auto rWinner = winnerTree[2 * i + 1];
        // lWinner beats rWinner
        if (compete(lWinner, rWinner, begins, ends, comp)) {
            winnerTree[i] = lWinner;
            loserTree[i] = rWinner;
        } else {
            winnerTree[i] = rWinner;
            loserTree[i] = lWinner;
        }
    }
    // loserTree[0] has the smallest element
    loserTree[0] = winnerTree[1];
}
} // namespace detail

template <typename InputIterator, typename OutputIterator, typename Compare>
void KwayMerge(InputIterator begins[], const InputIterator ends[], int k,
               OutputIterator result, Compare comp) {
    int loserTree[k];
    detail::buildTree(begins, ends, k, loserTree, comp);

    // loserTree[0] is a current winner
    auto &winner = loserTree[0];
    while (begins[winner] != ends[winner]) {
        // copy the smallest element to 'result'
        *(result++) = *(begins[winner]++);
        // update the tree
        for (int node = (winner + k) >> 1; node > 0; node >>= 1) {
            auto &loser = loserTree[node];
            // if the node beats the current winner, swap them
            if (detail::compete(loser, winner, begins, ends, comp)) {
                std::swap(loser, winner);
            }
        }
    }
}

template <typename InputIterator, typename OutputIterator, typename Compare>
void KwayMerge_heap(InputIterator begins[], const InputIterator ends[], int k,
                    OutputIterator result, Compare comp) {
    using Pair = std::pair<InputIterator, int>;
    const auto greater = [comp = comp](const Pair &a, const Pair &b) {
        // (a > b) = (b < a)
        return comp(*(b.first), *(a.first));
    };
    // the first element will be the smallest
    std::priority_queue<Pair, std::vector<Pair>, decltype(greater)> minHeap(
        greater);
    for (int i = 0; i < k; ++i) {
        if (begins[i] != ends[i])
            minHeap.emplace(begins[i], i);
    }

    while (!minHeap.empty()) {
        InputIterator dat;
        int i;
        std::tie(dat, i) = minHeap.top();
        minHeap.pop();
        // copy the smallest element to 'result'
        *(result++) = *(begins[i]++);
        if (begins[i] != ends[i]) {
            minHeap.emplace(begins[i], i);
        }
    }
}

} // namespace parallel_sortings

#endif // PARALLEL_SORTINGS_INCLUDE_KWAY_MERGE_HPP_
