#include "generator.hpp"
#include "sorter.hpp"

#include <iostream>

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "usage: %s <N> <count> <algoType>\n", argv[0]);
        std::exit(1);
    }
    using Size = uint32_t;

    const Size N = atoi(argv[1]);
    const int count = atoi(argv[2]);
    const int algoType = atoi(argv[3]);

    constexpr int NUM_INPUTS = 6;
    double times[NUM_INPUTS][count];
    for (int i = 0; i < NUM_INPUTS; ++i) {
        std::fill(times[i], times[i] + count, 0);
    }

    {
        // UniformInt
        auto A = new uint32_t[N];
#pragma omp parallel for
        for (Size i = 0; i < N; ++i) {
            A[i] = 0;
        }
        for (int seed = 0; seed < count; ++seed) {
            UniformInt(A, N, seed);
            times[0][seed] = sort_time(A, N, algoType);
        }
        delete[] A;
    }
    {
        // UniformFloat
        auto A = new float[N];
#pragma omp parallel for
        for (Size i = 0; i < N; ++i) {
            A[i] = 0;
        }
        for (int seed = 0; seed < count; ++seed) {
            UniformFlt(A, N, seed);
            times[1][seed] = sort_time(A, N, algoType);
        }
        delete[] A;
    }
    {
        // AlmostSorted
        auto A = new Size[N];
#pragma omp parallel for
        for (Size i = 0; i < N; ++i) {
            A[i] = 0;
        }
        for (int seed = 0; seed < count; ++seed) {
            AlmostSorted(A, N, seed);
            times[2][seed] = sort_time(A, N, algoType);
        }
        delete[] A;
    }
    {
        // Duplicate3
        auto A = new uint32_t[N];
#pragma omp parallel for
        for (Size i = 0; i < N; ++i) {
            A[i] = 0;
        }
        for (int seed = 0; seed < count; ++seed) {
            UniformInt(A, N, seed, 0u, 2u);
            times[3][seed] = sort_time(A, N, algoType);
        }
        delete[] A;
    }
    {
        // Pair
        auto keys = new uint64_t[N];
        auto A = new Pair<uint64_t, Size>[N];
#pragma omp parallel for
        for (Size i = 0; i < N; ++i) {
            keys[i] = A[i].key = 0;
        }
        for (int seed = 0; seed < count; ++seed) {
            UniformInt(keys, N, seed);
            for (Size i = 0; i < N; ++i) {
                A[i].key = keys[i];
            }
            times[4][seed] = sort_time(A, N, algoType);
        }
        delete[] keys;
        delete[] A;
    }
    {
        // Particle
        auto keys = new uint64_t[N];
        auto A = new Particle<uint64_t>[N];
#pragma omp parallel for
        for (Size i = 0; i < N; ++i) {
            keys[i] = A[i].key = 0;
        }
        for (int seed = 0; seed < count; ++seed) {
            UniformInt(keys, N, seed);
            for (Size i = 0; i < N; ++i) {
                A[i].key = keys[i];
            }
            times[5][seed] = sort_time(A, N, algoType);
        }
        delete[] keys;
        delete[] A;
    }

    /* Output the results of all 'count' times to stderr */
    std::cerr << getAlgoString(algoType) << " " << omp_get_max_threads()
              << "\n";
    fprintf(stderr,
            "## uniform_int uniform_flt almost_sort duplicate_3 Pair_(%3lu) "
            "Particle_(%3lu)\n",
            sizeof(Pair<uint64_t, Size>), sizeof(Particle<uint64_t>));

    for (int seed = 0; seed < count; ++seed) {
        fprintf(stderr, "%2u ", seed);
        for (int i = 0; i < NUM_INPUTS; ++i) {
            double t = times[i][seed];
            if (i < 4)
                fprintf(stderr, "%11.3f ", t);
            else if (i == 4)
                fprintf(stderr, "%10.3f ", t);
            else
                fprintf(stderr, "%14.3f ", t);
        }
        fprintf(stderr, "\n");
    }

    /* Output the average of the first 'count' times to stdout */
    printf("## uniform_int uniform_flt almost_sort duplicate_3 Pair_(%3lu) "
           "Particle_(%3lu)\n",
           sizeof(Pair<uint64_t, Size>), sizeof(Particle<uint64_t>));
    printf("%2d ", omp_get_max_threads());

    for (int i = 0; i < NUM_INPUTS; ++i) {
        double total_time = 0.0;
        for (int seed = 0; seed < count; ++seed) {
            double t = times[i][seed];
            total_time += t;
        }

        if (i < 4)
            printf("%11.3f ", total_time / count);
        else if (i == 4)
            printf("%10.3f ", total_time / count);
        else
            printf("%14.3f ", total_time / count);
    }
    printf("\n");
}
