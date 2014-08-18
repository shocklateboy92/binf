#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include<limits>

#include <boost/multi_array.hpp>

#include <assert.h>

typedef std::vector<double> vec;
typedef boost::multi_array<double, 2> mat;

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cerr << "Incorrect args!" << std::endl;
    }

    std::size_t num_states, num_symbols, num_observations;
    std::ifstream is(argv[1]);
    assert(is.is_open());

    is >> num_states;
    is >> num_symbols;
    is >> num_observations;

    mat transitions(boost::extents[num_states][num_states]);
    mat emissions(boost::extents[num_states][num_states]);
    vec observations(num_states);

    for (std::size_t row = 0; row < num_states; row++) {
        for (std::size_t col = 0; col < num_states; col++) {
            assert(is.good());
            is >> transitions[row][col];
        }
    }

    for (std::size_t row = 0; row < num_states; row++) {
        for (std::size_t col = 0; col < num_symbols; col++) {
            assert(is.good());
            is >> emissions[row][col];
        }
    }

    for (std::size_t row = 0; row < num_observations; row++) {
        assert(is.good());
        is >> observations[row];
    }

    mat V(boost::extents[num_observations][num_states]);

    for (std::size_t s = 0; s < num_states; s++) {
        V[0][s] = transitions[0][s] * emissions[s][observations[0]];
    }

    for (std::size_t t = 1; t < num_observations; t++) {
        for (std::size_t s1 = 1; s1 < num_states; s1++) {
            double max = std::numeric_limits<double>::min();
            std::size_t state_max = -1;
            for (std::size_t s2 = 0; s2 < num_states; s2++) {
                double prob = V[t-1][s2] *
                        transitions[s2][s1] *
                        emissions[s1][observations[t]];
                if (prob > max) {
                    max = prob;
                    state_max = s2;
                }
            }
            V[t][s1] = max;
        }
    }

    return 0;
}

