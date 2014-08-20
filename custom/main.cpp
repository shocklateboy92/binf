#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <unordered_map>
#include<limits>

#include <boost/multi_array.hpp>

#include <assert.h>

typedef std::vector<double> vec;
typedef boost::multi_array<double, 2> mat;
typedef std::vector<std::unordered_map<std::size_t, double>> s_mat;
typedef std::pair<std::size_t, double> s_val;

static const double SCALE = 4.0d;

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cerr << "Incorrect args!" << std::endl;
    }

    std::size_t num_states, num_symbols, num_observations;
    std::ifstream is(argv[1]);
    assert(is.is_open());

    // read model info
    is >> num_states;
    is >> num_symbols;
    is >> num_observations;

    s_mat transitions(num_states);
    mat emissions(boost::extents[num_states][num_states]);
    vec observations(num_states);

    // read transition probabilities between states
    for (std::size_t row = 0; row < num_states; row++) {
        for (std::size_t col = 0; col < num_states; col++) {
            assert(is.good());
            double val;
            is >> val;
            if (val) {
                transitions[row][col] = val;
            }
        }
    }

    // read emission probabilities of each state
    for (std::size_t row = 0; row < num_states; row++) {
        for (std::size_t col = 0; col < num_symbols; col++) {
            assert(is.good());
            is >> emissions[row][col];
        }
    }

    // read the observations
    for (std::size_t row = 0; row < num_observations; row++) {
        assert(is.good());
        is >> observations[row];
    }

    s_mat V(num_observations);

    // set state 0 as the initial state
    V[0][0] = 1.0d;

    // do viterbi iterative viterbi algorithm
    for (std::size_t t = 1; t < num_observations; t++) {
        for (s_val v : V[t-1]) {
            for (s_val transition : transitions[v.first]) {
                double v_next = v.second *
                        transition.second *
                        SCALE *
                        emissions[transition.first][observations[t-1]];
                if (V[t][transition.first] <= v_next) {
                    V[t][transition.first] = v_next;
                }
            }
        }
    }

    // print the most likely path
    for (std::unordered_map<std::size_t, double> v : V) {
        auto max = *std::max_element(v.begin(), v.end(),
                                    [](s_val a, s_val b){
                                        return a.second < b.second;
                                    });
        std::cout << max.first << ":" << max.second << std::endl;
    }

    return 0;
}

