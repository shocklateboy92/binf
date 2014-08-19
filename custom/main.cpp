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

struct state {
    std::size_t id;
    std::unordered_map<std::size_t, double> transitions;
//    std::vector<std::pair<std::size_t, double>> transitions;
};

struct partV {
    double prob;
    state *s;

    partV() : prob(0), s(nullptr) {}
    partV(double _prob, state *_s) : prob(_prob), s(_s) {}
};

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
    std::vector<state> states(num_states);

    for (std::size_t row = 0; row < num_states; row++) {
        for (std::size_t col = 0; col < num_states; col++) {
            assert(is.good());
//            is >> transitions[row][col];
            double val;
            is >> val;
            if (val) {
                states[row].transitions[col] = val;
            }
            states[row].id = row;
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

//    mat V(boost::extents[num_observations][num_states]);
    std::vector<std::unordered_map<std::size_t, partV>> V(num_observations);

    for (state &s : states) {
        partV v = { 1.0d / (double) num_states, &s};
        V[0][s.id] = v;
//        std::cout << v.prob << " ";
//        for (std::size_t i = 0; i < num_states; i++) {
//            V[0][i].prob = 1.0d / (double) num_states;
//        }
    }

    for (std::size_t t = 1; t < num_observations; t++) {
        for (auto v : V[t-1]) {
//            double max = std::numeric_limits<double>::min();
            for (std::pair<std::size_t, double> transition : states[v.first].transitions) {
                std::size_t sid = transition.first;
                double p2 = v.second.prob *
                        transition.second *
                        emissions[sid][observations[t]];
                if (V[t][transition.first].prob < p2) {
                    V[t][transition.first].prob = p2;
                    V[t][transition.first].s = &states[transition.first];
                }
            }
        }
    }

//    for (auto v : V) {
//        for (auto p : v) {
////            std::cout << p.first << ":" << p.second.prob << " ";
////            std::cout << p.count();
//        }
//        std::cout << v.size();
//        std::endl(std::cout);
//    }

    return 0;
}

