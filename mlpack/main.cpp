#include <iostream>
#include <fstream>
#include <mlpack/core.hpp>
#include <mlpack/methods/hmm/hmm.hpp>
#include <mlpack/core/dists/discrete_distribution.hpp>
#include <armadillo>
#include <armadillo_bits/arma_cmath.hpp>

//using namespace std;
using namespace mlpack;
using namespace arma;

int main(int argc, char *argv[])
{
    if (argc < 2) {
        return -1;
    }
    std::ifstream is(argv[1]);

    std::size_t num_states = 2;
    std::size_t num_symbols = 2;
    std::size_t num_observations = 2;
    is >> num_states;
    is >> num_symbols;
    is >> num_observations;

    arma::mat *transitions = new arma::mat(num_states, num_states);
    std::vector<distribution::DiscreteDistribution> emissions;
    arma::rowvec observations(num_observations);

    for (std::size_t row = 0; row < num_states; row++) {
        for (std::size_t col = 0; col < num_states; col++) {
            assert(is.good());
            is >> (*transitions)(row, col);
        }
    }

    for (std::size_t row = 0; row < num_states; row++) {
        assert(is.good());
        arma::vec probs(num_symbols);
        for (std::size_t col = 0; col < num_symbols; col++) {
            is >> probs(col);
        }
        emissions.push_back(distribution::DiscreteDistribution(probs));
    }

    for (std::size_t col = 0; col < num_observations; col++) {
        assert(is.good());
        is >> observations(col);
    }

    hmm::HMM<distribution::DiscreteDistribution> model(*transitions, emissions);
    arma::Col<std::size_t> results;
    double score = model.Predict(observations, results);
    std::cout << results;
    std::cout << score << std::endl;
    return 0;
}

