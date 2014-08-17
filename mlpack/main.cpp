#include <iostream>
#include <mlpack/core.hpp>
#include <mlpack/methods/hmm/hmm.hpp>
#include <mlpack/core/dists/discrete_distribution.hpp>
#include <armadillo>
#include <armadillo_bits/arma_cmath.hpp>

using namespace std;
using namespace mlpack;

int main()
{
    std::size_t num_states = 2;
    arma::vec initial = {0.5, 0.5};
    arma::mat transitions = {
        0.1d, 0.9d,
        0.4d, 0.6d
    };
    transitions.reshape(num_states, num_states);
    std::vector<distribution::DiscreteDistribution> emissions = {
        distribution::DiscreteDistribution(num_states),
        distribution::DiscreteDistribution(num_states)
    };

    hmm::HMM<distribution::DiscreteDistribution> model(transitions, emissions);
    arma::mat observations = {1, 0};
    arma::Col<std::size_t> results;
    double score = model.Predict(observations, results);
    std::cout << results;
    std::cout << score << std::endl;
    return 0;
}

