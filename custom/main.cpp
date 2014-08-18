#include <iostream>
#include <fstream>
#include <vector>

#include <assert.h>

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

    std::vector<std::vector<double>> transitions;
    transitions.reserve(num_states);
    std::vector<std::vector<double>> emissions;
    emissions.reserve(num_states);
    std::vector<std::size_t> observations(num_states);

    for (std::size_t row = 0; row < num_states; row++) {
        transitions[row] = std::vector<double>(num_states);
        for (std::size_t col = 0; col < num_states; col++) {
            assert(is.good());
            is >> transitions[row][col];
        }
    }

    for (std::size_t row = 0; row < num_states; row++) {
        emissions[row] = std::vector<double>(num_symbols);
        for (std::size_t col = 0; col < num_symbols; col++) {
            assert(is.good());
            is >> emissions[row][col];
        }
    }

    for (std::size_t row = 0; row < num_observations; row++) {
        assert(is.good());
        is >> observations[row];
    }

    return 0;
}

