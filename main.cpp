#include <iostream>
#include <random>
#include <vector>
#include <utility>
#include <fstream>
#include <string>
#include <tuple>
#include <array>
#include <algorithm>
#include <memory>
#include <sstream>


typedef std::pair<int, int> index;
typedef std::pair<double, double> data;
typedef std::vector<std::vector<int>> spins;


const int N = 10;
const int left_border = 0;
const int right_border = 9;
const int number_of_experiments = 100;


std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()


void lattice_filling (spins & lattice);

int energy_of_system (spins & lattice);

double temperature_definition (const double & E_mean, const int & spins_count);

void random_spin_flip (spins & lattice, const int & left, const int & right);

bool new_configuration (const double & dE, const double & T);

bool saturation (std::vector<int> & Energies, const int & error);

std::vector<int> established_regime (std::vector<int> & Energies, const int & error);

void data_file_creation (const std::string & name, std::vector<data> & exp_data);

void vector_clear_2d (spins & vector);

void mem_allocation_2d (spins & vector, const int & dim);

void plot (const std::string & name, const int & left, const int & right,
           const std::string & xlabel, const std::string & ylabel, const std::string & title);


int main () {
    spins lattice (N);
    mem_allocation_2d (lattice, N);
    lattice_filling(lattice);
    int E, E0 = energy_of_system(lattice);
    double T = temperature_definition(double(E0) / double(N), N);
    int i = 0;
    int j = 0;
    std::vector<int> Energies;
    Energies.emplace_back(E0);
    do {
        do {
            random_spin_flip(lattice, left_border, right_border);
            E = energy_of_system(lattice);
        } while (new_configuration(E0 - E, T));
        ++i;
        Energies.emplace_back(E);
    } while (!saturation(Energies, 5));
    return 0;
}


std::vector<int> established_regime (std::vector<int> & Energies, const int & error) {
    std::vector<int> stable;
    for (int i = 1; i < Energies.size(); ++i)
        if (std::abs(Energies[i] - Energies[i-1]) < error)
            stable.emplace_back(Energies[i]);
    return stable;
}


bool saturation (std::vector<int> & Energies, const int & error) {
    for (int i = 1; i < Energies.size(); ++i)
        if (std::abs(Energies[i] - Energies[i-1]) < error)
            return true;
    return false;
}


void lattice_filling (spins & lattice) {
    for (auto & i : lattice)
        for (int & j : i)
            j = 1;
}


int energy_of_system (spins & lattice) {
    int sum = 0;
    for (int i = 0; i < lattice.size(); ++i)
        for (int j = 0; j < lattice[i].size(); ++j) {
            int S_0 = lattice[i][j];
            int S_1 = (i == N-1) ? lattice[0][j] : lattice[i+1][j];
            int S_2 = (j == N-1) ? lattice[i][0] : lattice[i][j+1];
            sum += S_0*S_1 + S_0*S_2;
        }
    return sum;
}


double temperature_definition (const double & E_mean, const int & spins_count) {
    return -1.0 / std::atanh(E_mean / double(spins_count));
}


int Magnetic_moment (spins & lattice) {
    int sum = 0;
    for (auto & i : lattice)
        for (int j : i)
            sum += j;
    return sum;
}


bool new_configuration (const double & dE, const double & T) {
    std::uniform_int_distribution<> dis (0, 1);
    double ksi = dis(gen);
    return (dE < 0 || ksi < std::exp(-(dE/T)));
}


void random_spin_flip (spins & lattice, const int & left, const int & right) {
    std::uniform_int_distribution<> dis (left, right);
    int i = dis(gen);
    int j = dis(gen);
    lattice[i][j] *= -1;
}


// Did not allocate memory - died.
void mem_allocation_2d (spins & vector, const int & dim) {
    for (auto & i : vector)
        i.resize(dim);
}


void vector_clear_2d (spins & vector) {
    for (auto & i : vector)
        i.clear();
}


template <typename T>
std::string toString (T val) {
    std::ostringstream oss;
    oss << val;
    return oss.str();
}


template<typename T, size_t... Is>
std::string tuple_to_string_impl (T const& t, std::index_sequence<Is...>) {
    return ((toString(std::get<Is>(t)) + '\t') + ...);
}

template <class Tuple>
std::string tuple_to_string (const Tuple& t) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return tuple_to_string_impl(t, std::make_index_sequence<size>{});
}


void data_file_creation (const std::string & name, std::vector<data> & exp_data) {
    std::ofstream fout;
    fout.open(name, std::ios::out | std::ios::trunc);
    for (auto & i : exp_data)
        fout << tuple_to_string(i) << std::endl;
    fout.close();
}


void plot (const std::string & name, const int & left, const int & right,
           const std::string & xlabel, const std::string & ylabel, const std::string & title) {
    std::string range = "[" + toString(left) + ":" + toString(right) + "]";
    FILE *gp = popen("gnuplot  -persist", "w");
    if (!gp) throw std::runtime_error("Error opening pipe to GNUplot.");
    std::vector<std::string> stuff = {"set term jpeg size 700, 700",
                                      "set output \'" + name + ".jpg\'",
                                      "set title \'" + title + "\'",
                                      "set xlabel" + xlabel,
                                      "set ylabel" + ylabel,
                                      "set grid xtics ytics",
                                      "set xrange " + range,
                                      "set yrange " + range,
                                      "set key off",
                                      "set ticslevel 0",
                                      "set border 4095",
                                      "plot \'" + name + ".txt\' using 1:2 w lines",
                                      "set terminal wxt",
                                      "set output",
                                      "replot", "q"};
    for (const auto & it: stuff)
        fprintf(gp, "%s\n", it.c_str());
    pclose(gp);
}