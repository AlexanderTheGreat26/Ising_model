#include <iostream>
#include <random>
#include <vector>
#include <utility>
#include <fstream>
#include <string>
#include <tuple>
#include <array>
#include <algorithm>
#include <sstream>

typedef std::pair<int, int> index;
typedef std::pair<double, double> data;
typedef std::vector<std::vector<int>> spins;

const int N = 10;
const int left_border = 0;
const int right_border = 9;
const double q = 4;
const double J = 0;

std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()


void mem_allocation (spins & vector, const int & dim);

void lattice_filling (spins & lattice, const int & spin);

int energy_of_system (spins & lattice);

spins random_spin_flip (spins lattice);

double magnetic_moment (spins & lattice);

bool is_new_configuration (const double & dE, const double & T);

double temperature_definition (double & M);

std::vector<data> table_data_creation (std::vector<double> &E);

void data_file_creation (const std::string & name, std::vector<data> & exp_data);

void isohoric_conductivity (std::vector<double> & T, std::vector<double> & E);

void plot (const std::string & name, const std::string & xlabel, const std::string & ylabel,
           const std::string & title);


int main () {
    spins lattice(N), buff_lattice(N);
    mem_allocation(lattice, N);
    lattice_filling(lattice, 1);
    double E0 = energy_of_system(lattice);
    std::vector<double> E = {E0};
    std::vector<double> T;
    int j = 0;
    do {
        do {
            buff_lattice = std::move(random_spin_flip(lattice));
            E.emplace_back(energy_of_system(buff_lattice));
            double M = magnetic_moment(buff_lattice);
            T.emplace_back(temperature_definition(M));
            ++j;
        } while (!is_new_configuration(E[j] - E[j-1], T[j-1]) && j%5 != 0);
        lattice = buff_lattice;
    } while (j < 10000);

    isohoric_conductivity(T, E);
    std::vector<data> E_n = std::move(table_data_creation(E));
    data_file_creation("Energies", E_n);
    plot("Energies", "n", "E", "E = E(n)");
}


double mean_square (std::vector<double> & Energies) {
    int buf = 0;
    for (int i = 200; i < Energies.size(); ++i)
        buf += std::pow(Energies[i], 2);
    return buf / double (Energies.size() - 200);
}


double square_average (std::vector<double> & Energies) {
    int buf = 0;
    for (int i = 200; i < Energies.size(); ++i) // 200 is experimental data, lol.
        buf += Energies[i];
    return std::pow(double(buf) / double (Energies.size() - 200), 2);
}


void isohoric_conductivity (std::vector<double> & T, std::vector<double> & E) {
    double T2 = square_average(T);
    double c_v = 1.0 / T2 * (mean_square(E) - square_average(E));
    std::cout << "C_v\t=\t" << c_v << std::endl;
}


bool is_new_configuration (const double & dE, const double & T) {
    std::uniform_int_distribution<> dis (0, 1);
    double ksi = dis(gen);
    return (dE < 0 || ksi < std::exp(-(dE/T)));
}


double temperature_definition (double & M) {
    return std::fabs(q / atanh(M));
}


double magnetic_moment (spins & lattice) {
    int sum = 0;
    for (auto & i : lattice)
        for (int j : i)
            sum += j;
    return std::abs(double(sum) / double(std::pow(N, 2)));
}


spins random_spin_flip (spins lattice) {
    std::uniform_int_distribution<> dis (0, lattice.size()-1);
    int i = dis(gen);
    int j = dis(gen);
    lattice[i][j] *= -1;
    return lattice;
}


int energy_of_system (spins & lattice) {
    int sum = 0;
    for (int i = 0; i < lattice.size(); ++i)
        for (int j = 0; j < lattice[i].size(); ++j) {
            int S_0 = lattice[i][j];
            int S_1 = (i == N-1) ? lattice[0][j] : lattice[i+1][j];
            int S_2 = (j == N-1) ? lattice[i][0] : lattice[i][j+1];
            sum -= (S_0*S_1 + S_0*S_2);
        }
    return sum;
}


void lattice_filling (spins & lattice, const int & spin) {
    for (auto & i : lattice)
        for (int & j : i)
            j = spin;
}


void mem_allocation (spins & vector, const int & dim) {
    for (auto & i : vector)
        i.resize(dim);
}


std::vector<data> table_data_creation (std::vector<double> &E) {
    std::vector<data> result (E.size());
    for (int i = 0; i < E.size(); ++i)
        result[i] = std::make_pair(i, E[i]);
    return result;
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


void plot (const std::string & name, const std::string & xlabel, const std::string & ylabel,
           const std::string & title) {
    FILE *gp = popen("gnuplot  -persist", "w");
    if (!gp) throw std::runtime_error("Error opening pipe to GNUplot.");
    std::vector<std::string> stuff = {"set term jpeg size 700, 700",
                                      "set output \'" + name + ".jpg\'",
                                      "set title \'" + title + "\'",
                                      "set xlabel \'" + xlabel + "\'",
                                      "set ylabel \'" + ylabel + "\'",
                                      "set grid xtics ytics",
                                      "set key off",
                                      "set ticslevel 0",
                                      "set border 4095",
                                      "plot \'" + name + "\' using 1:2 w lines",
                                      "set terminal pop",
                                      "set output",
                                      "replot", "q"};
    for (const auto & it: stuff)
        fprintf(gp, "%s\n", it.c_str());
    pclose(gp);
}