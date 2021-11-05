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


void data_file_creation (const std::string & name, std::vector<data> & exp_data);

void vector_clear_2d (spins & vector);

void mem_allocation_2d (spins & vector, const int & dim);

void plot (const std::string & name, const int & left, const int & right);


int main () {
    spins lattice (N);

    return 0;
}



int Magnetic_moment (spins & lattice) {
    int sum = 0;
    for (auto & i : lattice)
        for (int & j : i)
            sum += j;
    return sum;
}


bool new_configuration (double & dE, double & ksi, const double & T) {
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


void plot (const std::string & name, const int & left, const int & right) {
    std::string range = "[" + toString(left) + ":" + toString(right) + "]";
    FILE *gp = popen("gnuplot  -persist", "w");
    if (!gp) throw std::runtime_error("Error opening pipe to GNUplot.");
    std::vector<std::string> stuff = {"set term jpeg size 700, 700",
                                      "set output \'" + name + ".jpg\'",
                                      "set title \'P(p)\'",
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