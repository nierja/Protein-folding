Hamiltonian Monte carlo Simulation of a coarse grained protein model. Given a fasta protein sequence, the program is able to determine energeticaly stable structures, similar to the native folding state. Hashing, implemented in `functions.cpp` helps to reduce time complexity of the algorithm to O(N), where N is the lenght of the protein chain. Simulation parameters to be set in `folding.cpp`. Program can be compiled with `g++ -g -std=c++17 -O3 -march=native *.cpp -o folding`.
