// functions.hpp

#ifndef FUNCTIONS_HPP  
#define FUNCTIONS_HPP 

#include <string>

using namespace std;

// ====================== potentials ========================
// some functions for potential energy calculation taken from pf.c

double u_lj(double r, double e, double s);
double du_lj(double r, double e, double s);
double u_angle(double r);
double du_angle(double r);
double u_torsion(double r);
double du_torsion(double r);
double u_ho(double r, double bond_rigidity, double bond_length);
double du_ho(double r, double bond_rigidity, double bond_length);
double u_r(double r, double kB_T);
double du_r(double r, double kB_T);
double u_r1(double r, double epsilon, double sigma, double cutoff);
double du_r1(double r, double epsilon, double sigma, double cutoff);
double u_r2(double r, double epsilon, double sigma, double cutoff);
double du_r2(double r, double epsilon, double sigma, double cutoff);

// ================== vector functions ======================

double distanc(double a[3], double b[3]);
void get_unitvec(double a[3], double b[3], double result[3]);
double unit_dotp(double unitvec_a[3], double unitvec_b[3]);
void unit_crsp(double unitvec_a[3], double unitvec_b[3], double result[3]);
double dot_product(double a[3], double b[3]);
void cross_product(double a[3], double b[3], double result[3]);
void norm_vector(double a[3]);
void get_vector(double a[3], double b[3], double ab[3]);
void subtract_vectors(double minuend[3], double subtrahend[3], double difference[3]);
void scale_vector(double a[3], double factor);

// ======================= statistics ========================

double get_variance(int array_len, double* array);
double block_method(int array_len, double* E);

// ======================== hashing ==========================

int next_prime(int n);
int get_key(double pos_x, double pos_y, double pos_z, int cell_size, int cube_size);
int hash_function(int hashtable_size, int key);

// ====================== miscelaneous =======================

string get_AB_string(string residue_string);
double uniform_random();


#endif