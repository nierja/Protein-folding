//functions.cpp
#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#define SQR(x) ((x)*(x))

#include <math.h>
#include <stdlib.h>
#include <random>
#include <assert.h>

#include "functions.hpp"

using namespace std;

// ====================== potentials =========================

// some potential energy functions taken from pf.c

double u_lj(double r, double e, double s) {
	double t = s / r;
	t *= t;
	t *= t * t;
	return e * (t * t - t);
}

double du_lj(double r, double e, double s) {
	double t = s / r;
	t *= t;
	t *= t * t;
	return -6.0 * e * (2.0 * t * t - t) / r;
}

double u_angle(double r) {
	return (-3*cos(2*r+(11*M_PI)/4)+cos(3*r)-sin(r)+3)/2;
}

double du_angle(double r) {
	return (6*sin(2*r+(11*M_PI)/4)-3*sin(3*r)-cos(r))/2;
}

double u_torsion(double r) {
	return (cos(3*r+1.0/2)-3*cos(2*r-13.0/10)+cos(r+1)+3)/2;
}

double du_torsion(double r) {
	return (-3*sin(3*r+1.0/2)+6*sin(2*r-13.0/10)-sin(r+1))/2;
}

double u_ho(double r, double bond_rigidity, double bond_length) {
	return (0.5 * bond_rigidity * SQR(r - bond_length));
}

double du_ho(double r, double bond_rigidity, double bond_length) {
	return bond_rigidity * (r - bond_length);
}

double u_r(double r, double kB_T)
{				// u(r)=a/r
	return kB_T / r;
}

double du_r(double r, double kB_T)
{				// du(r)/dr=-a/(r*r)
	return -kB_T / SQR(r);
}

double u_r1(double r, double epsilon, double sigma, double cutoff) {
	return epsilon * SQR((r - cutoff) / sigma) * ((r + cutoff / 2) / sigma);
}

double du_r1(double r, double epsilon, double sigma, double cutoff) {
	return (3 * epsilon * r * (r - cutoff)) / (sigma * SQR(sigma));
}

double u_r2(double r, double epsilon, double sigma, double cutoff) {
	return epsilon * (1 - (r / sigma)) * SQR((r / sigma) - (cutoff / sigma)) / 
	((r / sigma) + SQR(r / sigma));
}

double du_r2(double r, double epsilon, double sigma, double cutoff) {
	return -(epsilon*(r-cutoff)*(r*SQR(r)+(2*sigma+cutoff)*SQR(r)+
	(-SQR(sigma)-2*cutoff*sigma)*r-cutoff*SQR(sigma)))/(sigma*SQR(r)*SQR(r+sigma));
}

// ==================== vector functions ========================

double distanc(double a[3], double b[3]) {
	// returns distance between two 3D points
	double distance_squared = 0;
	for (int i = 0; i < 3; i++)
		distance_squared += SQR(a[i] - b[i]);
	return sqrt(distance_squared);
}

void get_unitvec(double a[3], double b[3], double result[3]) {
	// calculates unit vector between points a, b
	double ab_distance = distanc(a, b);
	for (int i = 0; i < 3; i++)
		result[i] = (b[i] - a[i]) / ab_distance;
	return;
}

double unit_dotp(double unitvec_a[3], double unitvec_b[3]) {
	// returns unit dot poduct of vectors a, b:  a*b
	// equivalent to cos(phi), where phi is the angle between a, b
	double result = 0;
	for (int i = 0; i < 3; i++)
		result += unitvec_a[i] * unitvec_b[i];
	if (result < -1)result = -1;
	if (result > 1)result = 1;
	return result;
}

void unit_crsp(double unitvec_a[3], double unitvec_b[3], double result[3]) {
	// returns unit cross poduct of vectors a, b:  a X b
	double cos_ab = unit_dotp(unitvec_a, unitvec_b);
	double sin_ab = sqrt(1 - SQR(cos_ab));
	result[0] = (unitvec_a[1]*unitvec_b[2] - unitvec_a[2]*unitvec_b[1])/sin_ab;
	result[1] = (unitvec_a[2]*unitvec_b[0] - unitvec_a[0]*unitvec_b[2])/sin_ab;
	result[2] = (unitvec_a[0]*unitvec_b[1] - unitvec_a[1]*unitvec_b[0])/sin_ab;
	return;
}

double dot_product(double a[3], double b[3]) {
	// returns dot product of two vectors
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void cross_product(double a[3], double b[3], double result[3]) {
	// returns cross product of two vectors
	result[0] = a[1] * b[2] - a[2] * b[1];
	result[1] = a[2] * b[0] - a[0] * b[2];
	result[2] = a[0] * b[1] - a[1] * b[0];
	return;
}

void norm_vector(double a[3]) {
	// norms vector a
	double norm_factor = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
	for (int i = 0; i < 3; i++)a[i] /= norm_factor;
	return;
}

void get_vector(double a[3], double b[3], double ab[3]) {
	// result = a - b
	for (int i = 0; i < 3; i++)ab[i] = b[i] - a[i];
	return;
}

void subtract_vectors(double minuend[3], double subtrahend[3], double difference[3]) {
	// difference = minuend - subtrahend;
	for (int i = 0; i < 3; i++)difference[i] = minuend[i] - subtrahend[i];
	return;
}

void scale_vector(double a[3], double factor) {
	// result = a * factor
	for (int i = 0; i < 3; i++)a[i] *= factor;
	return;
}

// ====================== statistics =======================

double get_variance(int array_len, double* array) {
	//computes variance of N=array_len values, stored in array
	double average, variance, sum = 0;
	for (int i = 0; i < array_len; i++) {
		sum += array[i];
	}
	average = sum / array_len; sum = 0;
	for (int i = 0; i < array_len; i++) {
		sum += SQR(array[i] - average);
	}
	variance = sum / array_len;
	return variance;
}

double block_method(int array_len, double* E) {
	//implements block method for computing error
	double y_error, new_error;
	double* new_E = (double*)malloc(array_len * sizeof(double));
	y_error = get_variance(array_len, E) / array_len;
	while (array_len > 2) {
		array_len /= 2;
		for (int i = 0; i < array_len; i++) {
			new_E[i] = ((E[2 * i] + E[2 * i + 1]) / 2);
		}
		new_error = get_variance(array_len, new_E);
		if (new_error - y_error < sqrt(2 / (array_len - 1)) * new_error)break;
		y_error = new_error;
	}
	return y_error;
}

// ======================== hashing ==========================

int next_prime(int n) {
	//returns lowest prime, that is higher than n
	int i, is_prime;
	while (1) {
		n++; i = 2; is_prime = 1;
		while (i <= sqrt(n)) {
			if (n % i == 0) {
				is_prime = 0; break;
			}
			i++;
		}
		if (is_prime)return n;
	}
}

int get_key(double pos_x, double pos_y, double pos_z, int cell_size, int cube_size) {
	//cell size == c (cutoff distance)
	//works only for grains, whose coordinates are in interval [-500,500]
	//example: grains coordinates:[11.25,13.46,17.82] -> int:key 255_256_258
	assert(-500.0 < pos_x && pos_x < 500.0);
	assert(-500.0 < pos_y && pos_y < 500.0);
	assert(-500.0 < pos_z && pos_z < 500.0);
	int x = pos_x; if (pos_x < 0)x--;
	int y = pos_y; if (pos_y < 0)y--;
	int z = pos_z; if (pos_z < 0)z--;

	return (((x + cube_size / 2) / cell_size) + \
		((y + cube_size / 2) / cell_size) * cube_size + \
		((z + cube_size / 2) / cell_size) * SQR(cube_size));
}

int hash_function(int hashtable_size, int key) {
	//returns hash of a given key == index in a hashtable of size hashtable_size
	int m = 0.6180339887 * key;
	return m % hashtable_size;
}

// ====================== miscelaneous =======================

string get_AB_string(string residue_string) {
	// translate fasta protein format into an "AB" string
    string hydrophobic_A = "IVLPCMAG";
    string result = "";
    for (int i = 0; i < residue_string.length(); i++) {
        auto it = hydrophobic_A.find(residue_string[i]);
        if (it != std::string::npos)
            result += "A";
        else
            result += "B";
    }
    return result;
}

double uniform_random() {
	// return uniform random number with mean 0
	return (1.0 + rand()) / (2.0 + RAND_MAX);
}
