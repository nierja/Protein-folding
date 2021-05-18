#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#define SQR(x) ((x)*(x))

// Turn ON/OFF target potentials
#define H_MC 1
#define ANGLE_BOND 0
#define ANGLE_BOND_SCALAR 0
#define ANGLE_BOND_DELTA 1
#define ANGLE_TORS 0
#define ANGLE_TORS_SCALAR 0
#define ANGLE_TORS_DELTA 1
#define NONBONDING 1
#define BONDING 1

#include "functions.hpp"
#include "classes.hpp"

#include <assert.h>
#include <fenv.h>
#include <vector>
#include <stdio.h>
#include <random>
#include <numeric>

using namespace std;
default_random_engine generator;
normal_distribution<double> distribution(0.0, 1.0);


//========================= constants and parameters ========================
double dt = 0.005;                           // velocity - Verlet step lenght[s * (m / e) * *(1 / 2)]
double m_grain = 1;                          // weight of single grain[m]
double e = 1;                                // LJ potential constant [e]
double s = 1;                                // LJ potential constant [s]
double cutoff = 5 * s;                       // cutoff distance[s]
int equi_steps = 0, steps = 5000, num_steps; // number of MC steps
double kB_T = 1;                             // kB_T == k_B * T[e]
double bond_length = 1;                      // [s]
double bond_rigidity = 10;					 // bond rigidity constant
//constexpr int MAX_NEIGHB_LEN = 92;		 // maximal number of grains stored inside one cell
//constexpr int MAX_NEIGHBOURS = 92;		 / maximal number of neighbouring grains
int cube_size = 1000;						 // base for calculating key of a grain
double delta_1 = 0;
double delta_2 = 0;
double k_1 = -1, k_2 = 0.5; 
double r_AA = 3.8;

//==================================================== AA sequences =====================================================
string fasta_4RXN = "MKKYTCTVCGYIYDPEDGDPDDGVNPGTDFKDIPDDWVCPLCGVGKDEFEEVEE";			// 54
string fasta_2GB1 = "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE";			// 56
string fasta_1IGD = "MTPAVTTYKLVINGKTLKGETTTKAVDAETAEKAFKQYANDNGVDGVWTYDDATKTFTVTE";	// 61
string fasta_2YGS = "MDAKARNCLLQHREALEKDIKTSYIMDHMISDGFLTISEEEKVRNEPTQQQRAAMLIK"		// 92
					"MILKKDNDSYVSFYNALLHEGYKDLAALLHDGIP";								
string fasta_1GAB = "TIDQWLLKNAKEDAIAELKKAGITSDFYFNAINKAKTVEEVNALKNEILKAHA";			// 53
string fasta_1BDD = "TADNKFNKEQQNAFYEILHLPNLNEEQRNGFIQSLKDDPSQSANLLAEAKKLNDAQAPKA";		// 60
string fasta_1E0L = "GATAVSEWTEYKTADGKTYYYNNRTLESTWEKPQELK";							// 37
string fasta_1CLB = "MKSPEELKGIFEKYAAKEGDPNQLSKEELKLLLQTEFPSLLKGGSTLDELFEELDKNGDGE"		// 75
					"VSFEEFQVLVKKIS";	
string fasta_1LQ7 = "GSRVKALEEKVKALEEKVKALGGGGRIEELKKKWEELKKKIEELGGGGEVKKVEEEVKKLE"		// 67
					"EEIKKL";	
string fasta_1E0G = "DSITYRVRKGDSLSSIAKRHGVNIKDVMRWNSDTANLQPGDKLTLFVK";					// 48
string fasta_3M0R = "MDETGKELVLALYDYQEKSPDEVTMKKGDILTLLNSTNKDWWKVEVNDRQGFVPAAYV";		// 58
string tors = "AAAA";																	// 4
string C_v = "MDAKARNCLLQHREALEKDI";														// 20


class Protein {
	/*class implementig protein model*/
public:
	int N = 0;									// protein chain length (coarse-grained)
	int sequence_length = 0;
	string sequence;							// A-B sequence of grain types
	vector <Grain> grains;						// container for individual grains
	vector <double> angles;						// N-2 bond and N-3 torsion angles
	vector<double> positions;					// stores 3N grains coordinates
	hashtable my_hashtable = hashtable(cutoff, cube_size);

	void update_hashtable() {
		// fills vector "positions" with 3N grains coordinates and passes them
		// into the load_hashtable method
		double min_x=grains[0].pos_x;
		double min_y=grains[0].pos_y;
		double min_z=grains[0].pos_z;

		positions.clear();
		for (int i = 0; i < N; i++) {
			if (grains[i].pos_x < min_x)min_x = grains[i].pos_x;
			if (grains[i].pos_y < min_y)min_y = grains[i].pos_y;
			if (grains[i].pos_z < min_z)min_z = grains[i].pos_z;
			positions.push_back(grains[i].pos_x);
			positions.push_back(grains[i].pos_y);
			positions.push_back(grains[i].pos_z);
		}
		//printf("\nloading hashtable...\n");
		my_hashtable.correction_x = min_x;
		my_hashtable.correction_y = min_y;
		my_hashtable.correction_z = min_z;
		my_hashtable.load_hashtable(positions, N);


		//printf("\n");

	}

	double distance(int grain_i, int grain_j) {
		// returns distance between grain_i, grain_j
		return sqrt(SQR(grains[grain_i].pos_x - grains[grain_j].pos_x) +
			SQR(grains[grain_i].pos_y - grains[grain_j].pos_y) +
			SQR(grains[grain_i].pos_z - grains[grain_j].pos_z));
	}

	double H() {
		/* V = sum(nonbonding E) + sum(bonding E) + sum(tors. angle E) + sum(bond. angle E)
		   T = sum(kinetic energy) */
		double r, V = 0, T = 0, E_cut = u_lj(cutoff, e, s);
		double a[3], b[3], c[3], d[3], ab[3], ba[3], bc[3], cd[3], result[3];
		int num_of_neighbours = 0, neighbour; int* neighbours;
#if NONBONDING		
		update_hashtable();
		//printf("getting neighbours...\n");
		for (int i = 0; i < N; i++) {
			//for (int neighbour = 0; neighbour < N; neighbour++) {
			neighbours = my_hashtable.get_neighbours(i);
			assert(neighbours[0] <= N);
			for (int j = 1; j <= neighbours[0]; j++) {
				neighbour = neighbours[j];
				//printf("H: grain %d and neighbour %d, r=%.3f, ", i, neighbour, distance(i, neighbour));
				//printf("%d.x = %.2f, %d.x = %.2f\t", i, grains[i].pos_x, neighbour, grains[neighbour].pos_x);
				if (i > neighbour) {
					//printf("%d > %d: YES\n", i, neighbour);
					r = distance(i, neighbour);
					if (r >= cutoff)continue;
					if (sequence[i] == 'A' && sequence[neighbour] == 'A') {
						V += u_lj(r, e, s) - E_cut;
						//V += u_r2(r);
					}
					else {
						V += 0.5 * (u_lj(r, e, s) - E_cut);
						//V += 0.5 * u_r2(r);
					}
				}
				//else printf("%d > %d: NO\n", i, neighbour);
			}
		}
		//printf("\n");
#endif
#if BONDING
		for (int i = 0; i < N - 1; i++) {
			r = distance(i, i + 1);
			V += u_ho(r, bond_rigidity, bond_length);
		}
#endif
#if ANGLE_BOND
		update_angles();
		for (int i = 0; i < N - 2; i++)
			V += u_angle(angles[i]);
#endif
#if ANGLE_BOND_SCALAR
		update_angles();
		for (int i = 0; i < N-2; i++) {
			a[0] = grains[i].pos_x; 	a[1] = grains[i].pos_y; 	a[2] = grains[i].pos_z;
			b[0] = grains[i + 1].pos_x; b[1] = grains[i + 1].pos_y; b[2] = grains[i + 1].pos_z;
			c[0] = grains[i + 2].pos_x; c[1] = grains[i + 2].pos_y; c[2] = grains[i + 2].pos_z;

			get_vector(a, b, ab);
			get_vector(b, c, bc);
			//printf("dotp: %lf\n", dot_product(ab, bc));
			
			V -= k_1 * dot_product(ab, bc);
		}
#endif
#if ANGLE_BOND_DELTA
		// decrements V with [k_1 * (sin(delta)cos(beta)-(cos(delta)(1-(dot_product(ab, bc)^2))))]
        update_angles();
        for (int i = 0; i < N-2; i++) {
            a[0] = grains[i].pos_x;     a[1] = grains[i].pos_y;     a[2] = grains[i].pos_z;
            b[0] = grains[i + 1].pos_x; b[1] = grains[i + 1].pos_y; b[2] = grains[i + 1].pos_z;
            c[0] = grains[i + 2].pos_x; c[1] = grains[i + 2].pos_y; c[2] = grains[i + 2].pos_z;

            get_vector(a, b, ab);
            get_vector(b, c, bc);

            V -= k_1 * ((sin(delta_1) * dot_product(ab, bc)) - (cos(delta_1) * (1-SQR(dot_product(ab, bc)))));
        }
#endif
#if ANGLE_TORS	
		update_angles();
		for (int i = 0; i < N - 3; i++)
			V -= 0.5 * u_torsion(angles[N - 2 + i]);
#endif
#if ANGLE_TORS_SCALAR
        update_angles();
        for (int i = 0; i < N-3; i++) {
            a[0] = grains[i].pos_x; a[1] = grains[i].pos_y; a[2] = grains[i].pos_z;
            b[0] = grains[i + 1].pos_x; b[1] = grains[i + 1].pos_y; b[2] = grains[i + 1].pos_z;
            c[0] = grains[i + 2].pos_x; c[1] = grains[i + 2].pos_y; c[2] = grains[i + 2].pos_z;
            d[0] = grains[i + 3].pos_x; d[1] = grains[i + 3].pos_y; d[2] = grains[i + 3].pos_z;

            get_vector(a, b, ab);
            get_vector(c, d, cd);

            V -= k_2 * dot_product(ab, cd);
        }
#endif
#if ANGLE_TORS_DELTA
		// decrements V with [k_1 * (sin(delta)cos(beta)-(cos(delta)(1-(dot_product(ab, bc)^2))))]
        update_angles();
        for (int i = 0; i < N-3; i++) {
            a[0] = grains[i].pos_x;     a[1] = grains[i].pos_y;     a[2] = grains[i].pos_z;
            b[0] = grains[i + 1].pos_x; b[1] = grains[i + 1].pos_y; b[2] = grains[i + 1].pos_z;
            c[0] = grains[i + 2].pos_x; c[1] = grains[i + 2].pos_y; c[2] = grains[i + 2].pos_z;
			d[0] = grains[i + 3].pos_x; d[1] = grains[i + 3].pos_y; d[2] = grains[i + 3].pos_z;

            get_vector(a, b, ab);
            get_vector(c, d, cd);

            V += k_2 * ((sin(delta_2) * dot_product(ab, cd)) - (cos(delta_2) * (1-SQR(dot_product(ab, cd)))));
        }
#endif

		for (int i = 0; i < N; i++) {
			//printf("%d.vel_x=%.2f,", i, grains[i].vel_x); assert(isfinite(grains[i].vel_x));
			T += 0.5 * m_grain * (SQR(grains[i].vel_x) + SQR(grains[i].vel_y) + SQR(grains[i].vel_z));
		}
		//printf("T: %f\n", T);
		//printf("V: %f\n", V);
		assert(isfinite(T + V));
		return T + V;
	}

	double r_gyr() {
		// returns mean of the second power of the gyration radius <R_gyr^2>
		double r_cm[] = {0, 0, 0};
		double r_i[] = {0, 0, 0};
		double r_gyr = 0, dist;
		for (int i = 0; i < N; i++) {
			r_cm[0] += grains[i].pos_x*r_AA;
			r_cm[1] += grains[i].pos_y*r_AA;
			r_cm[2] += grains[i].pos_z*r_AA;
		}
		for (int i = 0; i < 3; i++) {
			r_cm[i] /= N;
		}
		for (int i = 0; i < N; i++) {
			r_i[0] = grains[i].pos_x*r_AA;
			r_i[1] = grains[i].pos_y*r_AA;
			r_i[2] = grains[i].pos_z*r_AA;
			dist = distanc(r_i, r_cm);
			r_gyr += SQR(dist);
		}
		return r_gyr / N;
	}

	double r_ee() {
		// returns mean of the second power of end to end distance <R_ee^2>
		double sum = 0;
		for (int i = 0; i < N; i++) {
			sum += sqrt(SQR(grains[i].pos_x*r_AA) +
						SQR(grains[i].pos_y*r_AA) +
						SQR(grains[i].pos_z*r_AA));
		}
		return (sum / (N - 1));
	}
	
	Protein(string primary_struct) {
		//my_hashtable = hashtable(cutoff, cube_size);
		double init_coordinates[3], random_unitvec[3];
		double H_init, H_final;

		sequence = primary_struct;
		sequence_length = sequence.length();
		assert(sequence_length > 1);
		angles.reserve(sequence_length);	//to prevent reallocation
		grains.reserve(sequence_length);
		for (int i = 0; i < 2 * sequence_length - 5; i++) {
			angles.push_back(0);
		}
		// random shooting configuration
		for (int i = 0; i < sequence_length; i++) {
			if (i == 0) {
				grains.push_back(Grain());
				N++;
			}
			else {
				grains.push_back(Grain());
				N++;
				grains[i].pos_x = grains[i-1].pos_x + 1;
				grains[i].pos_y = grains[i-1].pos_y;
				grains[i].pos_z = grains[i-1].pos_z;
				for (int j = 0; j < 100; j++) {
					//save coordinates, measure energy
					init_coordinates[0] = grains[i].pos_x;
					init_coordinates[1] = grains[i].pos_y;
					init_coordinates[2] = grains[i].pos_z;
					H_init = H();
					for (int k = 0; k < 3; k++)
						random_unitvec[k] = distribution(generator);
					norm_vector(random_unitvec);
					for (int k = 0; k < 3; k++)
						random_unitvec[k] *= bond_length;
					grains[i].pos_x = grains[i-1].pos_x + random_unitvec[0];
					grains[i].pos_y = grains[i-1].pos_y + random_unitvec[1];
					grains[i].pos_z = grains[i-1].pos_z + random_unitvec[2];
					H_final = H();
					if (H_final > H_init) {
						grains[i].pos_x = init_coordinates[0];
						grains[i].pos_y = init_coordinates[1];
						grains[i].pos_z = init_coordinates[2];
					}
				}

			}
		}
		// starting configuration is a slightly deformed rod
		/*for (int i = 0; i < N; i++) {
			grains[i].pos_x = i/3+0.3*uniform_random();
			grains[i].pos_y = (i+1)/3+0.3*uniform_random();
			grains[i].pos_z = (i+2)/3+0.3*uniform_random();
		}*/
	}

	void restore() {
		double init_coordinates[3], random_unitvec[3];
		double H_init, H_final;
		N = 0;
		grains.clear();
		// random shooting configuration
		for (int i = 0; i < sequence_length; i++) {
			if (i == 0) {
				grains.push_back(Grain());
				N++;
			}
			else {
				N++;
				grains.push_back(Grain());
				grains[i].pos_x = grains[i-1].pos_x + 1;
				grains[i].pos_y = grains[i-1].pos_y;
				grains[i].pos_z = grains[i-1].pos_z;
				for (int j = 0; j < 100; j++) {
					//save coordinates, measure energy
					init_coordinates[0] = grains[i].pos_x;
					init_coordinates[1] = grains[i].pos_y;
					init_coordinates[2] = grains[i].pos_z;
					H_init = H();
					for (int k = 0; k < 3; k++)
						random_unitvec[k] = distribution(generator);
					norm_vector(random_unitvec);
					for (int k = 0; k < 3; k++)
						random_unitvec[k] *= bond_length;
					grains[i].pos_x = grains[i-1].pos_x + random_unitvec[0];
					grains[i].pos_y = grains[i-1].pos_y + random_unitvec[1];
					grains[i].pos_z = grains[i-1].pos_z + random_unitvec[2];
					H_final = H();
					if (H_final > H_init) {
						grains[i].pos_x = init_coordinates[0];
						grains[i].pos_y = init_coordinates[1];
						grains[i].pos_z = init_coordinates[2];
					}
				}

			}
		}
	}

	void to_xyz(const char * name) {
		//creates folding.xyz file with protein folding trajectory 
		double scale = 3.8;
		FILE* fp;
		fp = fopen(name, "a");
		fprintf(fp, "%d\n", N);
		fprintf(fp, "comment line\n");
		assert(grains.size() == N);
		for (int i = 0; i < N; i++) {
			if (sequence[i] == 'A')fprintf(fp, "CA %lf %lf %lf\n", grains[i].pos_x*scale,
			 		grains[i].pos_y*scale, grains[i].pos_z*scale);
			else fprintf(fp, "CA %lf %lf %lf\n", grains[i].pos_x*scale,
					grains[i].pos_y*scale, grains[i].pos_z*scale);
		}
		fclose(fp);
	}

	void to_csv(const char * name) {
		//creates .csv file with angle pairs
		FILE* fp;
		fp = fopen(name, "a");
		for (int i = 0; i < N-3; i++) {
			fprintf(fp, "%lf;%lf\n", angles[i]*(180/M_PI),
							 angles[N - 2 + i]*(180/M_PI));
		}
		fclose(fp);
	}

	void best_to_pdb(const char * name) {
		//creates .csv file with angle pairs
		FILE* fp;
		fp = fopen(name, "w");
		fprintf(fp, "CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1\n");
		for (int i = 0; i < N; i++) {
			fprintf(fp, "ATOM    %.3d  CA      X   1      %.3lf   %.3lf  %.3lf  0.00  0.00          CA\n",
			i+1, grains[i].pos_x*r_AA, grains[i].pos_y*r_AA, grains[i].pos_z*r_AA);
		}
		fprintf(fp, "CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1\n");
		fclose(fp);
	}


	void update_angles() {
		// stores 2N - 5 angles computed from XYZ coordinates into vector "angles"
		double a[3], b[3], c[3], d[3];			// point's coordinates
		double unitvec_ba[3], unitvec_bc[3], unitvec_cb[3], unitvec_cd[3];
		double n_1[3], n_2[3]; 
		double dot_prod;

		// update bond angles from interval (0, pi)
		for (int i = 0; i < N - 2; i++) {
			a[0] = grains[i].pos_x; a[1] = grains[i].pos_y; a[2] = grains[i].pos_z;
			b[0] = grains[i+1].pos_x; b[1] = grains[i+1].pos_y; b[2] = grains[i+1].pos_z;
			c[0] = grains[i+2].pos_x; c[1] = grains[i+2].pos_y; c[2] = grains[i+2].pos_z;
			get_unitvec(b, a, unitvec_ba);
			get_unitvec(b, c, unitvec_bc);
			dot_prod = unit_dotp(unitvec_ba, unitvec_bc);
			angles[i] = acos(dot_prod);
		}
		
		// update torsion angles from interval (-pi, pi)
		for (int i = 0; i < N - 3; i++) {
			a[0] = grains[i].pos_x; a[1] = grains[i].pos_y; a[2] = grains[i].pos_z;
			b[0] = grains[i+1].pos_x; b[1] = grains[i+1].pos_y; b[2] = grains[i+1].pos_z;
			c[0] = grains[i+2].pos_x; c[1] = grains[i+2].pos_y; c[2] = grains[i+2].pos_z;
			d[0] = grains[i + 3].pos_x; d[1] = grains[i + 3].pos_y; d[2] = grains[i + 3].pos_z;

			// computes unit vector n_1 perpendicular to the plane containing a, b, c
			get_unitvec(b, a, unitvec_ba);
			get_unitvec(b, c, unitvec_bc);
			unit_crsp(unitvec_ba, unitvec_bc, n_1);

			// computes unit vector n_2 perpendicular to the plane containing b, c, d
			get_unitvec(c, b, unitvec_cb);
			get_unitvec(c, d, unitvec_cd);
			unit_crsp(unitvec_cb, unitvec_cd, n_2);

			dot_prod = unit_dotp(n_1, n_2);
			// get the correct sign
			if (unit_dotp(n_1, unitvec_cd) < 0)
				angles[N - 2 + i] = -acos(dot_prod);
			else
				angles[N - 2 + i] = acos(dot_prod);
		}
	}
};


class Simulation : public Protein {
private:
	double* E = (double*)malloc(steps * sizeof(double));
	double* R_gyr = (double*)malloc(steps * sizeof(double));
	double* R_ee = (double*)malloc(steps * sizeof(double));
	double E_best = H();
	int accepted = 0;                           // number of accepted configurations
	int num_steps = 1;							// number of Verlet steps


	void generate_velocities() {
		// initial velocities are generated from normal distribution
		for (int i = 0; i < N; i++) {
			grains[i].vel_x = sqrt((2*kB_T) / (3*m_grain)) * distribution(generator);
			grains[i].vel_y = sqrt((2*kB_T) / (3*m_grain)) * distribution(generator);
			grains[i].vel_z = sqrt((2*kB_T) / (3*m_grain)) * distribution(generator);
		}
	}

	void update_coordinates() {
		//updates beads coordinates
		for (int i = 0; i < N; i++) {
			grains[i].pos_x += grains[i].vel_x * dt + 0.5 * grains[i].acc_x * SQR(dt);
			grains[i].pos_y += grains[i].vel_y * dt + 0.5 * grains[i].acc_y * SQR(dt);
			grains[i].pos_z += grains[i].vel_z * dt + 0.5 * grains[i].acc_z * SQR(dt);
		}
		//printf("update coordinates: ");
		for (int i = 0; i < N; i++) {
			//printf("%d,", i);
			assert(isfinite(grains[i].pos_x));
			assert(isfinite(grains[i].pos_y));
			assert(isfinite(grains[i].pos_z));
			/*if (grains[i].pos_x > 500)grains[i].pos_x = 450;
			if (grains[i].pos_y > 500)grains[i].pos_y = 450;
			if (grains[i].pos_z > 500)grains[i].pos_z = 450;
			if (grains[i].pos_x < -500)grains[i].pos_x = -450;
			if (grains[i].pos_y < -500)grains[i].pos_y = -450;
			if (grains[i].pos_z < -500)grains[i].pos_z = -450;*/

		}//printf("\n");
		return;
	}

	void update_velocity() {
		// updates beads velocity
		for (int i = 0; i < N; i++) {
			grains[i].vel_x += 0.5 * grains[i].acc_x * dt;
			grains[i].vel_y += 0.5 * grains[i].acc_y * dt;
			grains[i].vel_z += 0.5 * grains[i].acc_z * dt;
		}
		//printf("update velocity: ");
		for (int i = 0; i < N; i++) {
			//printf("%d,", i);
			assert(isfinite(grains[i].vel_x));
			assert(isfinite(grains[i].vel_y));
			assert(isfinite(grains[i].vel_z));
		}//printf("\n"); 
		return;
	}

	void get_coordinates(vector <double>& x, vector <double>& y, vector <double>& z) {
		// stores current grain's coordinates to vectors x, y, z
		x.clear();
		y.clear();
		z.clear();
		for (int i = 0; i < N; i++) {
			x.push_back(grains[i].pos_x);
			y.push_back(grains[i].pos_y);
			z.push_back(grains[i].pos_z);
		}return;
	}

	void load_coordinates(vector <double>x, vector <double>y, vector <double>z) {
		// updates protein coordinates with values stored in vectors x, y, z
		for (int i = 0; i < N; i++) {
			grains[i].pos_x = x[i];
			grains[i].pos_y = y[i];
			grains[i].pos_z = z[i];
		}return;
	}


	void update_acceleration() {
		// computes acceleration of each grain
		double r, f, a_x, a_y, a_z;					// distance, force and accelerations
		int k;										// hashtable size
		int neighbour;								// index of a neighbouring grain
		int* neighbours;							// points to an array returned by get_neighbours()
		double a[3], b[3], c[3], d[3], o[3];			// point's coordinates
		double unitvec_ab[3], unitvec_ba[3], unitvec_bc[3], unitvec_cb[3], unitvec_cd[3];
		double p_1[3], p_2[3], f_a[3], f_b[3], f_c[3], f_d[3];
		double oc[3], cd[3], ba[3], ab[3], cb[3], bc[3], t_1[3], t_2[3], t_3[3], t_4[3], t_c[3]; // t-x ... storage vectors
		double ab_1[3], ab_2[3], bc_1[3], bc_2[3];
		double acc_ax, acc_ay, acc_az, acc_cx, acc_cy, acc_cz, acc_dx, acc_dy, acc_dz;
		double k2 = 0.5, factor, dU;

		for (int i = 0; i < N; i++)
			grains[i].acc_x = grains[i].acc_y = grains[i].acc_z = 0;
#if NONBONDING
		update_hashtable();
		for (int i = 0; i < N; i++) {
			//for (int neighbour = 0; neighbour < N; neighbour++) {
			neighbours = my_hashtable.get_neighbours(i);
			assert(neighbours[0] <= N);
			for (int j = 1; j <= neighbours[0]; j++) {
				neighbour = neighbours[j];
				//printf("Acc: grain %d and neighbour %d, r=%.3f, ", i, neighbour, distance(i, neighbour));
				//printf("%d.x = %.2f, %d.x = %.2f\t", i, grains[i].pos_x, neighbour, grains[neighbour].pos_x);
				if (i > neighbour) {
					//printf("%d > %d: YES\n", i, neighbour);
					r = distance(i, neighbour);
					if (r >= cutoff)continue;
					//printf("*\n");
					if (sequence[i] == 'A' && sequence[neighbour] == 'A') {
						f = -du_lj(r, e, s);
						//f = -du_r2(r);
					}
					else {
						f = -0.5 * du_lj(r, e, s);
						//f = -0.5 * du_r2(r);
					}
					a_x = f * (grains[i].pos_x - grains[neighbour].pos_x) / (r * m_grain);
					a_y = f * (grains[i].pos_y - grains[neighbour].pos_y) / (r * m_grain);
					a_z = f * (grains[i].pos_z - grains[neighbour].pos_z) / (r * m_grain);
					//printf("force: %.3f.\tacc_x_ %.3f,\tycc_y: %.3f,\tycc_z: %.3f\n", f, a_x, a_y, a_z);
					grains[i].acc_x += a_x;
					grains[i].acc_y += a_y;
					grains[i].acc_z += a_z;
					grains[neighbour].acc_x -= a_x;
					grains[neighbour].acc_y -= a_y;
					grains[neighbour].acc_z -= a_z;
				}
				//else printf("%d > %d: NO\n", i, neighbour);
			}
		}
#endif
#if BONDING
		for (int i = 0; i < N; i++)
			for (int j = 0; j < 2; j++) {
				if (j == 0)k = i - 1;
				else k = i + 1;
				if (k < 0 || k > N - 1)continue;
				r = distance(i, k);
				f = -du_ho(r, bond_rigidity, bond_length);
				grains[i].acc_x += f * (grains[i].pos_x - grains[k].pos_x) / (r * m_grain);
				grains[i].acc_y += f * (grains[i].pos_y - grains[k].pos_y) / (r * m_grain);
				grains[i].acc_z += f * (grains[i].pos_z - grains[k].pos_z) / (r * m_grain);
			}
#endif
#if ANGLE_BOND

		update_angles();
		for (int i = 0; i < N - 2; i++) {
			a[0] = grains[i].pos_x; a[1] = grains[i].pos_y; a[2] = grains[i].pos_z;
			b[0] = grains[i + 1].pos_x; b[1] = grains[i + 1].pos_y; b[2] = grains[i + 1].pos_z;
			c[0] = grains[i + 2].pos_x; c[1] = grains[i + 2].pos_y; c[2] = grains[i + 2].pos_z;
			
			// grain a=i
			get_unitvec(b, a, unitvec_ba);
			get_unitvec(b, c, unitvec_bc);
			unit_crsp(unitvec_ba, unitvec_bc, p_1);
			unit_crsp(unitvec_ba, p_1, p_2);
			f = (-du_angle(angles[i]) / distance(i, i+1));
			acc_ax = f * p_2[0] / m_grain; grains[i].acc_x += acc_ax;
			acc_ay = f * p_2[1] / m_grain; grains[i].acc_y += acc_ay;
			acc_az = f * p_2[2] / m_grain; grains[i].acc_z += acc_az;

			// grain c=i+1
			get_unitvec(c, b, unitvec_cb);
			unit_crsp(unitvec_ba, unitvec_bc, p_1);
			unit_crsp(unitvec_cb, p_1, p_2);

			f = (-du_angle(angles[i]) / distance(i + 1, i + 2));
			acc_cx = f * p_2[0] / m_grain; grains[i + 2].acc_x += acc_cx;
			acc_cy = f * p_2[1] / m_grain; grains[i + 2].acc_y += acc_cy;
			acc_cz = f * p_2[2] / m_grain; grains[i + 2].acc_z += acc_cz;

			// grain b=i+1
			// !!! will colapse when a, c have different weights
			grains[i + 1].acc_x -= (acc_ax + acc_cx);
			grains[i + 1].acc_y -= (acc_ay + acc_cy);
			grains[i + 1].acc_z -= (acc_az + acc_cz);
		}

#endif
#if ANGLE_BOND_SCALAR
		update_angles();
		for (int i = 0; i < N - 2; i++) {
			a[0] = grains[i].pos_x; a[1] = grains[i].pos_y; a[2] = grains[i].pos_z;
			b[0] = grains[i + 1].pos_x; b[1] = grains[i + 1].pos_y; b[2] = grains[i + 1].pos_z;
			c[0] = grains[i + 2].pos_x; c[1] = grains[i + 2].pos_y; c[2] = grains[i + 2].pos_z;
			
			get_vector(a, b, ab);
			get_vector(b, a, ba);
			get_vector(b, c, bc);
			get_vector(c, b, cb);

			// grain a=i
			grains[i].acc_x -= cb[0] / m_grain;
			grains[i].acc_y -= cb[1] / m_grain;
			grains[i].acc_z -= cb[2] / m_grain;

			// grain c=i+2
			grains[i + 2].acc_x -= ab[0] / m_grain;
			grains[i + 2].acc_y -= ab[1] / m_grain;
			grains[i + 2].acc_z -= ab[2] / m_grain;

			// grain b=i+1
			grains[i + 1].acc_x -= (ba[0] + bc[0]) / m_grain;
			grains[i + 1].acc_y -= (ba[1] + bc[1]) / m_grain;
			grains[i + 1].acc_z -= (ba[2] + bc[2]) / m_grain;
		}
#endif
#if ANGLE_BOND_DELTA
        update_angles();
        for (int i = 0; i < N - 2; i++) {
			// f = [k_1 * (sin(delta)cos(beta)-(cos(delta)(1-(dot_product(ab, bc)^2))))];
			// w = cos(beta) ~ dot_product(ab, bc) -> df/dw = factor = 
			//   = k_1 * (2*cos(delta)*w + sin(delta);
			// dw/da = cb; dw/db = ba + bc; dw/dc = ab; 

            a[0] = grains[i].pos_x; a[1] = grains[i].pos_y; a[2] = grains[i].pos_z;
            b[0] = grains[i + 1].pos_x; b[1] = grains[i + 1].pos_y; b[2] = grains[i + 1].pos_z;
            c[0] = grains[i + 2].pos_x; c[1] = grains[i + 2].pos_y; c[2] = grains[i + 2].pos_z;

            get_vector(a, b, ab);
            get_vector(b, a, ba);
            get_vector(b, c, bc);
            get_vector(c, b, cb);
            factor = -k_1*((2*cos(delta_1)*dot_product(ab, bc)) + sin(delta_1));
			//dU = k_1*(sin(delta) * (dot_product(ab, bc)+0.01) - (cos(delta) * (1-SQR(dot_product(ab, bc)-0.01))))/0.02;
			//printf("factor: %lf\t num: %lf\n", factor, dU);

            // grain a=i
            grains[i].acc_x += (bc[0]*factor) / m_grain;
            grains[i].acc_y += (bc[1]*factor) / m_grain;
            grains[i].acc_z += (bc[2]*factor) / m_grain;

            // grain c=i+2
            grains[i + 2].acc_x += (ba[0]*factor) / m_grain;
            grains[i + 2].acc_y += (ba[1]*factor) / m_grain;
            grains[i + 2].acc_z += (ba[2]*factor) / m_grain;

            // grain b=i+1
            grains[i + 1].acc_x += ((ab[0] + cb[0])*factor) / m_grain;
            grains[i + 1].acc_y += ((ab[1] + cb[1])*factor) / m_grain;
            grains[i + 1].acc_z += ((ab[2] + cb[2])*factor) / m_grain;
        }
#endif
#if ANGLE_TORS	
		update_angles();
		for (int i = 0; i < N - 3; i++) {
			a[0] = grains[i].pos_x; a[1] = grains[i].pos_y; a[2] = grains[i].pos_z;
			b[0] = grains[i + 1].pos_x; b[1] = grains[i + 1].pos_y; b[2] = grains[i + 1].pos_z;
			c[0] = grains[i + 2].pos_x; c[1] = grains[i + 2].pos_y; c[2] = grains[i + 2].pos_z;
			d[0] = grains[i + 3].pos_x; d[1] = grains[i + 3].pos_y; d[2] = grains[i + 3].pos_z;
			for (int j = 0; j < 3; j++)o[j] = (b[j] + c[j]) / 2;	//o is center of bc bond

			// grain a - should be OK
			get_unitvec(b, a, unitvec_ba);
			get_unitvec(b, c, unitvec_bc);
			unit_crsp(unitvec_ba, unitvec_bc, p_1);

			f = (k2 * -du_torsion(angles[N - 2 + i])) / (distance(i, i + 1) * sin(angles[i]));
			for (int j = 0; j < 3; j++)f_a[j] = f * p_1[j];
			acc_ax = f_a[0] / m_grain; grains[i].acc_x += acc_ax;
			acc_ay = f_a[1] / m_grain; grains[i].acc_y += acc_ay;
			acc_az = f_a[2] / m_grain; grains[i].acc_z += acc_az;

			// grain d - should be OK
			get_unitvec(c, b, unitvec_cb);
			get_unitvec(c, d, unitvec_cd);
			unit_crsp(unitvec_cd, unitvec_cb, p_2);

			f = (k2 * -du_torsion(angles[N - 2 + i])) / (distance(i + 2, i + 3) * sin(angles[i + 1]));
			for (int j = 0; j < 3; j++)f_d[j] = f * p_2[j];
			acc_dx = f_d[0] / m_grain; grains[i + 3].acc_x += acc_dx;
			acc_dy = f_d[1] / m_grain; grains[i + 3].acc_y += acc_dy;
			acc_dz = f_d[2] / m_grain; grains[i + 3].acc_z += acc_dz;

			// grain c
			// compute t_c (ad. https://arxiv.org/pdf/1401.1181.pdf; eq. (46))
			get_vector(o, c, oc);
			get_vector(c, d, cd);
			get_vector(b, a, ba);
			
			for (int j = 0; j < 3; j++) {
				cd[j] *= 0.5;
				ba[j] *= 0.5;
			}
			cross_product(oc, f_d, t_1);
			cross_product(cd, f_d, t_2);
			cross_product(ba, f_a, t_3);
			for (int j = 0; j < 3; j++)
				t_c[j] = -(t_1[j] + t_2[j] + t_3[j]);
			cross_product(t_c, oc, t_4);
			for (int j = 0; j < 3; j++)
				f_c[j] = (1 / SQR(distanc(o, c)))*t_4[j];
			acc_cx = f_c[0] / m_grain; grains[i + 2].acc_x += acc_cx;
			acc_cy = f_c[1] / m_grain; grains[i + 2].acc_y += acc_cy;
			acc_cz = f_c[2] / m_grain; grains[i + 2].acc_z += acc_cz;

			// finally: grain b; will fail with different weights
			grains[i + 1].acc_x -= (acc_ax + acc_cx + acc_dx);
			grains[i + 1].acc_y -= (acc_ay + acc_cy + acc_dy);
			grains[i + 1].acc_z -= (acc_az + acc_cz + acc_dz);
		}
#endif
#if ANGLE_TORS_SCALAR
		update_angles();
		for (int i = 0; i < N-3; i++) {
			a[0] = grains[i].pos_x; a[1] = grains[i].pos_y; a[2] = grains[i].pos_z;
			b[0] = grains[i + 1].pos_x; b[1] = grains[i + 1].pos_y; b[2] = grains[i + 1].pos_z;
			c[0] = grains[i + 2].pos_x; c[1] = grains[i + 2].pos_y; c[2] = grains[i + 2].pos_z;
			d[0] = grains[i + 3].pos_x; d[1] = grains[i + 3].pos_y; d[2] = grains[i + 3].pos_z;

			get_vector(a, b, ab);
			get_vector(c, d, cd);
			scale_vector(ab, -k_2);
			scale_vector(cd, -k_2);
			
			// grain a=i
			grains[i].acc_x += cd[0] / m_grain;
			grains[i].acc_y += cd[1] / m_grain;
			grains[i].acc_z += cd[2] / m_grain;
			
			// grain b=i+1
			grains[i + 1].acc_x -= cd[0] / m_grain;
			grains[i + 1].acc_y -= cd[1] / m_grain;
			grains[i + 1].acc_z -= cd[2] / m_grain;

			// grain c=i+2
			grains[i + 2].acc_x += ab[0] / m_grain;
			grains[i + 2].acc_y += ab[1] / m_grain;
			grains[i + 2].acc_z += ab[2] / m_grain;

			// grain d=i+3
			grains[i + 3].acc_x -= ab[0] / m_grain;
			grains[i + 3].acc_y -= ab[1] / m_grain;
			grains[i + 3].acc_z -= ab[2] / m_grain;
		}
#endif
#if ANGLE_TORS_DELTA
        update_angles();
		for (int i = 0; i < N-3; i++) {
			a[0] = grains[i].pos_x; a[1] = grains[i].pos_y; a[2] = grains[i].pos_z;
			b[0] = grains[i + 1].pos_x; b[1] = grains[i + 1].pos_y; b[2] = grains[i + 1].pos_z;
			c[0] = grains[i + 2].pos_x; c[1] = grains[i + 2].pos_y; c[2] = grains[i + 2].pos_z;
			d[0] = grains[i + 3].pos_x; d[1] = grains[i + 3].pos_y; d[2] = grains[i + 3].pos_z;

			get_vector(a, b, ab);
            get_vector(c, d, cd);
            factor = -k_2*((2*cos(delta_2)*dot_product(ab, cd)) + sin(delta_2));
			//dU = k_2*(sin(delta) * (dot_product(ab, cd)+0.01) - (cos(delta) * (1-SQR(dot_product(ab, cd)-0.01))))/0.02;
			//printf("factor: %lf\t num: %lf\n", factor, dU);
			
			// grain a=i
			grains[i].acc_x -= (cd[0]*factor) / m_grain;
			grains[i].acc_y -= (cd[1]*factor) / m_grain;
			grains[i].acc_z -= (cd[2]*factor) / m_grain;
			
			// grain b=i+1
			grains[i + 1].acc_x += (cd[0]*factor) / m_grain;
			grains[i + 1].acc_y += (cd[1]*factor) / m_grain;
			grains[i + 1].acc_z += (cd[2]*factor) / m_grain;

			// grain c=i+2
			grains[i + 2].acc_x -= (ab[0]*factor) / m_grain;
			grains[i + 2].acc_y -= (ab[1]*factor) / m_grain;
			grains[i + 2].acc_z -= (ab[2]*factor) / m_grain;

			// grain d=i+3
			grains[i + 3].acc_x += (ab[0]*factor) / m_grain;
			grains[i + 3].acc_y += (ab[1]*factor) / m_grain;
			grains[i + 3].acc_z += (ab[2]*factor) / m_grain;
		}
#endif

		//printf("update acceleration: ");
		for (int i = 0; i < N; i++) {
			//printf("%d.acc_x=%f,", i, grains[i].acc_x);
			assert(isfinite(grains[i].acc_x));
			assert(isfinite(grains[i].acc_y));
			assert(isfinite(grains[i].acc_z));
		}
		//printf("\n");
	}

	void Verlet() {
		//printf("\n\n\tVERLET\n\n");
		for (int i = 0; i < num_steps; i++) {
			update_coordinates();
			update_velocity();
			update_acceleration();
			update_velocity();
		}
		return;
	}

	void Metropolis(int step_number) {
		vector <double> X, Y, Z;					// temporary storage of grain coordinates
		double H1, H2;
		FILE *energies_fp = fopen("energies.txt", "a");

		if (H_MC)generate_velocities();
		H1 = H();
		get_coordinates(X, Y, Z);
		Verlet();
		H2 = H();
		if (exp(-(H2 - H1) / kB_T) > uniform_random()) {
			accepted += 1;
			fprintf(energies_fp, "%lf\n", H2);
			if (H2 < E_best) {
				E_best = H2;
				to_xyz("best.xyz");
				best_to_pdb("best.pdb");
			}
			E[step_number] = H2;
		}
		else {
			load_coordinates(X, Y, Z);
			fprintf(energies_fp, "%lf\n", H1); 
			if (H1 < E_best) {
				E_best = H1;
				to_xyz("best.xyz");
				best_to_pdb("best.pdb");
			}
			E[step_number] = H1;
		}
		if ((step_number % 100) == 0) {
			to_xyz("folding.xyz");
		}
		fclose(energies_fp); 
	}

public:
	Simulation(string primary_struct) : Protein(primary_struct) {
	}

	void simulate(bool MC = 0, bool SA = 0) {
		// type of simulation is set by #define constants & function parameters
		double r, variance, sum_r_gyr, sum_r_ee;
		int watch_index = 1, counter = 0;			// index of grains printed to console
		FILE *variances_fp, *r_gyr_fp, *r_ee_fp;

		// delete previous output files and open new ones
		remove("best.xyz"); remove("variances.txt"); remove("angles.csv");
		remove("energies.txt"); remove("folding.xyz"); remove("r_gyr.txt");
		remove("r_ee.txt");
		variances_fp = fopen("variances.txt", "a");	
		r_gyr_fp = fopen("r_gyr.txt", "a");	
		r_ee_fp = fopen("r_ee.txt", "a");	
		to_xyz("best.xyz"); best_to_pdb("best.pdb");

		if (H_MC) {
				num_steps = 100;
			}
		
		if (MC) {
			for (int i = 0; i < equi_steps; i++) {
				Metropolis(i);
			}
			for (int i = 0; i < steps; i++) {
				Metropolis(i);
				r = distance(0, watch_index);
				if (i % 100 == 0) {
					to_csv("angles.csv");
				}
				if (i % 100 == 0) {
					printf("========================================================="
						"===================\n(H)MC step n. %d x: %lf, vel_x: %lf, acc_x:"
						" %lf, r %lf,\n         E_best: %lf, E_tot: %.12lf.\n", i, grains[watch_index].pos_x,
						grains[watch_index].vel_x, grains[watch_index].acc_x, r, E_best, H());
					printf("Hashing corrections: %lf %lf %lf\n", my_hashtable.correction_x, 
						my_hashtable.correction_y, my_hashtable.correction_z);
					printf("r_gyr = %lf\n", r_gyr());
					printf("r_ee = %lf\n", r_ee());
					fprintf(r_gyr_fp, "%lf\n", r_gyr()); 
					fprintf(r_ee_fp, "%lf\n", r_ee()); 
				}
			}
			printf("p_acc: %.12lf\n", (double)accepted / steps);
		}

		if (SA) {
			double start_temperature = 10 * kB_T;
			double final_temperature = 0.1 * kB_T;
			double cooling_constant = 0.98;

			kB_T = start_temperature;
			while (kB_T >= final_temperature) {	// SA loop
				for (int i = 0; i < steps; i++) {
					Metropolis(i);
					if (i % 100 == 0 && kB_T < 0.18) {
						to_csv("angles.csv");
					}
					R_gyr[i] = r_gyr();
					R_ee[i] = r_ee();
				}
				variance = get_variance(steps, E);
				//printf("%lf %lf\n", kB_T, variance);
				sum_r_gyr = sum_r_ee = 0;
				for (int i = 0; i < steps; i++) {
					sum_r_gyr += R_gyr[i];
					sum_r_ee += R_ee[i];
				}
				fprintf(variances_fp, "%lf %lf\n", kB_T, (1/SQR(kB_T))*variance);
				fprintf(r_gyr_fp, "%lf %lf %lf\n", kB_T, (sum_r_gyr/steps)); 
				fprintf(r_ee_fp, "%lf %lf\n", kB_T, (sum_r_ee/steps));
				printf("<R_gyr^2> = %lf\n", (sum_r_gyr/steps)); 
				printf("<R_ee^2> = %lf\n", (sum_r_ee/steps));
				printf("r_gyr_error: %lf\t r_ee_error: %lf\n", block_method(steps, R_gyr), block_method(steps, R_ee));
				kB_T *= cooling_constant;
				counter++;
				r = distance(0, watch_index);
				printf("==========================================================="
					   "==================\nSA step n. %d (kB_T=%.3lf) x: %lf, vel_x:"
					   " %lf, acc_x: %lf, \n        r: %lf, E_best: %lf, E_tot: %lf.\n",
					   counter, kB_T, grains[watch_index].pos_x, grains[watch_index].vel_x,
					   grains[watch_index].acc_x, r, E_best, H());
				//restore();
			}
			printf("p_acc: %.12lf\n", (double)accepted / (steps*counter));
		}

		fclose(variances_fp);
		fclose(r_gyr_fp);
		fclose(r_ee_fp);
		return;
	}
};


void xyz_to_angles(const char * file_name) {
	Simulation protein(get_AB_string(fasta_1CLB));
	char str[50], str2[50];
   	int number_of_atoms;
	double coordinates[3];
   	FILE *fp;

   	fp = fopen (file_name, "r");
	fscanf(fp, "%d\n", &number_of_atoms);
	fscanf(fp, "%s\n", str);
	printf("%s %s\n", str);
	printf("%d\n", protein.N);
	for (int j = 0; j < 1; j++) {
		for (int i = 0; i < protein.N; i++) {
			fscanf(fp, "%s %lf %lf %lf", str, &coordinates[0], &coordinates[1], &coordinates[2]);
			printf("%lf %lf %lf\n", coordinates[0], coordinates[1], coordinates[2]);
			protein.grains[i].pos_x = coordinates[0]*r_AA;
			protein.grains[i].pos_y = coordinates[1]*r_AA;
			protein.grains[i].pos_z = coordinates[2]*r_AA;
		}
		protein.update_angles();
		protein.to_csv("pdb_angles.csv");
	}
	fclose(fp);
	return;
}

int main(int argc, char* argv[]) {
	/*	depending on the type of simulation, may create following files:
		angles.csv		- collection of (bonding angle, tors. angle) pairs
							throughout simulation for creating "Ramchadran"
		best.pdb		- configuration with the lowest energy in pdb format
		best.xyz 		- collection of configurations with the least energy
							throughout simulation
		energies.txt 	- collected energies of configurations
		folding.xyz 	- trajectory of simulation
		r_ee.txt		- stores <R_ee^2>
		r_gyr.txt		- stores <R_gyr^2>
		variances.txt 	- only during SA, 
							collection of (1/SQR(kB_T))*variance for C_v(T).
	*/
	//feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );

	if (argc == 5) {
		printf("Taking command line arguments !\n\n");
		delta_1 = atof(argv[1]); printf("delta_1: %lf\n", delta_1);
		delta_2 = atof(argv[2]); printf("delta_2: %lf\n", delta_2);
		k_1 = atof(argv[3]); printf("k_1: %lf\n", k_1);
		k_2 = atof(argv[4]); printf("k_2: %lf\n", k_2);
	}

	Simulation protein(get_AB_string(fasta_1CLB));
	protein.simulate(0, 1);
	
	// Options for generating output:

	//xyz_to_angles("1clb.xyz");
	//system("vmd folding.xyz");
	//system("vmd best.xyz");
	//system("./TMalign best.pdb 3m0r.pdb");					// prints rmsd
	//system("gnuplot -p -e \"plot 'energies.txt'  w l\"");
	//system("gnuplot -p -e \"plot 'variances.txt' w l\"");
	//system("gnuplot -p -e \"plot 'r_ee.txt'      w l\"");
	//system("gnuplot -p -e \"plot 'r_gyr.txt'     w l\"");
	return 0;
}