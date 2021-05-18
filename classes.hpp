// classes.hpp

#ifndef CLASSES_HPP  
#define CLASSES_HPP 

#include <string>
#include <vector>

using namespace std;

constexpr int MAX_NEIGHB_LEN = 92;	     // maximal number of grains stored inside one cell
constexpr int MAX_NEIGHBOURS = 92;


class Grain {
	/* class for storing position, velocity
	  and acceleration of an indiviual grain */
public:
	double pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, acc_x, acc_y, acc_z;
	Grain(); 
};

class hashtable_item {
	/* building unit of a hashtable */
public:
	int free = 1;								// free == 0 -> contains some grains; else empty
	int key = 0;
	int num_of_values = 0;						// number of stored grains
	int value[MAX_NEIGHB_LEN];					// stores indexes of neighbouring grains
	hashtable_item();							//value[i] = -1 -> no grain index stored in value[i]
};

class hashtable {
	/* implemets hashtable for cell list hashing in O(n) time and space complexity */
public:
	double correction_x, correction_y, correction_z;
private:
	vector <double> grains;								// stores 3N grains coordinates
	vector<hashtable_item> hash_table;					// hashtable itself
	int k, N, num_of_colisions = 0, initialized = 0;
	int return_array[MAX_NEIGHB_LEN + 1];				// results from function find() stored inside
	int neighbours[MAX_NEIGHBOURS];						// results from function get_neighbours() stored inside
    double cutoff;
    int cube_size;

	void add_item(int grain_index);
	int* find(int key);

public:
    hashtable(double cutoff, int cube_size);
	void load_hashtable(vector <double> some_grains, int N);
	int* get_neighbours(int grain_index);
};


#endif