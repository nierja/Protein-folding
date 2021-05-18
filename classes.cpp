//classes.cpp
#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#define SQR(x) ((x)*(x))

#include "functions.hpp"
#include "classes.hpp"

#include <assert.h>
#include <math.h>

using namespace std;

// ==================================================================

Grain::Grain() {
		pos_x = 0;
		pos_y = 0;
		pos_z = 0;
		vel_x = 0;
		vel_y = 0;
		vel_z = 0;
		acc_x = 0;
		acc_y = 0;
		acc_z = 0;
		//m   = 1;
		//TODO: constants
	}

// ==================================================================

hashtable_item::hashtable_item() {
		for (int i = 0; i < MAX_NEIGHB_LEN; i++)value[i] = -1;
}


// ==================================================================
hashtable::hashtable(double cutoff, int cube_size) {
    this->cutoff=cutoff;
    this->cube_size=cube_size;
}

void hashtable::add_item(int grain_index) {
		//adds grain with given index to the hash_table
		//grains coordinates are stored in this->grains at positions [3*i], [3*i+1], [3*i+2]
		//hash(x,i) = (hash_function(x) + i) mod k

		//printf("adding grain %d;  calling get_key(%lf, %lf, %lf)\n", grain_index, grains[3 * grain_index],\
												grains[3 * grain_index + 1], grains[3 * grain_index + 2]);
		int key = get_key(grains[3 * grain_index]-correction_x, grains[3 * grain_index + 1]-correction_y,
		 grains[3 * grain_index + 2]-correction_z, cutoff, cube_size);
		int index = hash_function(k, key);
		while (true) {
			if (hash_table[index].free) {
				hash_table[index].free = 0;
				hash_table[index].key = key;
				hash_table[index].num_of_values++;
				hash_table[index].value[0] = grain_index;
				return;
			}
			if (!hash_table[index].free && hash_table[index].key == key) {
				int i = 0;
				while (hash_table[index].value[i] != -1)i++;
				hash_table[index].num_of_values++;
				hash_table[index].value[i] = grain_index;
				assert(hash_table[index].num_of_values <= MAX_NEIGHB_LEN);
				assert(i < MAX_NEIGHB_LEN);
				return;
			}
			index = (index + 1) % k; num_of_colisions++;
		}
	}

int* hashtable::find(int key) {
		//return_array[0] == number of grains with the corresponding key, 
		//return_array[1]...return_array[1+number of neighbours] contains indexes of those grains

		int index = hash_function(k, key);
		while (!hash_table[index].free) {
			if (key == hash_table[index].key) {
				return_array[0] = hash_table[index].num_of_values;
				assert(return_array[0] <= (MAX_NEIGHB_LEN + 1));
				for (int i = 0; i < hash_table[index].num_of_values; i++)
					return_array[i + 1] = hash_table[index].value[i];
				return return_array;
			}
			index = (index + 1) % k;
		}
		return_array[0] = 0;
		return return_array;
	}

void hashtable::load_hashtable(vector <double> some_grains, int N) {
		//loads grain's coordinates from vector "some_grains" into hashtable's public atribute - vector grains
		hash_table.clear();
		k = next_prime(2 * N);
		for (int i = 0; i < k; i++)
			hash_table.push_back(hashtable_item());
		this->N = N;
		this->grains = some_grains;
		//printf("grains size: %d\n", this->grains.size());		
		for (int i = 0; i < N; i++) {
			assert(isfinite(grains[3 * i]));
			assert(isfinite(grains[3 * i + 1]));
			assert(isfinite(grains[3 * i + 2]));
			add_item(i);
		}
		assert(grains.size() == 3 * N);
		assert(hash_table.size() == k);
		/*for (int i = 0; i < N; i++) {
			printf("key %d: %d\n", i, get_key(grains[3 * i],grains[3 * i+1],grains[3 * i+2]));
		}*/
		//print hashtable
		/*printf("hashtable: \t\t\t\t(0 = free; X = occupied)\n");
		for (int i = 0; i < k; i++) {
			if (i % 50 == 0)printf("\t\n%d\t", i);
			if (hash_table[i].free == 1)printf("0");
			else printf("%d", hash_table[i].num_of_values);
		}*/
		//printf("\n\nN: %d \tk: %d \tload factor: %.3f\n\n", N, k, (float)N / k);
	}

int* hashtable::get_neighbours(int grain_index) {
		//return_array[0] == number of neighbours of given grain + 1, 
		//return_array[1]...return_array[1+number of neighbours] are indexes of neighbouring grains, contains grain_index as well!

		//printf("getting neigh. of gr. %d; get_key(%lf, %lf, %lf)\n", grain_index,grains[3 * grain_index]\
														, grains[3 * grain_index + 1], grains[3 * grain_index + 2]);
		int key = get_key(grains[3 * grain_index]-correction_x, grains[3 * grain_index + 1]-correction_y,
		 grains[3 * grain_index + 2]-correction_z, cutoff, cube_size);
		int new_key, num_of_neighbours = 0, index = 1; int* h;
		for (int i = 0; i < MAX_NEIGHBOURS; i++)neighbours[i] = 0;
		for (int i = -1; i <= 1; i += 1)
			for (int j = -1; j <= 1; j += 1)
				for (int k = -1; k <= 1; k += 1) {
					new_key = key + i + j * cube_size + k * SQR(cube_size);
					//printf("i:%d, k:%d, nk:%d\n", grain_index, key, new_key);
					h = find(new_key);
					if (h[0] != 0) {
						//printf("h[0]=%d\n", h[0]);
						num_of_neighbours += h[0];
						for (int l = 0; l < h[0]; l++) {
							neighbours[index] = h[l + 1]; index++;
						}
					}
				}
		//printf("index: %d\n", index);
		//printf("num_of_neighbours: %d |\t", num_of_neighbours);
		neighbours[0] = num_of_neighbours;
		//for (int i = 0; i < index; i++)printf("%d ", neighbours[i]);
		//printf("\n\n");
		assert(index <= N + 1);
		assert(num_of_neighbours <= N);
		return neighbours;
	}