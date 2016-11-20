#ifndef CHROMOSOME_H_
#define CHROMOSOME_H_

#include "Util.h"

class Chromosome {
 public:
 	Chromosome(Random& rand, size_t length) {
 		solution = rand_vector(rand, length);
 	}
	vector<bool>& getSolution() {
		return solution;
	}
	Chromosome* getUplink() {
		return uplink;
	}
	Chromosome* getDownlink() {
		return downlink;
	}

	static void createLink(Chromosome* down, Chromosome* up) {
		down->uplink = up;
		up->downlink = down;
	}

 private:
 	vector<bool> solution;
	Chromosome* downlink;
	Chromosome* uplink;
};
#endif

