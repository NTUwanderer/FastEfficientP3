#ifndef CHROMOSOME_H_
#define CHROMOSOME_H_

#include "Util.h"

class Chromosome {
 public:
 	Chromosome(Random& rand, size_t length):
		uplink(0),
		downlink(0) {
 		solution = rand_vector(rand, length);
 	}

	const vector<bool> & getSolution() const {
		return solution;
	}

	vector<bool> & getSolution() {
		return solution;
	}

	bool&& operator[](std::size_t idx) {
		return solution[idx];
	}

	const bool& operator[](std::size_t idx) const {
		return solution[idx];
	}

	void setSolution(vector<bool>& _solution) {
		solution = _solution;
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

