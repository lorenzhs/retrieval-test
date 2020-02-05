#include "util/util.h"

using namespace std;

template<typename Strategy>
struct Retriever {
	using Hashable = typename Strategy::Hashable_t;
	int seed;
	typename Strategy::Solution sol;
	using Hash = typename Strategy::Hash;
	using Configuration = typename Strategy::Configuration;
	Configuration config;

	Retriever(const vector<Hashable> &keys, const vector<bool> &values, Configuration config) : seed(0), config(config) {
		assert(keys.size() == values.size());
		for(;;++seed) {
			Strategy S(keys.size(), config);
			for (int i = 0; i < keys.size(); ++i) {
				S.addElement(Hash(keys[i], seed), values[i]);
			}
			if (S.runConstruction()) {
				swap(sol, S.sol);
				break;
			}
		}
	}

	inline bool retrieve(const Hashable &s) const {
		return Strategy::retrieve(sol, Hash(s, seed));
	}

	size_t memoryUsage() const {
		return sizeof(this) + Strategy::numBits(sol);
	}
};
