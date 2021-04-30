#include "util/util.h"

using namespace std;

namespace walzer {

/* for benchmarking it sucks if one strategy
 * got "unlucky" and has to restart while others got lucky.
 * allow some cheating. */
int GAMED_SEED = 0;

template<typename Strategy_T>
struct Retriever {
	static const int MAX_TRIALS = 3;

	using Strategy = Strategy_T;
	using Hashable = typename Strategy::Hashable_t;

	int seed;
	typename Strategy::Solution sol;
	using Hash = typename Strategy::Hash;
	using Configuration = typename Strategy::Configuration;
	Configuration config;

	Retriever(Configuration config) : seed(GAMED_SEED), config(config) {}

	Retriever(const vector<Hashable>& keys, const vector<bool>& values,
                  Configuration config) : seed(GAMED_SEED), config(config) {
		Construct(keys, values);
	}

        template <typename data_t = bool>
	void Construct(const vector<Hashable> &keys, const vector<data_t> &values)  {
		assert(keys.size() == values.size());
		for(;;++seed) {
			if (seed >= MAX_TRIALS) {
				return;
			}
			Strategy S(keys.size(), config);
			for (size_t i = 0; i < keys.size(); ++i) {
				S.addElement(Hash(keys[i], seed), values[i]);
			}
			if (S.runConstruction()) {
                            std::swap(sol, S.sol);
				break;
			}
		}
	}

	bool hasSucceeded() const {
		return Strategy::numBits(sol) > 0;
	}

	inline auto retrieve(const Hashable &s) const {
		return Strategy::retrieve(sol, Hash(s, seed));
	}

	size_t memoryUsage() const {
		return sizeof(this) + Strategy::numBits(sol);
	}

	double getAvgSeed() const {
		return seed - GAMED_SEED; // one seed to rule them all
	}
};
}
