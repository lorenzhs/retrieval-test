#pragma once

#include "util/util.h"
#include "ChunkInfo.h"
using namespace std;

/* This is the Interface a ChunkedStrategy should contain. 
 * Just for reference here, I don't want virtual method invocations. */
struct ChunkedStrategyExample {
	using Hashable_t = string;
	struct Hash {
		Hash(const string& str, uint32_t seed);
		uint32_t chunkHash() const;
	};
	struct Configuration {};
	ChunkedStrategyExample(size_t m, Configuration config);
	void addElement(const Hash& H, bool rhs);
	bool runConstruction();
	void* data();
	struct Meta {
		uint32_t sizeLogical;
		uint32_t sizeBytes;
	} meta;

	inline static bool retrieve(const Hash &H, const void* ptr, uint32_t logicalOffset, uint32_t logicalSize);
};

template<typename Strategy_T, typename ChunkInfoManager = ChunkInfoPacked>
struct RetrieverChunked {
	static const int MAX_TRIALS = 50;

	using Strategy = Strategy_T;
	using Hashable = typename Strategy::Hashable_t;
	using Hash = typename Strategy::Hash;
	
	struct Configuration {
		typename Strategy::Configuration stratConfig;
		uint32_t chunkSize;
		friend ostream& operator<<(ostream& out, const Configuration& config) {
			return out << "(C = " << config.chunkSize << ", " << config.stratConfig << ")";
		}
	} config;
	
	ChunkInfoManager chunkInfo;
	vector<uint8_t> data;

	RetrieverChunked(const vector<Hashable>& keys, const vector<bool> &values, Configuration config) : config(config) {
		assert(keys.size() == values.size());
		using KeySet = vector<int>;
		using Hashes = vector<pair<Hash, bool>>;
		uint32_t nChunks = keys.size() / config.chunkSize;
		vector<KeySet> chunkKeys(nChunks);
		vector<Hashes> chunkHashes(nChunks);
		for (size_t i = 0; i < nChunks; ++i) {
			chunkKeys[i].reserve(config.chunkSize * 1.05);
			chunkHashes[i].reserve(config.chunkSize * 1.05);
		}

		for (size_t i = 0; i < keys.size(); ++i) {
			Hash H(keys[i],0);
			uint32_t c = reduce(H.chunkHash(), nChunks);
			chunkKeys[c].push_back(i);
			chunkHashes[c].push_back({ H,values[i]});
		}
		
		vector<uint8_t> seeds(nChunks);
		vector<Strategy> solvers; solvers.reserve(nChunks);
		for (uint32_t c = 0; c < nChunks; ++c) {
			KeySet& keySet(chunkKeys  [c]);
			Hashes& hashes(chunkHashes[c]);
			assert(keySet.size() == hashes.size());
			solvers.emplace_back(keySet.size(), config.stratConfig);
			Strategy& S = solvers[c];
			uint8_t seed = 0;
			for(;;) {
				if (seed >= MAX_TRIALS) {
					return;
				}
				for (auto& p : hashes) {
					S.addElement(p.first, p.second);
				}

				if (S.runConstruction()) {
					/* free some memory */
					keySet.clear(); hashes.clear();
					break;
				}
				/* reset */
				S = Strategy(chunkKeys[c].size(), config.stratConfig);
				++seed;
				for (int i = 0; i < keySet.size(); ++i) {
					hashes[i].first = Hash(keys[keySet[i]], seed);
				}
			}
			seeds[c] = seed;
		}

		vector<uint32_t> logicalSizes; logicalSizes.reserve(nChunks);
		uint32_t totSizeBytes = 0;
		for (Strategy& S : solvers) {
			logicalSizes.push_back(S.meta.sizeLogical);
			totSizeBytes += S.meta.sizeBytes;
		}
		chunkInfo = ChunkInfoManager(logicalSizes,seeds);

		data.resize(totSizeBytes);
		totSizeBytes = 0;
		for (uint32_t c = 0; c < nChunks; ++c) {
			uint32_t s = solvers[c].meta.sizeBytes;
			memcpy(data.data() + totSizeBytes, solvers[c].data(), s);
			totSizeBytes += s;
		}
	}

	bool hasSucceeded() const {
		return data.size() > 0;
	}

	inline bool retrieve(const Hashable& s) const {
		Hash H(s, 0);
		uint32_t c = reduce(H.chunkHash(), chunkInfo.nChunks());
		ChunkInfo ci(chunkInfo[c]);
		if (ci.seed) H = Hash(s, ci.seed);
		return Strategy::retrieve(H,data.data(),ci.logicalOffset, ci.logicalSize);
	}

	size_t memoryUsage() const {
		return 8*(sizeof(this) + chunkInfo.extraMemory() + data.size());
	}

	double getAvgSeed() const {
		return chunkInfo.getAvgSeed();
	}
};
