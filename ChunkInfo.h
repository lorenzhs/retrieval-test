#pragma once
#include<vector>
#include <cstdint>
#include<cassert>
#include<numeric>
#include<algorithm>
using namespace std;

#include "util/util.h"

namespace walzer {

// Classes to Capture Meta Information on Chunks

struct ChunkInfo {
	uint32_t logicalOffset;
	uint32_t logicalSize;
	uint8_t seed;
};

struct ChunkInfoPacked {
	static inline constexpr int seedBits = 6;
	static inline constexpr int offsetBits = 26;
	struct SizeAndSeed {
		uint32_t logicalOffset : offsetBits;
		uint32_t seed : seedBits;
	};
	ChunkInfoPacked() {}
	ChunkInfoPacked(const vector<uint32_t> &sizes, const vector<uint8_t> &seeds) : v(sizes.size()+1) {
		assert(sizes.size() == seeds.size());
		uint32_t nChunks = sizes.size();

		uint32_t sumSeeds = std::accumulate(seeds.begin(),seeds.end(), uint32_t(0), std::plus<uint32_t>());
		uint8_t  maxSeed  = *std::max_element(seeds.begin(), seeds.end());
		cout << "Average Seed = " << (double)sumSeeds / nChunks << ", Largest Seed = " << (int)maxSeed << endl;
		assert(maxSeed < (1 << seedBits));

		uint32_t currOffset = 0;
		for (uint32_t c = 0; c < nChunks; ++c) {
			v[c] = { currOffset,seeds[c] };
			currOffset += sizes[c];
		}
		/* sentinel*/
		assert(currOffset < (1 << offsetBits));
		v[nChunks] = { currOffset, 0 };
	}
	ChunkInfo operator[](size_t i) const {
		ChunkInfo ret;
		ret.seed = v[i].seed;
		ret.logicalOffset = v[i].logicalOffset;
		ret.logicalSize   = v[i+1].logicalOffset-v[i].logicalOffset;
		return ret;
	}
	size_t nChunks() const {
		return v.size() - 1;
	}
	size_t extraMemory() const {
		return v.size() * sizeof(SizeAndSeed);
	}
	double getAvgSeed() const {
		int res = 0;
		for (size_t i = 0; i < nChunks(); ++i) {
			res += (*this)[i].seed;
		}
		return (double)res / nChunks();
	}
private:
	vector<SizeAndSeed> v;
};

struct ChunkInfoCompressed {
	ChunkInfoCompressed() {}
	ChunkInfoCompressed(const vector<uint32_t>& sizes, const vector<uint8_t>& seeds) {
		assert(sizes.size() == seeds.size());
		nc = sizes.size();

		uint32_t sumSeeds = std::accumulate(seeds.begin(), seeds.end(), uint32_t(0), std::plus<uint32_t>());
		uint8_t  maxSeed = *std::max_element(seeds.begin(), seeds.end());
		cout << "Average Seed = " << (double)sumSeeds / nc << ", Largest Seed = " << (int)maxSeed << endl;
		seedBits = maxSeed ? (32 - clz(maxSeed)) : 0;

		uint32_t sumSizes = std::accumulate(sizes.begin(), sizes.end(), 0);
		avgLogicalSize = sumSizes / nc;

		int32_t expOffset = 0;
		int32_t actualOffset = 0;
		int32_t maxOffsetVarPlus = 0;
		int32_t maxOffsetVarMinus = 0;
		for (uint32_t s : sizes) {
			actualOffset += s;
			expOffset += avgLogicalSize;
			int32_t diff = actualOffset - expOffset;
			maxOffsetVarPlus = max(diff, maxOffsetVarPlus);
			maxOffsetVarMinus = min(diff, maxOffsetVarMinus);
		}
		nullOffset = -maxOffsetVarMinus;
		uint32_t var = maxOffsetVarPlus - maxOffsetVarMinus;
		offsetBits = var ? 32 - clz(var) : 0;
		metaBits = offsetBits + seedBits;

		cout << "offsetBits/seedBits: " << (uint32_t)offsetBits << "/" << (uint32_t)seedBits << endl;

		expOffset = 0;
		actualOffset = 0;
		BitStreamWriter bsw;
		for (uint32_t i = 0; i < nc; ++i) {
			bsw.write(actualOffset - expOffset + nullOffset, offsetBits);
			bsw.write(seeds[i], seedBits);
			actualOffset += sizes[i];
			expOffset += avgLogicalSize;
		}
		bsw.write(actualOffset - expOffset + nullOffset, offsetBits);
		bs = bsw.finalise();
	}
	ChunkInfo operator[](size_t i) const {
		ChunkInfo ret;
		uint64_t word = bs.read(metaBits * i, metaBits + offsetBits);
		ret.logicalOffset = word & ((1 << offsetBits) - 1);
		word >>= offsetBits;
		ret.seed = word & ((1 << seedBits) - 1);
		word >>= seedBits;
		ret.logicalSize = word;

		ret.logicalSize = (ret.logicalSize - ret.logicalOffset) + avgLogicalSize;
		ret.logicalOffset += avgLogicalSize * i - nullOffset;
		return ret;
	}
	size_t nChunks() const {
		return nc;
	}
	size_t extraMemory() const {
		return bs.extraMemory();
	}
	double getAvgSeed() const {
		int res = 0;
		for (size_t i = 0; i < nChunks(); ++i) {
			res += (*this)[i].seed;
		}
		return (double)res/nChunks();
	}
private:
	BitStream bs;
	uint32_t nc;
	uint32_t avgLogicalSize;
	uint32_t nullOffset; // dieser positive Wert entspricht "Null Abweichung zur Erwartung".
	uint8_t offsetBits, seedBits, metaBits;
};

}
