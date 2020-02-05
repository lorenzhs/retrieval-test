#include<vector>
#include<iostream>
#include<algorithm>
#include<cassert>
#include<array>
#include<tuple>
#include "util/util.h"
#include "util/line.h"

using namespace std;

/* Todo: Make l a template parameter
 * requires generic implementation of BlockSolver::Hash */

/* Todo: We assume that the logical size (number of bits in solution) is a multiple of 64.
 * Either make that more transparent or get rid of this assumption.
 * When not getting rid of it, make better use of it, for instance,
 * let the new logical size be the old one divided by 64. */

template<typename Hashable>
struct OneBlockStrategy {
	using Hashable_t = Hashable;
	static constexpr int l = 64;
    /* number of 64-bit words used per line */
    static constexpr size_t l_words = (l + 2*63)/64;
    static constexpr size_t l_bytes = (l + 7) / 8;
    using Pattern   = array<uint8_t,l_bytes>;
    using WORDS     = array<uint64_t,l_words>;
    
	/* Scalar product of two bitvectors.
	 * First one: Dense vector of 64-bit words "line".
	 * Second one: Given as array<uint64_t,l> and value start:
	               start/64 zero-words is implicitely prepended
				   and a suitable number of zero-words appended. */
    inline static bool mul(uint32_t start, const WORDS &words, const uint64_t* line) {
        int o = start / 64;
        uint64_t prod(0);
        for(int i = 0; i < l_words; ++i) {
            prod ^= words[i] & line[o+i];
        }
        return parityll(prod);
    }
    
	/* Assume Pattern is written into an array of 64 bit words
	   and then shifted by start positions. We return those
	   64-bit words that overlap with the shifted pattern. */
    inline static WORDS toWords(uint32_t start, const Pattern &p) {
        WORDS words; words.fill(0LL);
        memcpy(&words,&p,l_bytes);
        int o = start & 63;
        uint64_t mask = (1LL << o) - 1;
        uint64_t rest = words[0] & mask;
        words[0] &= ~mask;
        unsigned targetWord = (o < l ? l : o) / 64;
        unsigned targetPos  = (o < l ? l : o) & 63;
        words[targetWord] |= rest << targetPos;
        if (targetPos + min(o,l) > 64) {
            words[targetWord+1] = rest >> (64-(targetPos?targetPos:64)); // compiler error even though targetPos == 0 can't reach this */
        }
        return words;
    }
    
	/* Scalar Product of a bitvector given as 64-bit words "line"
	   and a bitvector given by a pattern that should be shifted by "start" positions to the right. */
    inline static bool mul(uint32_t start, const Pattern &p, const uint64_t* line) {
        return mul(start,toWords(start,p),line);
    }

	struct Equation {
		uint32_t start;
		WORDS words;
		bool rhs;
	};

	uint32_t nVar;
	vector<Equation> eqs;
	Line sol;

	/* Todo: Make Generic, i.e. work with l > 64 */
	struct Hash {
		uint32_t ch, start;
		Pattern p;
		Hash(const Hashable& s, uint32_t seed) {
			/* Todo: Get rid of this requirement */
			static_assert(sizeof(Hash) == 128 / 8);
			do_hash(s, seed, this);
		};
		uint32_t chunkHash() const { return ch; }
	};

	struct Configuration {
		double c;
		friend ostream& operator<<(ostream& out, const Configuration& config) {
			return out << "(c = " << config.c << ")";
		}
	};

	void* data() {
		return sol.v.data();
	};
	struct Meta {
		uint32_t sizeLogical;
		uint32_t sizeBytes;
	} meta;

	inline static bool retrieve(const Hash& H, const void* ptr, uint32_t logicalOffset, uint32_t logicalSize) {
		uint32_t nVar = logicalSize - l + 1;
		uint32_t start = H.start % nVar;
		return mul(start + logicalOffset, toWords(start,H.p), (uint64_t*)ptr);
	}
	
	OneBlockStrategy(size_t m, Configuration config) {
		nVar = m / config.c;
		meta.sizeLogical = nVar + l - 1;
		/* round up to full words: */
		meta.sizeLogical = (meta.sizeLogical + 63) & ~63L;
		meta.sizeBytes = meta.sizeLogical / 8;
		nVar = meta.sizeLogical - l + 1;
	}
    void addElement(const Hash& H, bool rhs) {
		uint32_t start = H.start % nVar;
        eqs.push_back({start,toWords(start,H.p),rhs});
    }
    
    bool runConstruction(bool silent = true) {
        {
            silent || cout << "Counting Sorting... "; cout.flush();
            vector<int> count(nVar);
            for(const Equation &eq: eqs) count[eq.start]++;
            for(int i = 1; i < nVar; count[i] += count[i-1], i++);
            vector<Equation> sorted(eqs.size());
            for(const Equation &eq: eqs) sorted[--count[eq.start]] = eq;
            swap(eqs,sorted);
            silent || cout << "DONE" << endl;
        }
        
        uint32_t additions = 0;
        vector<uint32_t> piv(eqs.size());
        for(int i = 0; i < eqs.size(); ++i) {
            /* determine word and offset of the first 1-entry (the pivot) */
            uint32_t k = 0;
            WORDS &W = eqs[i].words;
            while(k < l_words && !W[k]) k++;
            if (k == l_words) {
                silent || cout << "FAIL" << endl;
                return false;
            }
            piv[i]  = (eqs[i].start&~63) + k*64 + ctzll(W[k]);
            assert(piv[i] >= eqs[i].start);
            uint32_t pivotWord = piv[i] / 64;
            uint64_t pivotMask = 1ULL << (piv[i] & 63);
            
            /* eliminate from all subsequent rows */
            for(int j = i+1; j < eqs.size() && eqs[j].start <= piv[i]; ++j) {
                uint32_t wordOff = eqs[j].start / 64;
                WORDS &Wj = eqs[j].words;
                assert(pivotWord-wordOff >= 0);
                assert(pivotWord-wordOff < l_words);
                /* test if there is a 1-entry in pivot pos */
                uint32_t kj = pivotWord-wordOff;
                if (pivotMask & Wj[kj]) {
                    /* need to add W[k...l_words-1] to Wj[kj...] */
                    for(uint32_t kk = 0; kk < l_words-k; ++kk) {
                        assert(kj+kk < l_words);
                        Wj[kj+kk] ^= W[k+kk];
                        additions++;
                    }
                    /* don't forget the rhs */
                    eqs[j].rhs ^= eqs[i].rhs;
                }
            }
        }
        silent || cout << "SUCCESS, avg additions / per: " << (additions/(double)eqs.size()) << endl;
        
        /* back substitution */
        sol = Line(meta.sizeLogical);
        for(int i = eqs.size()-1; i >= 0; --i) {
            bool prod = mul(eqs[i].start,eqs[i].words,sol.v.data()) ^ eqs[i].rhs;
            sol.setBit(piv[i],prod);
        }
        
        return true;
    }
};