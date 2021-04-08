#include<vector>
#include<iostream>
#include<algorithm>
#include<cassert>
#include<array>
#include<tuple>
#include "util/util.h"
#include "util/line.h"
#include <cstring>
#include "util/xxhash.hpp"
using namespace std;

namespace walzer {
/*
    Variant of the OneBlockStrategy where
    • Each line contains 64*numWords many 1's
    • Obtained by copying one random 64 bit pattern numWords times.
*/
template<typename Hashable, int numWords = 1>
struct OneBlockStrategy {
    static constexpr char stratName[] = "OneBlock";
    using Hashable_t = Hashable;
    using Pattern = uint64_t;
    using WORDS = array<uint64_t, numWords + 1>;

    inline static WORDS toWords(uint32_t start, const Pattern& p) {
        WORDS words;
        words[0] = p;
        using hash_t = xxh::hash64_t;
        if constexpr (numWords > 1) {
            // generate 64 more random bits
            words[1] = xxh::xxhash3<64>({ p }, uint64_t(1));
            // generate the rest with double hashing. why? because it works.
            for (uint64_t i = 2; i < numWords; ++i) {
                words[i] = words[i - 1] + words[0];
            }
        }
        uint64_t mask = (1LL << (start & 63)) - 1;
        words[numWords] = words[0] & mask;
        words[0] &= ~mask;
        return words;
    }

    /* Scalar product of two bitvectors.
     * First one: Dense vector of 64-bit words "line".
     * Second one: Given as array<uint64_t,l> and value start:
                   start/64 zero-words is implicitely prepended
                   and a suitable number of zero-words appended. */
    inline static bool mul(uint32_t start, const WORDS& words, const uint64_t* line) {
        int o = start / 64;
        uint64_t prod(0);
        for (int i = 0; i <= numWords; ++i) {
            prod ^= words[i] & line[o + i];
        }
        return parityll(prod);
    }

    /* Scalar Product of a bitvector given as 64-bit words "line"
       and a bitvector given by a repeated pattern shifted by "start" positions to the right.
       Except, it is not quite shifted. */
    inline static bool mul(uint32_t start, const Pattern& p, const uint64_t* line) {
        return mul(start, toWords(start, p), line);
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
        uint32_t nVar = logicalSize - 64 * numWords + 1;
        uint32_t start = reduce(H.start, nVar);
        return mul(start + logicalOffset, H.p, (uint64_t*)ptr);
    }

    OneBlockStrategy(size_t m, Configuration config) {
        nVar = m / config.c;
        meta.sizeLogical = nVar + 64 * numWords - 1;
        /* round up to full words: */
        meta.sizeLogical = (meta.sizeLogical + 63) & ~63L;
        meta.sizeBytes = meta.sizeLogical / 8;
        nVar = meta.sizeLogical - 64 * numWords + 1;
    }
    void addElement(const Hash& H, bool rhs) {
        uint32_t start = reduce(H.start, nVar);
        eqs.push_back({ start,toWords(start,H.p),rhs });
    }

    bool runConstruction(bool silent = true) {
        {
            silent || cout << "Counting Sorting... "; cout.flush();
            vector<int> count(nVar);
            for (const Equation& eq : eqs) count[eq.start]++;
            for (int i = 1; i < nVar; count[i] += count[i - 1], i++);
            vector<Equation> sorted(eqs.size());
            for (const Equation& eq : eqs) sorted[--count[eq.start]] = eq;
            swap(eqs, sorted);
            silent || cout << "DONE" << endl;
        }

        uint32_t additions = 0;
        vector<uint32_t> piv(eqs.size());
        for (int i = 0; i < eqs.size(); ++i) {
            /* determine word and offset of the first 1-entry (the pivot) */
            uint32_t k = 0;
            WORDS& W = eqs[i].words;
            while (k < numWords + 1 && !W[k]) k++;
            if (k == numWords + 1) {
                silent || cout << "FAIL" << endl;
                return false;
            }
            piv[i] = (eqs[i].start & ~63) + k * 64 + ctzll(W[k]);
            assert(piv[i] >= eqs[i].start);
            uint32_t pivotWord = piv[i] / 64;
            uint64_t pivotMask = 1ULL << (piv[i] & 63);

            /* eliminate from all subsequent rows */
            for (int j = i + 1; j < eqs.size() && eqs[j].start <= piv[i]; ++j) {
                uint32_t wordOff = eqs[j].start / 64;
                WORDS& Wj = eqs[j].words;
                assert(pivotWord - wordOff >= 0);
                assert(pivotWord - wordOff <= numWords);
                /* test if there is a 1-entry in pivot pos */
                uint32_t kj = pivotWord - wordOff;
                if (pivotMask & Wj[kj]) {
                    /* need to add W[k...numWords] to Wj[kj...] */
                    for (uint32_t kk = 0; kk <= numWords - k; ++kk) {
                        assert(kj + kk <= numWords);
                        Wj[kj + kk] ^= W[k + kk];
                        additions++;
                    }
                    /* don't forget the rhs */
                    eqs[j].rhs ^= eqs[i].rhs;
                }
            }
        }
        silent || cout << "SUCCESS, avg additions / per: " << (additions / (double)eqs.size()) << endl;

        /* back substitution */
        sol = Line(meta.sizeLogical);
        for (int i = eqs.size() - 1; i >= 0; --i) {
            bool prod = mul(eqs[i].start, eqs[i].words, sol.v.data()) ^ eqs[i].rhs;
            sol.setBit(piv[i], prod);
        }

        return true;
    }
};
}
