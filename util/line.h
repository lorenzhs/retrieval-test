#ifndef LINE_H
#define LINE_H

#include<vector>
#include<cassert>
using namespace std;

/* todo: properly understand move semantics */

template<typename T>
vector<T> operator^(const vector<T> &v1, const vector<T> &v2) {
    vector<T> v3(v1);
    for(int i = 0; i < (int)v2.size(); ++i) {
        v3[i] ^= v2[i];
    }
    return std::move(v3);
}
template<typename T>
vector<T>& operator^=(vector<T> &v1, const vector<T> &v2) {
    for(int i = 0; i < (int)v2.size(); ++i) {
        v1[i] ^= v2[i];
    }
    return v1;
}

struct Line;
void swap(Line &l1, Line &l2);

struct Line {
    vector<uint64_t> v;
    int nBits;
    
    Line(int nBits, uint64_t bits) : nBits(nBits) {
        assert(nBits <= 64);
        resize(nBits);
        v[0] = bits;
    };
    
    
    /* move constructor */
    Line(Line &&other) noexcept {
        swap(*this, other);
    };
    Line& operator=(Line &&other) noexcept {
        swap(*this, other);
        return *this;
    }
    Line() : v(2), nBits(0) {}
    Line(int nBits, vector<uint64_t> &&v) : v(v), nBits(nBits) {};
    Line(int nBits) : v((nBits >> 6) + 2), nBits(nBits) {};
    inline void resize(int nBitsNew) {
        nBits = nBitsNew;
        v.resize((nBits >> 6) + 2);
        /* in case of shrinking, make zeroes */
        v[nBits >> 6] &= ((1LL << (nBits & 63)) - 1);
        v[(nBits >> 6)+1] = 0;
    }
    inline int size() const {
        return nBits;
    }
    bool any() const {
        for(int i = 0; i <= (nBits >> 6); ++i) {
            if(v[i]) return true;
        }
        return false;
    }
    inline bool getBit(int i) const {
        return (v[i >> 6] >> (i & 63)) & 1;
    }
    inline void setBit(int i, bool val) {
        uint64_t mask = 1LL << (i & 63);
        v[i >> 6] = (v[i >> 6] & ~mask) | ((uint64_t)val << (i & 63));
    }
    inline void addToBit(int i, bool val) {
        v[i >> 6] ^= (uint64_t)val << (i & 63);
    }
    inline bool operator[](int i) const {
        return getBit(i);
    }    
    Line copy() const {
        Line L;
        L.v = v;
        L.nBits = nBits;
        return L;
    }
    void copyFrom(const Line &other) {
        v = other.v;
        nBits = other.nBits;
    }
    bool operator==(const Line &other) const {
        if(other.nBits != nBits) return false;
        for(int i = 0; i <= (nBits >> 6); ++i) {
            if (other.v[i] != v[i])
                return false;
        }
        return true;
    }
    Line operator^(const Line &other) const {
        assert(other.nBits == nBits);
        return Line { nBits, v ^ other.v};
    }
    Line& operator^=(const Line &other) {
        if (nBits < other.nBits)
            resize(other.nBits);
        assert(v.size() >= other.v.size());
        v ^= other.v;
        return *this;
        // assert(other.nBits == nBits);
        // v ^= other.v;
        // return *this;
    }
    void addFrom(const Line &other, int knownZeroes) {
        knownZeroes >>= 6;
        for(int i = knownZeroes; i < (int)v.size(); ++i) {
            v[i] ^= other.v[i];
        }
    }
    Line& addShifted(const Line &other, int shift) {
        assert(other.size() + shift <= size());
        
        uint8_t* target = (uint8_t*)v.data();
        uint8_t* source = (uint8_t*)other.v.data();
        
        /* full bytes are easy */
        target += shift >> 3; shift &= 7;
        
        int bytes = (other.size() >> 3) + 8; // ein reserve uint64_t existiert
        while(bytes >= 15) {
            *(uint64_t*)target     ^= (*(uint64_t*)source) << shift;
            *(uint64_t*)(target+8) ^= (*(uint64_t*)(source+7)) >> (8-shift);
            target += 15, source += 15;
            bytes -= 15;
        }
        /* bytes viele darf ich noch kopieren,
         * die letzten 7 sind optional */
        if (bytes >= 8) {
            *(uint64_t*)target     ^= (*(uint64_t*)source) << shift;
        }
        return *this;
    }
    
    bool scalarProductShifted(const Line &other, int shift) const {
        assert(other.size() + shift <= size());
        
        const uint8_t* target = (uint8_t*)v.data();
        const uint8_t* source = (uint8_t*)other.v.data();
        uint64_t result(0);
        /* full bytes are easy */
        target += shift >> 3; shift &= 7;
        
        int bytes = (other.size() >> 3) + 8; // ein reserve uint64_t existiert
        while(bytes >= 15) {
            result ^= *(uint64_t*)target & (*(uint64_t*)source << shift);
            result ^= *(uint64_t*)(target+8) & (*(uint64_t*)(source+7) >> (8-shift));
            target += 15, source += 15;
            bytes -= 15;
        }
        /* bytes viele darf ich noch kopieren,
         * die letzten 7 sind optional */
        if (bytes >= 8) {
            result ^= *(uint64_t*)target & (*(uint64_t*)source << shift);
        }
        return parityll(result);
    }
    
    
    struct BitsView {
        int byteOffset;
        int shift;
        uint64_t mask;
        BitsView(int firstBit, int numBits) : byteOffset(firstBit >> 3), shift(firstBit & 7), mask((1 << numBits)-1) {}
    };
    
    struct BitView {
        int byteOffset;
        uint8_t mask;
        BitView(int offset) : byteOffset(offset >> 3), mask(1 << (offset & 7)) {}
    };
    
    inline uint64_t get(const BitsView &view) const {
        return ((*((uint64_t*)(((uint8_t*)v.data())+view.byteOffset))) >> view.shift) & view.mask;
    }
    inline void set(const BitsView &view, uint64_t profile) {
        uint64_t &word(*(uint64_t*)(((uint8_t*)v.data())+view.byteOffset));
        word = (word & ~view.mask) | (profile << view.shift);
    }
    inline bool get(const BitView &view) const {
        return (*(((uint8_t*)v.data())+view.byteOffset)) & view.mask;
    }
    inline void set(const BitView &view, bool bit) {
        uint8_t &byte(*(((uint8_t*)v.data())+view.byteOffset));
        byte &= ~view.mask;
        byte |= bit * view.mask;
    }
};

ostream& operator<<(ostream &out, const Line &l) {
    for(int i = 0; i < l.size(); ++i) {
        out << l[i];
    }
    return out;
}

void swap(Line &l1, Line &l2) {
    swap(l1.v,l2.v);
    swap(l1.nBits,l2.nBits);
}

bool scalarProduct(const Line &i1, const Line &i2, int knownZeroes = 0) {
    int fromWord = knownZeroes >> 6;
    int toWord  = (min(i1.nBits, i2.nBits) >> 6)+1;
    assert((int)min(i1.v.size(),i2.v.size()) >= toWord);
    uint64_t prod = 0;
    for(int i = fromWord; i < toWord; ++i) {
        prod ^= i1.v[i] & i2.v[i];
    }
    return parityll(prod);
}

#endif