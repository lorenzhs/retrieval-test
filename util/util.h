#ifndef UTIL_H
#define UTIL_H

#include<random>
#include<vector>
#include<map>
#include<cassert>
#include<unordered_set>
#include<string>
#include<array>
#include<memory>
#include<bitset>
#include<fstream>
#include<chrono>
#include<algorithm>
#include "MurmurHash3.h"
#include "MurmurHash3.cpp"

using namespace std;
using std::chrono::duration_cast;
typedef std::chrono::high_resolution_clock Clock;
typedef std::chrono::milliseconds MilliSeconds;

struct Timer {
private:
	vector<uint32_t> times;
	using Time = Clock::time_point;
	Time currStart;
public:
	Timer() {};
	void start() {
		currStart = Clock::now();
	}
	uint32_t stop() {
		uint32_t d = duration_cast<std::chrono::microseconds>(Clock::now() - currStart).count();
		times.push_back(d);
		return d;
	}
	uint32_t median_mus() {
		int n = times.size();
		sort(times.begin(), times.end());
		if (n % 2) {
			return times[n / 2];
		}
		else {
			return (times[n / 2] + times[(n - 1) / 2]) / 2;
		}
	}
};

template<typename Lambda>
uint32_t benchmark_mus(Lambda &&f) {
	auto start = Clock::now();
	f();
	auto end = Clock::now();
	return duration_cast<std::chrono::microseconds>(end - start).count();
}

template<typename Lambda>
uint32_t benchmark_mus_median(Lambda f, int times) {
	vector<uint32_t> vals(times);
	for (uint32_t& t : vals) {
		t = benchmark_mus(f);
	}
	sort(vals.begin(), vals.end());
	if (times % 2) {
		return vals[times / 2];
	} else {
		return (vals[times / 2] + vals[(times - 1) / 2]) / 2;
	}
}

struct PerformanceLog {
	static map<string, PerformanceLog> logs;
	map<string, vector<double>> data;
	void log(string key, double value) {
		data[key].push_back(value);
	}
};

double avg(const vector<double>& v) {
	if (!v.size()) return NAN;
	double sum = 0.0;
	for (double d : v) sum += d;
	return sum / v.size();
}
ostream& operator<<(ostream& out, const PerformanceLog& log) {
	for (const pair<string, vector<double>>& p : log.data) {
		out << p.first << ": " << avg(p.second) << endl;
	}
	return out;
}

struct BitStream {
	BitStream() {}
	BitStream(vector<uint64_t>&& v) : v(v) {
		this->v.push_back(0ULL); // we need to own another word
	}
	vector<uint64_t> v;
	uint64_t read(uint32_t o, uint8_t l) const {
		uint32_t byte = o >> 3;
		uint8_t  bit = o & 7;
		uint64_t word = *(uint64_t*)((uint8_t*)v.data() + byte);
		return (word >> bit) & ((1ULL << l) - 1);
	}
	size_t extraMemory() const {
		return v.size() * sizeof(uint64_t);
	}
};

struct BitStreamWriter {
	static constexpr uint64_t HALFWORD = (1ULL << 32) - 1;
	uint64_t buff;
	uint8_t pos;
	bool halfWordWritten;
	vector<uint64_t> v;
	BitStreamWriter() : buff(0), pos(0), halfWordWritten(false) {}
	void write(uint32_t x, uint8_t n) {
		buff |= (uint64_t)x << pos;
		pos += n;
		if (pos >= 32) {
			if (halfWordWritten) {
				v.back() |= (buff & HALFWORD) << 32;
			}
			else {
				v.push_back(buff & HALFWORD);
			}
			halfWordWritten = !halfWordWritten;
			buff >>= 32;
			pos -= 32;
		}
	}
	BitStream finalise() {
		if (halfWordWritten) {
			v.back() |= buff << 32;
		}
		else {
			v.push_back(buff);
		}
		return BitStream(std::move(v));
	}
};

struct UniqueRandomInt64 {
	static mt19937_64 rnd;
	static unordered_set<uint64_t> usedRandomInt64;

	uint64_t val;
	UniqueRandomInt64() {
		do {
			val = rnd();
		} while (usedRandomInt64.count(val));
		usedRandomInt64.insert(val);
	}
};

unordered_set<uint64_t> UniqueRandomInt64::usedRandomInt64;
mt19937_64 UniqueRandomInt64::rnd;

void do_hash(const string& s, uint32_t seed, void* buf) {
	MurmurHash3_x64_128(s.c_str(), s.size(), seed, buf);
}


#if USE_MULTIPLY_SHIFT_HASHING

inline uint32_t hash_pair_multiply_shift(uint64_t x, uint64_t a1, uint64_t a2, uint64_t b) {
	return ((a1 + x) * (a2 + (x >> 32)) + b) >> 32;
}

void do_hash(const UniqueRandomInt64& s, uint32_t seed, void* buf) {
	static vector<uint64_t> pair_multiply_shift_data;
	static mt19937_64 gen;
	uint32_t* buf32 = (uint32_t*)buf;
	seed *= 12;
	while (pair_multiply_shift_data.size() < seed + 12) {
		pair_multiply_shift_data.push_back(gen());
	}
	for (int i = 0; i < 4; ++i) {
		const int64_t& a1 = pair_multiply_shift_data[seed++];
		const int64_t& a2 = pair_multiply_shift_data[seed++];
		const int64_t& b = pair_multiply_shift_data[seed++];
		*buf32 = hash_pair_multiply_shift(s.val, a1, a2, b);
		buf32++;
	}
}

#else
void do_hash(const UniqueRandomInt64& s, uint32_t seed, void* buf) {
	MurmurHash3_x64_128(&s, sizeof(UniqueRandomInt64), seed, buf);
}
#endif


//using FourHashes = array<uint24_t, 4>;
template<int k,typename Hashable>
struct DoubleHashSequence {
	inline static const unsigned length = k;
	array<uint32_t, max(4,k)> hashes; // 4 needed to full 128 bits
	uint32_t operator[](uint32_t i) const {
		return hashes[i];
	}
	DoubleHashSequence(const Hashable& s, uint32_t seed) {
		do_hash(s,seed,&hashes);
		if constexpr (k >= 4) {
			uint32_t val = hashes[2];
			for (int i = 4; i < k; ++i) {
				hashes[i] = (val += hashes[3]);
			}
		}
	}
};

template<typename Hashable>
struct DynamicDoubleHashSequence {
private:
	array<uint32_t, 4> buff;
	static vector<uint32_t> vec;
public:
	/* Hashes will only be valid until someone else calls getHashes
	 * very much not thread-safe */
	const vector<uint32_t>& getMoreHashes(int d) const {
		vec.resize(max(d,3));
		vec[0] = buff[1];
		vec[1] = buff[2];
		vec[2] = buff[3];
		uint32_t val = buff[2];
		for (int i = 3; i < d; ++i) {
			vec[i] = (val += buff[3]);
		}
		return vec;
	}
	uint32_t getFirstHash() const {
		return buff[0];
	}
	DynamicDoubleHashSequence(const Hashable& s, uint32_t seed) {
		do_hash(s, seed, &buff);
	}
	/* should not ordinarily be used as it computes all hashes */
	uint32_t operator[](uint32_t i) const {
		return getMoreHashes(i+1)[i];
	}
};

template<typename T>
vector<uint32_t> DynamicDoubleHashSequence<T>::vec;

struct KeyFile {
private: vector<string> strings;
public:
	const vector<string>& keys() const {
		return strings;
	}
	KeyFile(const string& fname, int limit = 0) {
		std::ios::sync_with_stdio(false);
		strings.reserve(limit);
		std::ifstream file(fname);
		string buf;
		while (std::getline(file, buf) && limit--)
		{
			strings.push_back(std::move(buf));
		}
	}
};

/*
struct FileIterator {
    string command;
    int limit;
private:
    FileIterator(const string &fname, const string &c, int limit = 0) : limit(limit) {
        command = c+fname + " 2>&1";
    }
public:
    static FileIterator zippedFileIterator(const string &fname, int limit = 0) {
        return FileIterator(fname, "7z e -so ", limit);
    }
    static FileIterator textFileIterator(const string &fname, int limit = 0) {
        return FileIterator(fname, "cat ", limit);
    }
    template<typename Callable>
    void foreachLine(Callable f) const {
        std::array<char, 2050> buffer; // längster String hat Länge 2047
        buffer[sizeof(buffer)-2] = '\0';
        std::string line;
        std::shared_ptr<FILE> pipe(_popen(command.c_str(), "r"), _pclose);
        if (!pipe) throw std::runtime_error("popen() failed!");
        int count = limit;
        while (!feof(pipe.get()) && (!limit || count--)) {
            if (fgets(buffer.data(), sizeof(buffer), pipe.get()) != nullptr) {
                f(buffer.data());
                assert(buffer[sizeof(buffer)-2] == '\0'); // should be long enough
            }
        }
        assert(count == -1 || count == 0);
    }
};*/

template<int l>
inline unsigned extract(const void *array, int blockID) {
    static_assert(l <= sizeof(uint32_t)*8-8); /*not wider than 3 bytes */
    return (*(uint32_t*)((uint8_t*)array+l*blockID/8) >> (l*blockID % 8)) & ((1 << l)-1);
}
template<> inline unsigned extract<8>(const void *array, int blockID) {
    return ((uint8_t*)array)[blockID];
}
template<> inline unsigned extract<16>(const void *array, int blockID) {
    return ((uint16_t*)array)[blockID];
}
template<> inline unsigned extract<32>(const void *array, int blockID) {
    return ((uint32_t*)array)[blockID];
}

template<typename IntType>
bool scalarProduct(IntType i1, IntType i2) = delete;

#ifdef __GNUC__
	inline bool parityll(uint64_t x) { return __builtin_parityll(x); }
	inline bool parityl(uint32_t x) { return __builtin_parityl(x); }
	template<> inline bool scalarProduct<unsigned long     >(unsigned long      i1, unsigned long      i2) { return __builtin_parityl(i1 & i2); }
	template<> inline bool scalarProduct<unsigned char     >(unsigned char      i1, unsigned char      i2) { return __builtin_parity(i1 & i2); }
	template<> inline bool scalarProduct<unsigned long long>(unsigned long long i1, unsigned long long i2) { return __builtin_parityll(i1 & i2); }
	
	inline uint32_t ctz(uint32_t x) { return __builtin_ctz(x); }
	inline uint32_t clz(uint32_t x) { return __builtin_clz(x); }
	inline uint32_t ctzll(uint64_t x) { return __builtin_ctzll(x); }

	template<int l>
	int popcount(const bitset<l>& s) {
		static_assert(l <= 32);
		return __builtin_popcount(s.to_ulong());
	}
#else
	#include <intrin.h>
	inline bool parityll(uint64_t x) { return __popcnt64(x) & 1; }
	inline bool parityl(uint32_t x) { return __popcnt(x) & 1; }
	template<> inline bool scalarProduct<unsigned char     >(unsigned char      i1, unsigned char      i2) { return __popcnt16(i1 & i2) & 1; }
	template<> inline bool scalarProduct<unsigned long     >(unsigned long      i1, unsigned long      i2) { return __popcnt  (i1 & i2) & 1; }
	template<> inline bool scalarProduct<unsigned long long>(unsigned long long i1, unsigned long long i2) { return __popcnt64(i1 & i2) & 1; }

	#pragma intrinsic(_BitScanForward)
	inline uint32_t ctz(uint32_t x) { unsigned long ret;  _BitScanForward(&ret, x); return ret; }
	inline uint32_t clz(uint32_t x) { unsigned long ret;  _BitScanReverse(&ret, x); return 31-ret; }
	inline uint32_t ctzll(uint64_t x) { unsigned long ret;  _BitScanForward64(&ret, x); return ret; }

	template<int l>
	int popcount(const bitset<l>& s) {
		static_assert(l <= 32);
		return __popcnt(s.to_ulong());
	}
#endif

#endif