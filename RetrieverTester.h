#include "util/util.h"
#include<chrono>
#include<iostream>
#include<sstream>

#ifdef __GNUC__
	#include <cxxabi.h>
#endif

using std::chrono::duration_cast;
typedef std::chrono::high_resolution_clock Clock;
typedef std::chrono::milliseconds MilliSeconds;

/* Todo: Keys other than string
		 More than one bit retrieval */
template<typename Retriever>
struct RetrieverTester {
	using Hashable = typename Retriever::Hashable;
	RetrieverTester() {}
	using Hash = typename Retriever::Hash;
	
	static bool toRetrieve(const string& s) { return s.length() & 1; }
	static bool toRetrieve(const UniqueRandomInt64& s) { return s.val & 1; }

	bool verify(const Retriever &R, const vector<Hashable> &keys, const vector<bool> &values) {
        cout << "Verifying ... ";
        int errors = 0;
		assert(keys.size() == values.size());
		for(int i = 0; i < keys.size(); ++i)
            if (R.retrieve(keys[i]) != values[i])
                errors++;
        errors ? (cout << errors << " errors found!") : cout << "Done", cout << endl;
        return errors != 0;
    }
    
    void testQueries(const Retriever& R, const vector<Hashable>& keys, PerformanceLog &log, int times = 1) {
        uint32_t buff = 0;
		Timer hashTimer;
		for (int i = 0; i < times; ++i) {
			hashTimer.start();
			for (const Hashable& s : keys) {
				Hash H(s, 0);
				buff ^= *(uint32_t*)& H;
			}
			hashTimer.stop();
		}
        
		/* in cache queries */
		Timer queryCachedTimer;
		uint32_t sum = 0;
		{
			for (int i = 0; i < times; ++i) {
				queryCachedTimer.start();
				// introduce data dependency between queries
				for(int j = 0; j < 2; ++j)
				for (int i = 0; i < keys.size(); i+=2) {
					sum += R.retrieve(keys[i ^ (sum & 1)]);
				}
				/*
				for (const Hashable& s : keys) {
					sum += R.retrieve(s);
				}*/
				queryCachedTimer.stop();
			}
		}		
        
		/* out of cache tests */
		Timer queryUncachedTimer;
		{
			constexpr int nCopies = 1 << 10;
			vector<Retriever> copies(nCopies, R);

			for (int i = 0; i < times; ++i) {
				queryUncachedTimer.start();
				for (int i = 0; i < keys.size(); ++i) {
					sum += copies[(sum+i) & (nCopies - 1)].retrieve(keys[i]);
				}
				queryUncachedTimer.stop();
			}
		}


		uint32_t tHash_mus = hashTimer.median_mus();
		uint32_t tQueryCached_mus = queryCachedTimer.median_mus();
		uint32_t tQueryUncached_mus = queryUncachedTimer.median_mus();
		double queryCachedNS = tQueryCached_mus * 1000.0 / keys.size();
		double queryUncachedNS = tQueryUncached_mus * 1000.0 / keys.size();
		double hashPerKeyNS = tHash_mus * 1000.0 / keys.size();
        cout << "QueryCost[ns] " << queryCachedNS << "/" << queryUncachedNS <<  endl;
		cout << "HashCost[ns]  " << hashPerKeyNS;
		cout << ((buff - sum) ? " " : "") << endl; /* avoid optimising away */

		log.log("Query Time Cached", queryCachedNS);
		log.log("Query Time Uncached", queryUncachedNS);
		log.log("Hash Time", hashPerKeyNS);
    }

	void replaceAll(string &name, const string& search, const string& replace) {
		size_t pos = 0;
		while ((pos = name.find(search, pos)) != std::string::npos) {
			name.replace(pos, search.length(), replace);
			pos += replace.length();
		}
	}
	void prettifyName(std::string &name) {
		replaceAll(name,"struct ", "");
		replaceAll(name, "UniqueRandomInt64", "Int64");
		replaceAll(name, "Strategy", "");
		replaceAll(name, "ChunkInfo", "");
	}

	void testAll(const vector<Hashable>& keys, typename Retriever::Configuration config, int times = 1, string texConfig = "") {
		stringstream ss;
		#ifdef __GNUC__
			int status;
			char *realname = abi::__cxa_demangle(typeid(Retriever).name(), 0, 0, &status);
			ss << realname << config;
			free(realname);
		#else
			ss << typeid(Retriever).name() << config;
		#endif
		string name = ss.str();
		prettifyName(name);
		if (!PerformanceLog::logs.count(name)) PerformanceLog::creationOrder.push_back(name);
		PerformanceLog& log = PerformanceLog::logs[name];
		log.setShortName(Retriever::Strategy::stratName);
		log.setTexConfig(texConfig);
		cout << "Testing " << name << endl;
		vector<bool> values; values.reserve(keys.size());
		for(const Hashable&s : keys) {
			values.push_back(toRetrieve(s));
		}

		Timer constructionTimer;
		for (int i = 0; i < times - 1; ++i) {
			constructionTimer.start();
			Retriever R(keys, values, config);
			constructionTimer.stop();
		}
		
		constructionTimer.start();
		Retriever R(keys, values, config);
		constructionTimer.stop();

		if (!R.hasSucceeded()) {
			cout << "GIVING UP" << endl;
			return;
		}

		uint32_t dc = constructionTimer.median_mus();
		double timePerKeyMS = (double)dc / values.size();
		double bitsPerKey = ((double)R.memoryUsage() / values.size());
		log.log("Construction Time per Key", timePerKeyMS);
		log.log("Bits per Key", (double)R.memoryUsage() / values.size());
		log.log("Average Seed", R.getAvgSeed());
		cout << "Construction Time ([mus/key] / [mus]) " << timePerKeyMS << " / " << dc << endl;
		cout << "Bits per Key used: " << bitsPerKey << endl;
		verify(R, keys, values);
		testQueries(R, keys, log, max(5,times));
	}
};