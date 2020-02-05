#include "util/util.h"
#include<chrono>
#include<iostream>
#include<sstream>

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
        
		Timer queryTimer;
		uint32_t sum = 0;
		for (int i = 0; i < times; ++i) {
			queryTimer.start();
			for (const Hashable& s : keys) {
				sum += R.retrieve(s);
			}
			queryTimer.stop();
		}
        
		uint32_t tHash_mus = hashTimer.median_mus();
		uint32_t tQuery_mus = queryTimer.median_mus();
		double queryPerKeyNS = tQuery_mus * 1000.0 / keys.size();
		double hashPerKeyNS = tHash_mus * 1000.0 / keys.size();
        cout << "QueryCost[ns] " << queryPerKeyNS << endl;
		cout << "HashCost[ns]  " << hashPerKeyNS;
		cout << ((buff - sum) ? " " : "") << endl; /* avoid optimising away */

		log.log("Query Time", queryPerKeyNS);
		log.log("Hash Time", hashPerKeyNS);
    }

	void testAll(const vector<Hashable>& keys, typename Retriever::Configuration config, int times = 1) {
		stringstream ss; ss << config;
		string name = typeid(Retriever).name() + ss.str();
		PerformanceLog& log = PerformanceLog::logs[name];
		cout << "Testing " << name << config << endl;
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

		uint32_t dc = constructionTimer.median_mus();
		double timePerKeyMS = (double)dc / values.size();
		double bitsPerKey = ((double)R.memoryUsage() / values.size());
		log.log("Construction Time per Key", timePerKeyMS);
		log.log("Bits per Key", (double)R.memoryUsage() / values.size());
		cout << "Construction Time ([mus/key] / [mus]) " << timePerKeyMS << " / " << dc << endl;
		cout << "Bits per Key used: " << bitsPerKey << endl;
		verify(R, keys, values);
		testQueries(R, keys, log, times);
	}
};