//#define USE_MULTIPLY_SHIFT_HASHING 1
#define USE_XXHASH 1

#include "util/util.h"
#include "RetrieverTester.h"
#include "Retriever.h"
#include "RetrieverChunked.h"
#include "LMSSStrategy.h"
#include "GOVStrategy.h"
#include "CoupledStrategy.h"
#include "FuseStrategy.h"
#include "OneBlockStrategy.h"
#include "TwoBlockStrategy.h"
#include "BPZStrategy.h"
#include <cstring>

using namespace walzer;

map<string, PerformanceLog> PerformanceLog::logs;
vector<string> PerformanceLog::creationOrder;

int main() {
	size_t n = 1E7;

	/* Tests with strings also implemented:
	KeyFile file("../data/eu-2015-host.ids", n);
	const vector<string>& keys = file.keys();
	using Hashable = string;*/

	vector<UniqueRandomInt64> keys(n);
	using Hashable = UniqueRandomInt64;

	/* Interleaving the tests of the different strategies,
	 * because apparently my CPU gets slower in the later trials,
	 * maybe because it gets hot. */
	int outerRepeats = 10;
	int innerRepeats = 1;

	#define TEST_BPZ
	#define TEST_GOV3
	#define TEST_GOV4
	#define TEST_LMSS
	#define TEST_TWOBLOCK
	#define TEST_ONEBLOCK
	#define TEST_COUPLED

	//#define TEST_FUSE

	char teXParamsBuf[100];

	for (int i = 0; i < outerRepeats; ++i) {
		cout << "#####################################" << endl;
		cout << "# TESTING ROUND " << (i+1) << "/" << outerRepeats << endl;
		cout << "#####################################" << endl;
		// strategies by others
		#ifdef TEST_BPZ
		{
			double c = 0.81;
			sprintf(teXParamsBuf, "c = %.2f", c);
			using Strat = BPZStrategy<Hashable>;
			Strat::Configuration config{ c };
			RetrieverTester<Retriever<Strat>> tester;
			tester.testAll(keys, config, innerRepeats, teXParamsBuf);
		}
		#endif
		#ifdef TEST_GOV3
		{
			constexpr int k = 3; double c = 0.91; uint32_t C = 10000;
			sprintf(teXParamsBuf, "k = %d, c = %.2f, C = %d", k, c, C);
			using Strat = GOVStrategy<k, Hashable>;
			using Ret   = RetrieverChunked<Strat, ChunkInfoPacked>;
			Strat::Configuration sConfig{ c };
			Ret::Configuration config{ sConfig, C };
			RetrieverTester<Ret> tester;
			tester.testAll(keys, config, innerRepeats, teXParamsBuf);
		}
		{
			constexpr int k = 3; double c = 0.90; uint32_t C = 1000;
			sprintf(teXParamsBuf, "k = %d, c = %.2f, C = %d", k, c, C);
			using Strat = GOVStrategy<k, Hashable>;
			using Ret = RetrieverChunked<Strat, ChunkInfoCompressed>;
			Strat::Configuration sConfig{ c };
			Ret::Configuration config{ sConfig, C };
			RetrieverTester<Ret> tester;
			tester.testAll(keys, config, innerRepeats, teXParamsBuf);
		}
		#endif
		#ifdef TEST_GOV4
		{
			constexpr int k = 4; double c = 0.97; uint32_t C = 10000;
			sprintf(teXParamsBuf, "k = %d, c = %.2f, C = %d", k, c, C);
			using Strat = GOVStrategy<k, Hashable>;
			using Ret   = RetrieverChunked<Strat, ChunkInfoPacked>;
			Strat::Configuration sConfig{ c };
			Ret::Configuration config{ sConfig, C };
			RetrieverTester<Ret> tester;
			tester.testAll(keys, config, innerRepeats, teXParamsBuf);
		}
		{
			constexpr int k = 4; double c = 0.96; uint32_t C = 1000;
			sprintf(teXParamsBuf, "k = %d, c = %.2f, C = %d", k, c, C);
			using Strat = GOVStrategy<k, Hashable>;
			using Ret = RetrieverChunked<Strat, ChunkInfoCompressed>;
			Strat::Configuration sConfig{ c };
			Ret::Configuration config{ sConfig, C };
			RetrieverTester<Ret> tester;
			tester.testAll(keys, config, innerRepeats, teXParamsBuf);
		}
		#endif
		#ifdef TEST_LMSS
		{
			constexpr int D = 12; double c = 0.9;
			sprintf(teXParamsBuf, "D = %d, c = %.2f", D, c);
			using Strat = RetrieverLMSS<Hashable>;
			Strat::Configuration config{ c, D };
			RetrieverTester<Retriever<Strat>> tester;
			tester.testAll(keys, config, innerRepeats, teXParamsBuf);
		}
		{
			constexpr int D = 150; double c = 0.99;
			sprintf(teXParamsBuf, "D = %d, c = %.2f", D, c);
			using Strat = RetrieverLMSS<Hashable>;
			Strat::Configuration config{ c, D };
			RetrieverTester<Retriever<Strat>> tester;
			tester.testAll(keys, config, innerRepeats, teXParamsBuf);
		}
		{
			GAMED_SEED = 2;
			constexpr int D = 800; double c = 0.997;
			sprintf(teXParamsBuf, "D = %d, c = %.3f", D, c);
			using Strat = RetrieverLMSS<Hashable>;
			Strat::Configuration config{ c, D };
			RetrieverTester<Retriever<Strat>> tester;
			tester.testAll(keys, config, innerRepeats, teXParamsBuf);
			GAMED_SEED = 0;
		}
		#endif
		#ifdef TEST_COUPLED
		{
			constexpr int k = 3; double z = 120; double c = 0.905;
			sprintf(teXParamsBuf, "k = %d, z = %.0f, c = %.2f", k, z, c);
			using Coupled3 = CoupledStrategy<k, Hashable>;
			Coupled3::Configuration config{ c, 1 / z };
			RetrieverTester<Retriever<Coupled3>> tester;
			tester.testAll(keys, config, innerRepeats, teXParamsBuf);
		}
		{
			constexpr int k = 4; double z = 120; double c = 0.96;
			sprintf(teXParamsBuf, "k = %d, z = %.0f, c = %.2f", k, z, c);
			using Coupled4 = CoupledStrategy<k, Hashable>;
			Coupled4::Configuration config{ c, 1 / z };
			RetrieverTester<Retriever<Coupled4>> tester;
			tester.testAll(keys, config, innerRepeats, teXParamsBuf);
		}
		{
			constexpr int k = 7; double z = 120; double c = 0.979; // needs a bit more slack, i.e. eps' = eps*(1.5) or so.
			sprintf(teXParamsBuf, "k = %d, z = %.0f, c = %.2f", k, z, c);
			using Coupled7 = CoupledStrategy<k, Hashable>;
			Coupled7::Configuration config{ c, 1 / z };
			RetrieverTester<Retriever<Coupled7>> tester;
			tester.testAll(keys, config, innerRepeats, teXParamsBuf);
		}
		#endif
		#ifdef TEST_FUSE
		{
			using Fuse3 = FuseStrategy<3, Hashable>;
			Fuse3::Configuration config{ 0.91,100 };
			RetrieverTester<Retriever<Fuse3>> tester;
			tester.testAll(keys, config, innerRepeats);
		}
		{
			using Fuse4 = FuseStrategy<4, Hashable>;
			Fuse4::Configuration config{ 0.96,200 };
			RetrieverTester<Retriever<Fuse4>> tester;
			tester.testAll(keys, config, innerRepeats);
		}
		{
			using Fuse7 = FuseStrategy<7, Hashable>;
			Fuse7::Configuration config{ 0.9849,500 };
			RetrieverTester<Retriever<Fuse7>> tester;
			tester.testAll(keys, config, innerRepeats);
		}
		#endif
		#ifdef TEST_TWOBLOCK
		{
			constexpr int l = 16; double c = 1 / 1.0005; uint32_t C = 10000;
			sprintf(teXParamsBuf, "\\ell = %d, c = %.4f, C = %d", l, c, C);
			using Strat = TwoBlockStrategy<l, Hashable>;
			using Ret   = RetrieverChunked<Strat, ChunkInfoCompressed>;
			Strat::Configuration sConfig{ c };
			Ret::Configuration config{ sConfig, C };
			RetrieverTester<Ret> tester;
			tester.testAll(keys, config, innerRepeats, teXParamsBuf);
		}
		{
			constexpr int l = 16; double c = 1 / 1.0005; uint32_t C = 20000;
			sprintf(teXParamsBuf, "\\ell = %d, c = %.4f, C = %d", l, c, C);
			using Strat = TwoBlockStrategy<l, Hashable>;
			using Ret = RetrieverChunked<Strat, ChunkInfoCompressed>;
			Strat::Configuration sConfig{ c };
			Ret::Configuration config{ sConfig, C };
			RetrieverTester<Ret> tester;
			tester.testAll(keys, config, innerRepeats, teXParamsBuf);
		}
		#endif

		#ifdef TEST_ONEBLOCK
		{
			// working comfortably
			constexpr int L = 64; double c = 0.93; uint32_t C = 10000;
			sprintf(teXParamsBuf, "L = %d, c = %.2f, C = %d", L, c, C);
			using Strat = OneBlockStrategy<Hashable,L/64>;
			using Ret   = RetrieverChunked<Strat>;
			Strat::Configuration sConfig{ c };
			Ret::Configuration config{ sConfig, C };
			RetrieverTester<Ret> tester;
			tester.testAll(keys, config, innerRepeats, teXParamsBuf);
		}
		{
			// struggling, large average seeds
			constexpr int L = 64; double c = 0.96; uint32_t C = 10000;
			sprintf(teXParamsBuf, "L = %d, c = %.2f, C = %d", L, c, C);
			using Strat = OneBlockStrategy<Hashable, L/64>;
			using Ret = RetrieverChunked<Strat>;
			Strat::Configuration sConfig{ c };
			Ret::Configuration config{ sConfig, C };
			RetrieverTester<Ret> tester;
			tester.testAll(keys, config, innerRepeats, teXParamsBuf);
		}
		{
			constexpr int L = 128; double c = 0.98; uint32_t C = 20000;
			sprintf(teXParamsBuf, "L = %d, c = %.2f, C = %d", L, c, C);
			using Strat = OneBlockStrategy<Hashable,L/64>;
			using Ret   = RetrieverChunked<Strat>;
			Strat::Configuration sConfig{ c };
			Ret::Configuration config{ sConfig, C };
			RetrieverTester<Ret> tester;
			tester.testAll(keys, config, innerRepeats, teXParamsBuf);
		}
		{
			constexpr int L = 3*64; double c = 0.99; uint32_t C = 50000;
			sprintf(teXParamsBuf, "L = %d, c = %.2f, C = %d", L, c, C);
			using Strat = OneBlockStrategy<Hashable, L/64>;
			using Ret = RetrieverChunked<Strat>;
			Strat::Configuration sConfig{ c };
			Ret::Configuration config{ sConfig, C };
			RetrieverTester<Ret> tester;
			tester.testAll(keys, config, innerRepeats, teXParamsBuf);
		}
		{
			GAMED_SEED = 1;
			constexpr int L = 16 * 64; double c = 0.996;
			uint32_t C = n;
			sprintf(teXParamsBuf, "L = %d, c = %.4f, C = n", L, c);
			using Strat = OneBlockStrategy<Hashable, L/64>;
			using Ret = RetrieverChunked<Strat>;
			Strat::Configuration sConfig{ c };
			Ret::Configuration config{ sConfig, (uint32_t)C }; // single chunk
			RetrieverTester<Ret> tester;
			tester.testAll(keys, config, innerRepeats, teXParamsBuf);
			GAMED_SEED = 0;
		}
		#endif
	}

	for (const auto& namedLog : PerformanceLog::logs) {
		cout << "###" << namedLog.first << "###" << endl;
		cout << namedLog.second << endl;
	}
	cout << "\\def\\retrieverData{%" << endl;
	bool first = true;
	for (string& name : PerformanceLog::creationOrder) {
		PerformanceLog& log = PerformanceLog::logs[name];
		if (!first) cout << ",%" << endl; first = false;
		cout << "\\" << log.shortName
			 << "/" << "{" << log.texConfig << "}"
			 << "/" << avg(log.data["Construction Time per Key"])
			 << "/" << avg(log.data["Bits per Key"])
			 << "/" << avg(log.data["Query Time Cached"])
			 << "/" << avg(log.data["Query Time Uncached"])
			// << "/" << "{" << namedLog.first << "}"
			;
	}
	cout << "}" << endl;
	getchar();
}
