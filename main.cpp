#define USE_MULTIPLY_SHIFT_HASHING 1

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

map<string, PerformanceLog> PerformanceLog::logs;

int main() {
	int n = 1E7;

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
	
	for (int i = 0; i < outerRepeats; ++i) {
		// START
		{
			using Strat = GOVStrategy<3, Hashable>;
			using Ret   = RetrieverChunked<Strat, ChunkInfoPacked>; // ChunkInfoCompressed ?
			Strat::Configuration sConfig{ 0.91 };
			Ret::Configuration config{ sConfig, 10000 };
			RetrieverTester<Ret> tester;
			tester.testAll(keys, config, innerRepeats);
		}
		{
			using Strat = GOVStrategy<4, Hashable>;
			using Ret   = RetrieverChunked<Strat, ChunkInfoPacked>; // ChunkInfoCompressed
			Strat::Configuration sConfig{ 0.97 };
			Ret::Configuration config{ sConfig, 10000 };
			RetrieverTester<Ret> tester;
			tester.testAll(keys, config, innerRepeats);
		}
		{
			using Strat = TwoBlockStrategy<16, Hashable>;
			using Ret   = RetrieverChunked<Strat, ChunkInfoCompressed>;
			Strat::Configuration sConfig{ 1 / 1.0005 };
			Ret::Configuration config{ sConfig, 10000 };
			RetrieverTester<Ret> tester;
			tester.testAll(keys, config, innerRepeats);
		}
		{
			// less dense!
			using Strat = OneBlockStrategy<Hashable>;
			using Ret   = RetrieverChunked<Strat>;
			Strat::Configuration sConfig{ 0.93 };
			Ret::Configuration config{ sConfig, 10000 };
			RetrieverTester<Ret> tester;
			tester.testAll(keys, config, innerRepeats);
		}
		{
			using Strat = OneBlockStrategy<Hashable>;
			using Ret   = RetrieverChunked<Strat>;
			Strat::Configuration sConfig{ 0.95 };
			Ret::Configuration config{ sConfig, 10000 };
			RetrieverTester<Ret> tester;
			tester.testAll(keys, config, innerRepeats);
		}
		{
			// denser!
			using Strat = OneBlockStrategy<Hashable>;
			using Ret   = RetrieverChunked<Strat>;
			Strat::Configuration sConfig{ 0.97 };
			Ret::Configuration config{ sConfig, 10000 };
			RetrieverTester<Ret> tester;
			tester.testAll(keys, config, innerRepeats);
		}
		// Unchunked strategies
		{
			using Strat = BPZStrategy<Hashable>;
			Strat::Configuration config{ 0.81 };
			RetrieverTester<Retriever<Strat>> tester;
			tester.testAll(keys, config, innerRepeats);
		}
		{
			using Strat = RetrieverLMSS<Hashable>;
			Strat::Configuration config{ 0.9, 12 };
			RetrieverTester<Retriever<Strat>> tester;
			tester.testAll(keys, config, innerRepeats);
		}
		{
			using Strat = RetrieverLMSS<Hashable>;
			Strat::Configuration config{ 0.99, 150 };
			RetrieverTester<Retriever<Strat>> tester;
			tester.testAll(keys, config, innerRepeats);
		}
		{
			using Coupled3 = CoupledStrategy<3, Hashable>;
			Coupled3::Configuration config{ 0.91, 0.01 * (3 - 1) };
			RetrieverTester<Retriever<Coupled3>> tester;
			tester.testAll(keys, config, innerRepeats);
		}
		{
			using Fuse3 = FuseStrategy<3, Hashable>;
			Fuse3::Configuration config{ 0.91,100 };
			RetrieverTester<Retriever<Fuse3>> tester;
			tester.testAll(keys, config, innerRepeats);
		}
		{
			using Coupled4 = CoupledStrategy<4, Hashable>;
			Coupled4::Configuration config{ 0.96,0.005 * (4 - 1) };
			RetrieverTester<Retriever<Coupled4>> tester;
			tester.testAll(keys, config, innerRepeats);
		}
		{
			using Fuse4 = FuseStrategy<4, Hashable>;
			Fuse4::Configuration config{ 0.96,200 };
			RetrieverTester<Retriever<Fuse4>> tester;
			tester.testAll(keys, config, innerRepeats);
		}
		{
			using Coupled7 = CoupledStrategy<7, Hashable>;
			Coupled7::Configuration config{ 0.9849,0.002 * (7 - 1) * 1.15 }; // needs a bit more slack, i.e. eps' = eps*(1.5) or so.
			RetrieverTester<Retriever<Coupled7>> tester;
			tester.testAll(keys, config, innerRepeats);
		}
		{
			using Fuse7 = FuseStrategy<7, Hashable>;
			Fuse7::Configuration config{ 0.9849,500 };
			RetrieverTester<Retriever<Fuse7>> tester;
			tester.testAll(keys, config, innerRepeats);
		}
	}

	for (const auto& namedLog : PerformanceLog::logs) {
		cout << "###" << namedLog.first << "###" << endl;
		cout << namedLog.second << endl;
	}
	getchar();
}