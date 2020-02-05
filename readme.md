The purpose of this code is to compare the performance of various retrieval data structures.
The following variations have been implemented:

* BPZStrategy.h [1]
* CoupledStrategy.h [2]
* FuseStrategy.h [3]
* GOVStrategy.h [4]
* LMSSStrategy.h [5]
* OneBlockStrategy.h [6]
* TwoBlockStrategy.h [7]

These Strategies are to be used as template parameter to a *retriever*:
* Retriever.h: Lightweight Wrapper, works with [1,2,3,5]
* RetrieverChunked.h: Partition input into chunks, uses Strategy on each chunk. Works with [4,6,7].

ChunkInfo.h contains two Strategies for compressing the meta information of the chunks, if RetrieverChunked is used. RetrieverTester benchmarks a retriever.

I am not much of a coder, so I don't know much about best practices or how to properly package my code for others to use. There are probably portability issues. On windows I use `g++ -O3 -std=c++17 main.cpp`. There is a dependency to boost.

Feel free to anything you want from this folder for your own benchmarks.

* [1]: Botelho, Pagh, Ziviani, 2007, https://doi.org/10.1007/978-3-540-73951-7_13
* [2]: Walzer2020, https://arxiv.org/abs/2001.10500
* [3]: Dietzfelbinger, Walzer 2019, http://dx.doi.org/10.4230/LIPIcs.ESA.2019.38
* [4]: Genuzio, Ottaviano, Vigna, 2016, https://doi.org/10.1007/978-3-319-38851-9_23
* [5]: Luby, Mitzenmacher, Shokrollahi, Spielman, 2001, https://doi.org/10.1109/18.910575
* [6]: Dietzfelbinger, Walzer 2019, http://dx.doi.org/10.4230/LIPIcs.ESA.2019.39
* [7]: Dietzfelbinger, Walzer 2019, http://dx.doi.org/10.4230/LIPIcs.STACS.2019.24
