#pragma once

#include "util/util.h"
#include<vector>
#include<stack>
#include<cassert>
#include<array>
#include<iostream>

using namespace std;

template<int k,typename Hashable>
struct CoupledStrategy {
	static constexpr char stratName[] = "Coupled";
	using Hashable_t = Hashable;
	using Hash = DoubleHashSequence<k + 1,Hashable>;

	struct Solution {
		vector<bool> v;
		uint32_t n = 0, N = 0;
	};
	static size_t numBits(const Solution& sol) {
		return sol.v.size();
	}

	uint32_t n, N;
	Solution sol;

	struct Configuration {
		double c;
		double eps;
		friend ostream& operator<<(ostream& out, const Configuration& config) {
			return out << "(c = " << config.c << ", z = " << 1/config.eps << ")";
		}
	};

	typedef array<uint32_t, k> Edge;
	vector<Edge> edges;
	vector<bool> values;

	CoupledStrategy(size_t m, Configuration config) {
		n = m / config.c;
		N = config.eps * n;
		edges.reserve(m);
		values.reserve(m);
	}

	void addElement(const Hash& H, bool rhs) {
		uint32_t base = reduce(H[0], n);
		Edge e;
		for (int kk = 0; kk < k; ++kk) {
			e[kk] = base + reduce(H[kk + 1],N);
		}
		edges.push_back(e);
		values.push_back(rhs);
	}

	bool runConstruction() {
		struct NodeInfo {
			int degree;
			int incidenceXOR;
		};
		vector<NodeInfo> nodes(n+N, { 0,0 });
		for (int i = 0; i < (int)edges.size(); ++i) {
			for (int v : edges[i]) {
				nodes[v].degree++;
				nodes[v].incidenceXOR ^= i;
			}
		}
		stack<int> deg1Verts;
		for (int i = 0; i < nodes.size(); ++i) {
			if (nodes[i].degree == 1) {
				deg1Verts.push(i);
			}
		}

		vector<pair<uint32_t, uint32_t>> peelingOrder;
		while (!deg1Verts.empty()) {
			int v = deg1Verts.top(); deg1Verts.pop();
			if (nodes[v].degree) {
				assert(nodes[v].degree == 1);
				int j = nodes[v].incidenceXOR;
				peelingOrder.push_back({ j,v });
				Edge& e = edges[j];
				for (int i : e) {
					nodes[i].incidenceXOR ^= j;
					nodes[i].degree--;
					if (nodes[i].degree == 1) {
						deg1Verts.push(i);
					}
				}
			}
		}
		int remaining = 0;
		for (NodeInfo& ni : nodes) {
			remaining += ni.degree;
		}

		if (remaining) {
			cout << "FAILURE! (Remaining degrees " << remaining << ")" << endl;
			return false;
		}
		else {
			//cout << "Success!" << endl;
			vector<bool>& v(sol.v);
			v.resize(n+N);
			sol.N = N;
			sol.n = n;
			for (auto it = peelingOrder.rbegin(); it != peelingOrder.rend(); ++it) {
				Edge& e = edges[it->first];
				bool val = values[it->first];
				for (uint32_t i : e) {
					val ^= v[i];
				}
				v[it->second] = val;
			}
			return true;
		}
	}

	inline static bool retrieve(const Solution& sol, const Hash& H) {
		/* this heavily assumes little endian and stuff */
		/* this heavily assumes l = 8 and stuff */
		uint32_t base = reduce(H[0],sol.n);
		bool res = false;
		for (int kk = 1; kk <= k; ++kk) {
			res ^= sol.v[base + reduce(H[kk],sol.N)];
		}
		return res;
	}
};
