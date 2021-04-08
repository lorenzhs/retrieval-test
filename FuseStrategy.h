#pragma once

#include "util/util.h"
#include<vector>
#include<stack>
#include<cassert>
#include<array>
#include<iostream>

using namespace std;

namespace walzer {

template<int k,typename Hashable>
struct FuseStrategy {
	static constexpr char stratName[] = "Fuse";
	using Hashable_t = Hashable;
	using Hash = DoubleHashSequence<k+1,Hashable>;

	struct Solution {
		vector<bool> v;
		uint32_t span = 0, l = 0;
	};
	static size_t numBits(const Solution& sol) {
		return sol.v.size();
	}

	uint32_t n,span,l;
	Solution sol;

	struct Configuration {
		double c;
		uint32_t l;
		friend ostream& operator<<(ostream& out, const Configuration& config) {
			return out << "(c = " << config.c << ", l = " << config.l << ")";
		}
	};

	typedef array<uint32_t, k> Edge;
	vector<Edge> edges;
	vector<bool> values;

	FuseStrategy(size_t m, Configuration config) {
		l = config.l;
		n = m / config.c;
		if (n % l) n += l - n % l;
		span = n / l;
		n += (k - 1) * span;
		edges.reserve(m);
		values.reserve(m);
	}

	void addElement(const Hash& H, bool rhs) {
		uint32_t base = reduce(H[0], l) * span;
		Edge e;
		for (int kk = 0; kk < k; ++kk, base += span) {
			e[kk] = base + reduce(H[kk+1],span);
		}
		edges.push_back(e);
		values.push_back(rhs);
	}

	bool runConstruction() {
		struct NodeInfo {
			int degree;
			int incidenceXOR;
		};
		vector<NodeInfo> nodes(n, { 0,0 });
		for (int i = 0; i < (int)edges.size(); ++i) {
			for (int v : edges[i]) {
				nodes[v].degree++;
				nodes[v].incidenceXOR ^= i;
			}
		}
		stack<int> deg1Verts;
		for (int i = 0; i < n; ++i) {
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
				Edge & e = edges[j];
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
			v.resize(n);
			sol.span = span;
			sol.l = l;
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

	inline static bool retrieve(const Solution & sol, const Hash & H) {
		/* this heavily assumes little endian and stuff */
		/* this heavily assumes l = 8 and stuff */
		uint32_t base = reduce(H[0], sol.l) * sol.span;
		bool res = false;
		for (int kk = 1; kk <= k; ++kk, base += sol.span) {
			res ^= sol.v[base + reduce(H[kk],sol.span)];
		}
		return res;
	}
};
}
