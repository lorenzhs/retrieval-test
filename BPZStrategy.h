#include "util/util.h"
#include<vector>
#include<stack>
#include<cassert>
#include<array>
#include<iostream>

using namespace std;

template<typename Hashable>
struct BPZStrategy {
	using Hashable_t = Hashable;
    using Hash = DoubleHashSequence<4,Hashable>; /* last one is padding */
	using Solution = vector<bool>;
	static size_t numBits(const Solution& sol) {
		return sol.size();
	}

	uint32_t n;
	Solution sol;
	
	struct Configuration {
		double c;
		friend ostream& operator<<(ostream& out, const Configuration& config) {
			return out << "(c = " << config.c << ")";
		}
	};
	inline static const Configuration Default = { 0.81 };

	typedef array<uint32_t, 3> Edge;
	vector<Edge> edges;
	vector<bool> values;

	BPZStrategy(size_t m, Configuration config) {
		n = m / config.c;
		if (n % 7) n = (n + 8) & ~7;
		edges.reserve(n);
		values.reserve(n);
	}

	void addElement(const Hash& H, bool rhs) {
		edges.emplace_back(Edge{ H[0] % n,H[1] % n,H[2] % n });
		values.emplace_back(rhs);
	}

	bool runConstruction() {
		struct NodeInfo {
            int degree;
            int incidenceXOR;
        };
        vector<NodeInfo> nodes(n,{0,0});
        for(int i = 0; i < (int)edges.size(); ++i) {
            for (int v : edges[i]) {
                nodes[v].degree++;
                nodes[v].incidenceXOR ^= i;
            }
        }
        stack<int> deg1Verts;
        for(int i = 0; i < n; ++i) {
            if (nodes[i].degree == 1) {
                deg1Verts.push(i);
            }
        }
        
        vector<pair<uint32_t,uint32_t>> peelingOrder;
        while(!deg1Verts.empty()) {
            int v = deg1Verts.top(); deg1Verts.pop();
            if (nodes[v].degree) {
                assert(nodes[v].degree == 1);
                int j = nodes[v].incidenceXOR;
                peelingOrder.push_back({j,v});
                Edge &e = edges[j];
                for(int i : e) {
                    nodes[i].incidenceXOR ^= j;
                    nodes[i].degree--;
                    if (nodes[i].degree == 1) {
                        deg1Verts.push(i);
                    }
                }
            }
        }
        int remaining = 0;
        for(NodeInfo &ni : nodes) {
            remaining += ni.degree;
        }
        
        if (remaining) {
            cout << "FAILURE! (Remaining degrees " << remaining << ")" << endl;
			return false;
        } else {
            //cout << "Success!" << endl;
            sol.resize(n);
            for(auto it = peelingOrder.rbegin(); it != peelingOrder.rend(); ++it) {
                Edge &e =   edges[it->first];
                bool val = values[it->first];
                for(uint32_t i : e) {
                    val ^= sol[i];
                }
                sol[it->second] = val;
            }
			return true;
        } 
    }

	inline static bool retrieve(const Solution &sol, const Hash &H) {
		size_t n = sol.size();
		return sol[H[0] % n] ^ sol[H[1] % n] ^ sol[H[2] % n];
	}
};
