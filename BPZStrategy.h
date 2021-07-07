#include "util/util.h"

#include <array>
#include <cassert>
#include <iostream>
#include <stack>
#include <vector>

// not faster in my experiments, maybe smaller memory footprint, unclear.
// #define ONE_STACK_PEELING

using namespace std;

namespace walzer {
template <typename Hashable, typename data_t = bool>
struct BPZStrategy {
    static constexpr char stratName[] = "BPZ";
    using Hashable_t = Hashable;
    using Hash = DoubleHashSequence<4, Hashable>; /* last one is padding */
    using Solution = vector<data_t>;
    static size_t numBits(const Solution& sol) {
        return sol.size() * (std::is_same_v<data_t, bool> ? 1 : sizeof(data_t) * 8);
    }

    uint32_t n;
    Solution sol;

    struct Configuration {
        double c;
        friend ostream& operator<<(ostream& out, const Configuration& config) {
            return out << "(c = " << config.c << ")";
        }
    };
    inline static const Configuration Default = {0.81};

    typedef array<uint32_t, 3> Edge;
    vector<Edge> edges;
    vector<data_t> values;

    BPZStrategy(size_t m, Configuration config) {
        n = m / config.c;
        if (n & 7)
            n = (n + 8) & ~7; // extend to byte
        edges.reserve(n);
        values.reserve(n);
    }

    void addElement(const Hash& H, data_t rhs) {
        edges.emplace_back(Edge{reduce(H[0], n), reduce(H[1], n), reduce(H[2], n)});
        values.emplace_back(rhs);
    }

    bool runConstruction() {
        struct NodeInfo {
            int degree;
            int incidenceXOR;
        };
        vector<NodeInfo> nodes(n, {0, 0});
        for (int i = 0; i < (int)edges.size(); ++i) {
            for (int v : edges[i]) {
                nodes[v].degree++;
                nodes[v].incidenceXOR ^= i;
            }
        }
#ifdef ONE_STACK_PEELING
        vector<pair<uint32_t, uint32_t>> peelingOrder;
        peelingOrder.reserve(1.1 * edges.size());

        for (int v = 0; v < n; ++v) {
            if (nodes[v].degree == 1) {
                peelingOrder.push_back({nodes[v].incidenceXOR, v});
            }
        }

        int done = 0, scanned = 0;
        while (scanned < peelingOrder.size()) {
            auto p = peelingOrder[scanned++];
            uint32_t v = p.second;
            uint32_t j = p.first;
            if (nodes[p.second].degree) {
                peelingOrder[done++] = p;
                assert(nodes[v].degree == 1);
                assert(nodes[v].incidenceXOR == j);
                Edge& e = edges[j];
                for (int i : e) {
                    nodes[i].incidenceXOR ^= j;
                    nodes[i].degree--;
                    if (nodes[i].degree == 1) {
                        peelingOrder.push_back({nodes[i].incidenceXOR, i});
                    }
                }
            }
        }
        peelingOrder.resize(done);
#else
        stack<int> deg1Verts;
        for (uint32_t i = 0; i < n; ++i) {
            if (nodes[i].degree == 1) {
                deg1Verts.push(i);
            }
        }

        vector<pair<uint32_t, uint32_t>> peelingOrder;
        peelingOrder.reserve(edges.size());
        while (!deg1Verts.empty()) {
            int v = deg1Verts.top();
            deg1Verts.pop();
            if (nodes[v].degree) {
                assert(nodes[v].degree == 1);
                int j = nodes[v].incidenceXOR;
                peelingOrder.push_back({j, v});
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
#endif
        int remaining = 0;
        for (NodeInfo& ni : nodes) {
            remaining += ni.degree;
        }

        if (remaining) {
            cout << "FAILURE! (Remaining degrees " << remaining << ")" << endl;
            return false;
        } else {
            // cout << "Success!" << endl;
            sol.resize(n);
            for (auto it = peelingOrder.rbegin(); it != peelingOrder.rend(); ++it) {
                Edge& e = edges[it->first];
                data_t val = values[it->first];
                for (uint32_t i : e) {
                    val ^= sol[i];
                }
                sol[it->second] = val;
            }
            return true;
        }
    }

    inline static data_t retrieve(const Solution& sol, const Hash& H) {
        size_t n = sol.size();
        return sol[reduce(H[0], n)] ^ sol[reduce(H[1], n)] ^ sol[reduce(H[2], n)];
    }
};
} // namespace walzer
