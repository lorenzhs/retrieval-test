#include "util/util.h"

#include <array>
#include <cassert>
#include <iostream>
#include <stack>
#include <stdexcept>
#include <vector>

using namespace std;

namespace walzer {

template <typename Hashable>
struct RetrieverLMSS {
    static constexpr char stratName[] = "LMSS";
    using Hashable_t = Hashable;
    using Hash = DynamicDoubleHashSequence<Hashable>;

    vector<uint32_t> dist;
    void buildDistribution(int D) {
        uint32_t MAX = std::numeric_limits<uint32_t>::max();
        vector<double> lamb(D + 2);
        double H = 0.0;
        for (int i = 1; i <= D; H += 1.0 / i++)
            ;
        for (int i = 2; i <= D + 1; lamb[i] = 1.0 / H / (i - 1), i++)
            ;
        double sum = 0.0;
        for (int i = 2; i <= D + 1; sum += lamb[i] / i, i++)
            ;
        dist = vector<uint32_t>(D + 2);
        for (int i = 2; i <= D + 1;
             dist[i] = dist[i - 1] + MAX * lamb[i] / i / sum, i++)
            ;
        assert(MAX - dist[D + 1] < 10000);
        dist[D + 1] = MAX;
    }
    double getAvgDegree() {
        double sum = 0.0;
        double MAX = (double)std::numeric_limits<uint32_t>::max();
        for (int i = 0; i < dist.size(); ++i) {
            sum += 1.0 - dist[i] / MAX;
        }
        return sum + 3;
    }
    static uint32_t deg(const vector<uint32_t>& dist, uint32_t rnd) {
        uint32_t res = 0;
        for (; dist[res] < rnd; ++res)
            ;
        return res;
    }

    struct Configuration {
        double c;
        uint32_t D;
        friend ostream& operator<<(ostream& out, const Configuration& config) {
            return out << "(c = " << config.c << ", D = " << config.D << ")";
        }
    };

    struct Solution {
        vector<uint32_t> dist;
        vector<bool> v;
        uint32_t N;
    } sol;
    static size_t numBits(const Solution& sol) {
        return sol.dist.size() * sizeof(sol.dist[0]) * 8 + sol.v.size();
    }

    uint32_t n, N;
    using Edge = vector<uint32_t>;
    vector<Edge> edges;
    vector<bool> values;

    RetrieverLMSS(size_t m, Configuration config) { /* # ele/chunk to aim for */
        n = m / config.c;
        buildDistribution(config.D);
        N = n * (1 - 1.0 / config.D / config.D);
    }
    void addElement(const Hash& H, bool rhs) {
        /* note: I copy the vector of hashes here */
        int d = deg(dist, H.getFirstHash()) + 3;
        vector<uint32_t> vec = H.getMoreHashes(d);
        for (size_t i = 0; i < vec.size(); ++i) {
            vec[i] = i < 3 ? reduce(vec[i], n - N) + N : reduce(vec[i], N);
        }
        edges.push_back(std::move(vec));
        values.push_back(rhs);
    }

    bool runConstruction() {
        struct NodeInfo {
            uint32_t degree;
            uint32_t incidenceXOR;
        };
        vector<NodeInfo> nodes(n, {0, 0});
        for (int i = 0; i < (int)edges.size(); ++i) {
            for (int v : edges[i]) {
                nodes[v].degree++;
                nodes[v].incidenceXOR ^= i;
            }
        }
        stack<int> deg1Verts;
        for (int i = 0; (uint32_t)i < n; ++i) {
            if (nodes[i].degree == 1) {
                deg1Verts.push(i);
            }
        }

        vector<pair<uint32_t, uint32_t>> peelingOrder;
        while (!deg1Verts.empty()) {
            int v = deg1Verts.top();
            deg1Verts.pop();
            if (nodes[v].degree) {
                assert(nodes[v].degree == 1);
                int j = nodes[v].incidenceXOR;
                peelingOrder.push_back({j, v});
                Edge& e = edges[j];
                for (uint32_t i : e) {
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
            throw std::runtime_error("Failed to construct");
            return false;
        } else {
            // cout << "Success!" << endl;

            vector<bool> v(n);
            for (auto it = peelingOrder.rbegin(); it != peelingOrder.rend(); ++it) {
                Edge& e = edges[it->first];
                bool val = values[it->first];
                for (uint32_t i : e) {
                    val ^= v[i];
                }
                v[it->second] = val;
            }
            sol.v = std::move(v);
            sol.dist = std::move(dist);
            sol.N = N;
            return true;
        }
    }

    inline static bool retrieve(const Solution& sol, const Hash& H) {
        /* this heavily assumes little endian and stuff */
        /* this heavily assumes l = 8 and stuff */
        const vector<uint32_t>& vec =
            H.getMoreHashes(deg(sol.dist, H.getFirstHash()) + 3);
        bool res = false;
        for (uint32_t i = 0; i < 3; ++i) {
            uint32_t p = reduce(vec[i], sol.v.size() - sol.N) + sol.N;
            res ^= sol.v[p];
        }
        for (uint32_t i = 3; i < vec.size(); ++i) {
            uint32_t p = reduce(vec[i], sol.N);
            res ^= sol.v[p];
        }
        return res;
    }
};
} // namespace walzer
