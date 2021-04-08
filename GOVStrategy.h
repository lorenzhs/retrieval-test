#include<vector>
#include<iostream>
#include<algorithm>
#include<cassert>
#include<array>
#include<bitset>
#include<stack>
#include<tuple>
#include "util/line.h"
#include "util/util.h"
#include "util/pqck.h"
#include <boost/numeric/ublas/storage.hpp>

using boost::numeric::ublas::bounded_array;
using namespace std;

namespace walzer {

template<int k,typename Hashable>
struct GOVStrategy {
    static constexpr char stratName[] = "GOV";
	using Hashable_t = Hashable;
	struct Hash : DoubleHashSequence<k+1,Hashable> {
		Hash(const Hashable& str, uint32_t seed) : DoubleHashSequence<k+1,Hashable>(str, seed) {};
		uint32_t chunkHash() { return (*this)[k];  };
	};
	struct Configuration {
		double c;
		friend ostream& operator<<(ostream& out, const Configuration& config) {
			return out << "(c = " << config.c << ")";
		}
	};

	Line sol;
	uint8_t* data() { return (uint8_t*)sol.v.data(); };
	struct Meta {
		uint32_t sizeLogical;
		uint32_t sizeBytes;
	} meta;

	inline static bool retrieve(const Hash& H, const uint8_t* ptr, uint32_t logicalOffset, uint32_t logicalSize) {
		bool res = false;
		for (int i = 0; i < k; ++i) {
			uint32_t pp = reduce(H[i], logicalSize) + logicalOffset;
			res ^= bool(*(ptr + (pp >> 3)) & (1 << (pp & 0b111)));
		}
		return res;
	};

	typedef vector<Line> Matrix;

    typedef enum { EQ_SPARSE, EQ_SOLVED, EQ_DENSE } EqState;
    struct Equation {
        bounded_array<uint32_t,k> idleVars;
        EqState state;
        Line activeCoeff;
        bool rhs;
        void removeIdleVar(size_t vid) {
            auto it = std::find(idleVars.begin(),idleVars.end(),vid);
            assert(it != idleVars.end());
            assert(*it == vid);
            size_t s = idleVars.size();
            *it = idleVars[s-1];
            idleVars.resize(s-1);
        }
    };

    typedef enum { VAR_SOLVED, VAR_IDLE, VAR_DENSE } VarState;
    struct Variable {
        vector<uint32_t> eqs;
        VarState state;
    };

    vector<Equation> eqs;
    vector<Variable> vars;

    uint32_t nvars;
	GOVStrategy(size_t m, Configuration config) {
		nvars = m / config.c;
		nvars = (nvars + 7) & ~0b111; /* make divisible by 8 */
        vars.resize(nvars, {{},VAR_IDLE});
        eqs.reserve(m);
		meta.sizeBytes = nvars / 8;
		meta.sizeLogical = nvars;
    }

	void addElement(const Hash& H, bool rhs) {
		array<uint32_t, k> var;
        for (int i = 0; i < k; ++i) var[i] = reduce(H[i], nvars);
		/* make sure each varg occurs only once */
        sort(var.begin(),var.end());
        bounded_array<uint32_t,k> varUniq(0);
        for(unsigned i = 0; i < var.size(); ++i) {
            size_t s = varUniq.size();
            if (s > 0 && varUniq[s-1] == var[i])
                varUniq.resize(s-1);
            else {
                varUniq.resize(s+1);
                varUniq[s] = var[i];
            }
        }
        for(uint32_t vid : varUniq) {
            vars[vid].eqs.push_back(eqs.size());
        }
        eqs.emplace_back(Equation{varUniq,EQ_SPARSE,Line(),rhs});
    }

	bool runConstruction();
private:
	bool fourRussianGauss(/*in */Matrix &M,
           /*in */vector<bool> &rhss,
           /*out*/Line &solution);
};

template<int l,typename Hashable>
bool GOVStrategy<l,Hashable>::runConstruction() {
    /* Get Variables Sorted by #incident Equations */
    vector<size_t> varsSorted; varsSorted.reserve(vars.size()); {
        using p32 = pair<uint32_t,uint32_t>;
        vector<p32> degToVar; degToVar.reserve(vars.size());
        for(uint32_t i = 0; i < vars.size(); ++i) {
            degToVar.emplace_back(p32{(uint32_t)vars[i].eqs.size(),i});
        }
        sort(degToVar.begin(),degToVar.end(),std::greater<p32>());
        for(p32 &pair : degToVar) {
            varsSorted.push_back(pair.second);
        }
    }

    /* Lazy Gauss Solver */
    stack<size_t> lowDegEqs;
    for(size_t i = 0; i < eqs.size(); ++i) {
        if (eqs[i].idleVars.size() <= 1)
            lowDegEqs.push(i);
    }

    size_t nextActivation = 0;
    vector<size_t> activationOrder;
    vector<size_t> solveOrder;
    for(;;) {
        while(!lowDegEqs.empty()) {
            size_t eqID = lowDegEqs.top(); lowDegEqs.pop();
            Equation &eq(eqs[eqID]);
            if(eq.idleVars.size() == 0) {
                eq.state = EQ_DENSE;
                continue;
            }
            assert(eq.idleVars.size() == 1);
            size_t vid = eq.idleVars[0];
            Variable &v(vars[vid]);
            for(size_t otherEqID : v.eqs) {
                if(otherEqID == eqID) continue;
                Equation &otherEq(eqs[otherEqID]);
                otherEq.removeIdleVar(vid);
                otherEq.rhs ^= eq.rhs;
                otherEq.activeCoeff ^= eq.activeCoeff;
                if (otherEq.idleVars.size() == 1) {
                    lowDegEqs.push(otherEqID);
                }
            }
            v.state = VAR_SOLVED;
            eq.state = EQ_SOLVED;
            solveOrder.push_back(eqID);
        }
        if (nextActivation == varsSorted.size()) break;
        size_t vid = varsSorted[nextActivation++];
        Variable &v(vars[vid]);
        if (v.state == VAR_SOLVED) continue;
        assert(v.state == VAR_IDLE);
        if(v.eqs.size() == 0) break;
        v.state = VAR_DENSE;
        activationOrder.push_back(vid);
        for(size_t eqID : v.eqs) {
            Equation &eq(eqs[eqID]);
            eq.activeCoeff.resize(activationOrder.size());
            eq.activeCoeff.setBit(activationOrder.size()-1,1);
            eq.removeIdleVar(vid);
            if (eq.idleVars.size() == 1) {
                lowDegEqs.push(eqID);
            }
        }
    }

    /* Solve the Core using the Method of Four Russians */
    Matrix M;
    vector<bool> rhss;
    for(Equation &eq: eqs) {
        if (eq.state == EQ_DENSE) {
            eq.activeCoeff.resize(activationOrder.size());
            M.emplace_back(std::move(eq.activeCoeff));
            rhss.push_back(eq.rhs);
        }
    }

    Line activeVarSol(activationOrder.size());
    if(!fourRussianGauss(M,rhss,activeVarSol)) return false;
    sol.resize(nvars);

    /* attach solutions to corresponding variables */
    for(size_t i = 0; i < activationOrder.size(); ++i) {
        size_t vid = activationOrder[i];
        sol.setBit(vid,activeVarSol[i]);
    }
    /* back substitution */
    for(int i = solveOrder.size()-1; i >= 0; --i) {
        size_t eqID = solveOrder[i];
        Equation &eq(eqs[eqID]);
        bool bit = scalarProduct(eq.activeCoeff,activeVarSol)^eq.rhs;
        sol.setBit(eq.idleVars[0],bit);
    }

    return true;
}

/* Note: GOV do not implement the Method of Four Russians
 * but this speeds up construction, so I added it. */
template<int l, typename Hashable>
bool GOVStrategy<l,Hashable>::fourRussianGauss(/*in */Matrix &M,
           /*in */vector<bool> &rhss,
           /*out*/Line &solution) {
    int n = M.size();
    if (n == 0) return true;
    int m = M[0].size();


    int r = 0, c = 0;
    auto chooseK = [&]() -> int {
        return (m - c <= 3 ? 1 : log2(m - c) - log2(log2(m - c)));
    };

    vector<int> pivots;
    while ((r < n) && (c < m)) {
        int k = chooseK();
        assert(k >= 1 && k <= m-c);
        vector<Line> T(1 << k);
        vector<bool> Trhss(1 << k);
        int d = 0;
        Line::BitsView view(c,k);
        for (int rr = r; rr < n; ++rr) {
            uint32_t p = M[rr].get(view);
            if (!p) continue;

            if (T[p].size()) {
                M[rr].addFrom(T[p],c);
                assert(M[rr].get(view) == 0);
                rhss[rr] = rhss[rr] ^ Trhss[p];
            } else {
                for (uint32_t pp = 1; (int)pp < (1 << k); ++pp) {
                    uint32_t newp = p ^ pp;
                    if (T[pp].size() && !T[newp].size()) {
                        T[newp].copyFrom(T[pp]);
                        T[newp].addFrom(M[rr],c);
                        Trhss[newp] = rhss[rr] ^ Trhss[pp];
                    }
                }
                swap(T[p], M[rr]); // now M[rr] is empty
                assert(!M[rr].size());
                swap(Trhss[p], rhss[rr]);
                swap(M[rr], M[r + d]); // now M[i + d] is empty
                assert(!M[r+d].size());
                swap(rhss[rr], rhss[r + d]);
                d++;
            }
        }
        /* get basis of the T equations */
        vector<int> pivotPattern(k,-1);
        for (uint32_t p = 1; (int)p < (1 << k); ++p)
            if (T[p].size())
                pivotPattern[ctz(p)] = p;

        [[maybe_unused]] int oldr = r;
        for(int cc = 0; cc < k; ++cc) {
            int p = pivotPattern[cc];
            if (p != -1) {
                pivots.push_back(c+cc);
                swap(rhss[r],Trhss[p]);
                assert(!M[r].size());
                swap(M[r++],T[p]);
            }
        }
        c += k;
        assert(oldr + d == r);
        assert((int)pivots.size() == r);
    }

    /* matrix is now in normal form */
    /* remaining equations (beyond r) have zero left hand sides
     * and should have zero right hand sides as well */
    #ifdef ABORT_IF_DEPENDENT
        if (r < n) return false;
    #else
        for(int rr = r; rr < n; ++rr) {
            /* this is an or */
            if (rhss[rr]) return false;
        }
    #endif
    /* system is solvable; infer one solution,
     * initialise variables from left to right. */
    for(int rr = r-1; rr >= 0; --rr) {
        solution.setBit(pivots[rr], scalarProduct(solution, M[rr], pivots[rr]) ^ rhss[rr]);
    }

    return true;
}
}
