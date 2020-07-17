#include<vector>
#include<iostream>
#include<iomanip>
#include<queue>
#include<algorithm>
#include<cassert>
#include<array>
#include<bitset>
#include<stack>
#include<tuple>
#include "util/line.h"
#include "util/util.h"
#include "util/pqck.h"

#ifndef DEBUG
    #define DEBUG(x)   // do nothing
#endif
#ifndef DEBUG_DETAILED
    #define DEBUG_DETAILED(x) // do nothing
#endif
#ifndef PROFILE
    #define PROFILE(x)
#endif

using namespace std;

/* Todo: Do something reasonable for l > 16 */
template<int l, typename Hashable>
struct TwoBlockStrategy {
    static constexpr char stratName[] = "TwoBlock";
	using Hashable_t = Hashable;
	struct Hash {
		uint32_t ch,b1,b2;
		uint16_t p1, p2;
		static_assert(l <= 16);
		Hash(const Hashable& s, uint32_t seed) {
			static_assert(sizeof(Hash) == 128 / 8);
			do_hash(s, seed, this);
		};
		uint32_t chunkHash() const { return ch; };
	};
	struct Configuration {
		double c;
		friend ostream& operator<<(ostream& out, const Configuration& config) {
			return out << "(c = " << config.c << ")";
		}
	};
	
	struct Meta {
		uint32_t sizeLogical;
		uint32_t sizeBytes;
	} meta;

	TwoBlockStrategy(size_t m, Configuration config) : activationQueue(0) {
		meta.sizeLogical = (m / config.c + l-1) / l;
		/* todo: get rid of this assumption */
		static_assert(l % 8 == 0);
		meta.sizeBytes = (meta.sizeLogical * l) / 8;
		activationQueue = PQCK(meta.sizeLogical);
		vargs.resize(meta.sizeLogical);
		eqs.reserve(m);
	}
	void addElement(const Hash& H, bool rhs) {
		int vgid1 = reduce(H.b1, meta.sizeLogical);
		int vgid2 = reduce(H.b2, meta.sizeLogical);
		bitset<l> p1(H.p1);
		bitset<l> p2(H.p2);
		/* make sure each varg occurs only once */
		if (vgid1 == vgid2 || p1.none()) {
			vgid1 = vgid2;
			p1 ^= p2;
			p2.reset();
		}
		int eqID = eqs.size();
		if (p1.none())
			vgid1 = -1;
		else {
			vargs[vgid1].lines.push_back(p1);
			vargs[vgid1].eqIDs.push_back(eqID);
		}
		if (p2.none())
			vgid2 = -1;
		else {
			vargs[vgid2].lines.push_back(p2);
			vargs[vgid2].eqIDs.push_back(eqID);
		}

		eqs.emplace_back(Equation{ {vgid1,vgid2},EQ_SPARSE,Line(), rhs });
	}
	bool runConstruction();

	vector<uint8_t> sol;
	void* data() {
		/* this copying step could be avoided, but is probably not all that expensive. */
		if (!sol.size()) {
			sol.resize(meta.sizeBytes);
			uint8_t* ptr = sol.data();
			int bit = 0;
			for (auto& vg : vargs) {
				for (int i = 0; i < l; ++i) {
					*ptr |= vg.sol[i] << bit++;
					ptr += bit >> 3; bit &= 7;
				}
			}
		}
		return sol.data();
	};
	inline static bool retrieve(const Hash& H, const void* ptr, uint32_t logicalOffset, uint32_t logicalSize) {
		static_assert(min(sizeof(unsigned),sizeof(H.p1)) * 8 >= l);
		unsigned r  = extract<l>(ptr, logicalOffset + reduce(H.b1, logicalSize)) & H.p1;
				 r ^= extract<l>(ptr, logicalOffset + reduce(H.b2, logicalSize)) & H.p2;
		return parityl(r);
	}

private:
	typedef vector<Line> Matrix;
    typedef int EqID;
    typedef int VGID;
    typedef struct {
        int g;
        int o;
    } VID;
    
    /* PeelInfo is dictates that the variable with id vid
     * must be initalised to fulfill a certain equation. */
    struct PeelInfo{
        int eqID;
        bitset<l> p;
    };
    stack<PeelInfo> peelLog;
    
    void lazyGauss();
        bitset<l> solveEquation(int eqID);
        void printState();
    void backSubstitution();
    bool solveCore();
    bool gauss(/*in */Matrix &M,
               /*in */vector<bool> &rhss,
               /*out*/Line &solution);
    bool fourRussianGauss(/*in */Matrix &M,
           /*in */vector<bool> &rhss,
           /*out*/Line &solution);
    
    /* Equations are sparse in the beginning.
     * Solved Equations need never be considered again,
     * an entry in the peelLog guarantees that they will be fulfilled.
     * Dense Equations will be contained in the core */
    typedef enum { EQ_SPARSE, EQ_SOLVED, EQ_DENSE } EqState;
    struct Equation {
        int degree() { return (vgids[0] != -1) + (vgids[1] != -1); }
        void removeVGID(VGID vgid) {
            if (vgid == vgids[0]) {
                vgids[0] = vgids[1], vgids[1] = -1;
            } else {
                assert(vgids[1] == vgid);
                vgids[1] = -1;
            }
            if (degree() == 0) state = EQ_DENSE;
        }
        VGID vgids[2];
        EqState state;
        Line activeCoeff;
        bool rhs;
    };
    
    struct VGroup {
        VGroup() : activePos(), peeled(), activated(false) {};
        float priority() {
            if (activated || peeled.all()) return 0;
            return eqIDs.size() / (float)popcount<l>(~peeled);
        }
        int findLine(EqID eqID) {
            int pos = std::find(eqIDs.begin(), eqIDs.end(), eqID) - eqIDs.begin();
            assert(pos < (int)eqIDs.size());
            return pos;
        }
        vector<EqID> eqIDs;
        vector<bitset<l>> lines;
        bitset<l> sol;
        array<int,l> activePos;
        bitset<l> peeled;
        bool activated;
    };
    stack<int> lowDegEqs;
    PQCK activationQueue;
    Line activeVarSol;
    vector<VID> activeVarToVID;
    vector<Equation> eqs; /* constantly modified */
    vector<VGroup> vargs;
    
    /* profiling */
    PROFILE(
        typedef std::chrono::high_resolution_clock Clock;
        typedef std::chrono::microseconds Duration;
        typedef Clock::time_point TimePoint;
        Duration totalTime; Duration coreTime;
        int coreSystemEqs; int coreSystemVars;
    )

    void clearAuxData() {
        this->data(); // rescues all data, if present
        assert(peelLog.size() == 0);
        assert(lowDegEqs.size() == 0);
        activationQueue = PQCK(0);
        activeVarSol = Line();
        activeVarToVID.clear();
        eqs.clear();
        vargs.clear();
    }
};

template<int l,typename Hashable>
bool TwoBlockStrategy<l,Hashable>::runConstruction() {
    PROFILE(TimePoint start = Clock::now();)
        DEBUG(cout << "% Starting Blocked Structured Solver with m = " << eqs.size() << " equations on n = " << l * vargs.size() << " variables." << endl;)
        lazyGauss();
        DEBUG(cout << "\\end{state}\\endinput" << endl << "Peeled " << peelLog.size() << " variables" << endl;)
        bool solved = solveCore(); /* calls gauss */
    if (!solved) { 
        DEBUG(cout << "Core System is not solvable. Aborting." << endl;)
    } else {
        DEBUG(cout << "Done. Doing backsubstitution... ";)
            backSubstitution();
        DEBUG(cout << "Done." << endl;)
        DEBUG_DETAILED(
            cout << "solution: ";
            for(VGroup &vg : vargs) {
                for(int i = l-1; i >= 0; --i) {
                    cout << vg.sol[i];
                } cout << "|";
            } cout << endl;
        )
        /* solution ready in this->vargs[...].sol
         * but copy it into a vector and free all the memory */
        //clearAuxData();
    }
    PROFILE(
        TimePoint end = Clock::now();
        totalTime = std::chrono::duration_cast<Duration>(end - start);
    )
    return solved;
}

template<int l, typename Hashable>
bitset<l> TwoBlockStrategy<l,Hashable>::solveEquation(int eqID) {
    Equation &eq(eqs[eqID]);
    assert(eq.degree() == 1 && eq.vgids[0] != -1);
    /* picks a variable in that group
     * solves for it,
     * remove the variable from the group
     * returns the id of the variable and the pattern of the eq */
    VGID vgid(eq.vgids[0]);
    VGroup &vg(vargs[vgid]);
    
    /* find and remove the corresponding line in the varGroup */
        int pos = vg.findLine(eqID);
        bitset<l> p(vg.lines[pos]);
        assert(p.any());
        int bitIdx = ctz(p.to_ulong());
        DEBUG_DETAILED(cout << "% for bit " << bitIdx << "(vgid " << vgid<< ")" << endl;)
        int bit    = 1 << bitIdx;
        // remove
        vg.lines[pos] = vg.lines.back(); vg.lines.pop_back();
        vg.eqIDs[pos] = vg.eqIDs.back(); vg.eqIDs.pop_back();
    
    /* add the equation to each equation involving the chosen bit */
    for(int i = 0; i < (int)vg.eqIDs.size(); ++i) {
        if (!(vg.lines[i][bitIdx])) continue;
        vg.lines[i] ^= p;
        Equation &otherEq(eqs[vg.eqIDs[i]]);
        otherEq.rhs ^= eq.rhs;
        otherEq.activeCoeff ^= eq.activeCoeff;
        
        /* remove the equation, if needed */
        if (vg.lines[i].none()) {
            otherEq.removeVGID(vgid);
            lowDegEqs.push(vg.eqIDs[i]);
            vg.lines[i] = vg.lines.back(); vg.lines.pop_back();
            vg.eqIDs[i] = vg.eqIDs.back(); vg.eqIDs.pop_back();
            i--;
        }
    }
    vg.peeled |= bit;
    activationQueue.changeKey(vgid,vg.priority());
    
    return p;
}

template<int l, typename Hashable>
void TwoBlockStrategy<l, Hashable>::printState() {
    vector<int> compressedPos(vargs.size());
    int nactive = 0;
    for(int i = 0; i < (int)vargs.size(); ++i) {
        if (!vargs[i].activated && vargs[i].eqIDs.size()) {
            compressedPos[i] = nactive++;
        }
    }
    cout << "\\begin{state}\\begin{tabular}{";
    cout << "cc"; for(int i = 0; i < nactive; cout << (i++ ? ":" : "|") << "c");
    cout << "|c}" << endl;
    
    string S(nactive*(l+1),' ');
    for(int i = 0; i < nactive; ++i) {
        S[i*(l+1)+l] = '&';
    }
    
    cout << "&\\densecol";
    for(int i = 0; i < (int)vargs.size(); ++i) {
        if (!vargs[i].activated && vargs[i].eqIDs.size()) {
            cout << "&\\varg{" << i << "}";
        }
    }
    cout << "&\\rhscol\\\\" << endl;
    
    for(int state : {EQ_DENSE,EQ_SPARSE}) {
        cout << "\\hline" << endl;
        for(int eqID = 0; eqID < (int)eqs.size(); ++eqID) {
            Equation &eq(eqs[eqID]);
            if (eq.state != state) continue;
            cout << "\\eqid{" << eqID << "} &";
            /* first the active vars */
            eq.activeCoeff.resize(activeVarToVID.size());
            for(int i = 0; i < (int)activeVarToVID.size(); ++i) {
                cout << eq.activeCoeff[i];
            } cout << "&";
            /* now the sparse part */
            string line = S;
            for(int j = 0; j < 2; ++j) {
                if (eq.vgids[j] != -1) {
                    VGID vgid(eq.vgids[j]);
                    VGroup &vg(vargs[vgid]);
                    stringstream ss;
                    ss << vg.lines[vg.findLine(eqID)];
                    string sss = ss.str();
                    for(int i = 0; i < l; ++i) {
                        if (vg.peeled[l-1-i]) {
                            sss[i] = '*';
                        }
                    }
                    memcpy(&line[(l+1)*compressedPos[vgid]],&sss[0],l);
                }
            }
            cout << line << " \\rhs{" << eq.rhs << "}\\\\" << endl;
        }
    }
    cout << "\\hline\\end{tabular}" << endl << endl;
};

template<int l, typename Hashable>
void TwoBlockStrategy<l, Hashable>::lazyGauss() {
    for(EqID eqID = 0; eqID < (int)eqs.size(); ++eqID) {
        if (eqs[eqID].degree() < 2) {
            lowDegEqs.push(eqID);
        }
    }
    
    for(VGID vgid = 0; vgid < (int)vargs.size(); ++vgid) {
        activationQueue.insert(vgid,vargs[vgid].priority());
    }
    
    DEBUG_DETAILED(bool change = true;)
    for(;;) {
        DEBUG_DETAILED(if(change) printState(); change = false;)
        if(!lowDegEqs.empty()) {
            int eqID = lowDegEqs.top(); lowDegEqs.pop();
            Equation &eq(eqs[eqID]);
            if (eq.state == EQ_SOLVED || eq.state == EQ_DENSE)
                continue;
            if (eq.degree() == 0) {
                eq.state = EQ_DENSE;
                continue;
            }
            DEBUG_DETAILED(change = true; cout << "\\solve{" << eqID << "}\\end{state}" << endl;)
            bitset<l> p = solveEquation(eqID);
            peelLog.push({eqID, p});
            eq.state = EQ_SOLVED;
        } else {
            if (!activationQueue.size()) {
                break;
            }
            int vgid = activationQueue.pop();
            VGroup &vg(vargs[vgid]);
            DEBUG_DETAILED(if (vg.eqIDs.size()) change = true && cout << "\\activate{" << vgid << "}\\end{state}" << endl;)
            assert(vg.activated == false);
            // only activate variables that are actually used:
            bitset<l> used(0);
            for(bitset<l> line : vg.lines) used |= line;
            assert((used & vg.peeled).none());
            bitset<l> activate = used;
            int firstVarOffset = activeVarToVID.size();
            for(int i = 0; i < l; ++i) {
                if (activate[i]) {
                    vg.activePos[i] = activeVarToVID.size();
                    activeVarToVID.push_back(VID{vgid,i});
                }
            }
            
            for(int i = 0; i < (int)vg.eqIDs.size(); ++i) {
                EqID eqID(vg.eqIDs[i]);
                Equation &eq(eqs[eqID]);
                eq.activeCoeff.resize(activeVarToVID.size());
                /* put that part in the bitvector */
                int offset = firstVarOffset;
                for(int j = 0; j < l; ++j) {
                    if (activate[j]) {
                        eq.activeCoeff.setBit(offset++,vg.lines[i][j]);
                    }
                }
                eq.removeVGID(vgid);
                lowDegEqs.push(eqID);
            }
            vg.activated = true;
            vg.eqIDs.clear();
            vg.lines.clear();
        }
    }
}

template<int l, typename Hashable>
void TwoBlockStrategy<l, Hashable>::backSubstitution() {
    /* work backwards through equations declared solved */
    /* compute unique value for the single variable that fulfills the eq */


    /* assume activeSols is initialised and the information
     * from there already propagated to the respective vargs. */
    
    while(!peelLog.empty()) {
        PeelInfo pi(peelLog.top()); peelLog.pop();
        Equation &eq(eqs[pi.eqID]);
        VGroup   &vg(vargs[eq.vgids[0]]);
        
        bool activeProd = scalarProduct(eq.activeCoeff, activeVarSol);
        bool sparseProd = scalarProduct(pi.p.to_ulong(), vg.sol.to_ulong());
        bool toggle = activeProd ^ sparseProd ^ eq.rhs;
        int bitIdx = ctz(pi.p.to_ulong());
        DEBUG_DETAILED(
            cout << "Equation " << pi.eqID << (toggle ? " not " : " already ") << "fulfilled. ";
            cout << "concerns bit " << bitIdx << " of varg " << eq.vgids[0] << endl;
        )
        vg.sol ^=  (toggle << bitIdx);
    }
}

template<int l, typename Hashable>
bool TwoBlockStrategy<l, Hashable>::solveCore() {
    Matrix M;
    vector<bool> rhss;
    for(Equation &eq: eqs) {
        if (eq.state == EQ_DENSE) {
            eq.activeCoeff.resize(activeVarToVID.size());
            M.emplace_back(std::move(eq.activeCoeff));
            rhss.push_back(eq.rhs);
        }
    }
    
    DEBUG(cout << "Solving core system of size " << M.size() << " x " << (M.size() ? M[0].size() : 0) << endl;)
    PROFILE((coreSystemEqs = M.size(), coreSystemVars = activeVarToVID.size());)
    activeVarSol.resize(activeVarToVID.size());
    PROFILE(TimePoint coreStart = Clock::now();)
    bool res = fourRussianGauss(M,rhss,activeVarSol);
    PROFILE(TimePoint coreEnd = Clock::now();
            coreTime  = std::chrono::duration_cast<Duration>(coreEnd - coreStart);)
    if (!res) return false;
    
    DEBUG_DETAILED(cout << "Core Solution (activation order): ";)

    /* attach solutions to corresponding variables */
    for(int i = 0; i < (int)activeVarToVID.size(); ++i) {
        VID vid = activeVarToVID[i];
        VGroup &vg(vargs[vid.g]);
        assert(!vg.sol[vid.o]);
        assert(!vg.peeled[vid.o]);
        vg.sol[vid.o] = activeVarSol[i];
        DEBUG_DETAILED(cout << activeVarSol[i];)
    }
    DEBUG_DETAILED(cout << endl;)
    return true;
}

template<int l, typename Hashable>
bool TwoBlockStrategy<l, Hashable>::fourRussianGauss(/*in */Matrix &M,
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
        
        int oldr = r;
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

template<int l, typename Hashable>
bool TwoBlockStrategy<l, Hashable>::gauss(/*in */Matrix &M,
           /*in */vector<bool> &rhss,
           /*out*/Line &solution) {
    int n = M.size();
    if (n == 0) return true;
    int m = M[0].size();

    vector<int> pivots;
    int c = 0, r = 0;
    while(c < m && r < n) {
        int rr = r;
        for(; rr < n && !M[rr][c]; ++rr);
        if (rr == n) {
            c++;
            continue;
        }
        pivots.push_back(c);
        swap(M[rr], M[r]);
        bool tmp = rhss[rr];
        rhss[rr] = rhss[r];
        rhss[r]  = tmp;
        for(rr = r+1; rr < n; ++rr) {
            if (M[rr][c]) {
                M[rr] ^= M[r];
                rhss[rr] = rhss[rr] ^ rhss[r];
            }
        }
        c++,r++;
    }

    /* matrix is now in normal form */
    /* remaining equations (beyond r) have zero left hand sides
     * and should have zero right hand sides as well */
    for(int rr = r; rr < n; ++rr) {
        /* this is an or */
        if (rhss[rr]) return false;
    }
    
    /* system is solvable; infer one solution,
     * initialise variables from left to right. */
    for(int rr = r-1; rr >= 0; --rr) {
        solution.setBit(pivots[rr], scalarProduct(solution, M[rr], pivots[rr]) ^ rhss[rr]);
    }
    return true;
}
