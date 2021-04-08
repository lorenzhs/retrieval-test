#pragma once
#include<queue>
#include<cassert>
using namespace std;

namespace walzer {

/* very simple class that implements a priority queue
 * with a change key function. Keys are changed simply
 * by adding the element into the queue again. */
/* it is assumed that the elements 0,…,maxEle-1 are initially
 * in the queue with priority initialPriority each */

/* "[P]riority [Q]ueue with [C]hange [K]ey" */
class PQCK {
private:
    int nDistinctEntries;
    vector<float> lastKnownKey;
    typedef pair<float,int> PairType;
    typedef priority_queue<PairType> QType;
    QType q;
public:
    int size() { return nDistinctEntries; }
    PQCK(int maxEle) : nDistinctEntries(0), lastKnownKey(maxEle,-1.0) {}
    void insert(int ele, float key) {
        assert(lastKnownKey[ele] == -1.0);
        nDistinctEntries++;
        lastKnownKey[ele] = key;
        q.push({key,ele});
    }
    void changeKey(int ele, float newKey) {
        assert(lastKnownKey[ele] != -1.0);
        lastKnownKey[ele] = newKey;
        q.push({newKey,ele});
    }
    int pop() {
        assert(nDistinctEntries);
        int ele; float prio;
        do {
            assert(!q.empty());
            ele  = q.top().second;
            prio = q.top().first;
            q.pop();
        } while(lastKnownKey[ele] != prio);
        lastKnownKey[ele] = -1.0;
        nDistinctEntries--;
        return ele;
    }
};

}
