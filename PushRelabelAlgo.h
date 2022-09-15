#include "baseClass.h"
#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_unordered_set.h"
#include <atomic>
#include <vector>
#include <unordered_set>
using namespace std;

class RushRelabelAlgo{
    public:
    Clock c;
     void pushRelabel(MaxFlowStd *mx,MaxFlowRes *result);
    private:
    //excess val of vertex
    int* ex;
    atomic_int* h;
    //matrix representation of residual network
    int** res;
    //Accumulated excess flow during one iteration
    atomic_int* acumEx;
    int* localLabel;
    int* localEx;
    vector<int> *reverseMap;
    vector<tbb::concurrent_vector<int>> discoveredVertices;
    unordered_set<int> activeSet;
    //Dynamic flow in graph
    int** flows;
    //statistic for the bfs condition
    int* work;
    void init(MaxFlowStd* mx);
    void preflow(MaxFlowStd* mx);
    void globalRelabel(int numV, int sc, int sk);
    void tryGlobalRelabel(MaxFlowStd *mx, MaxFlowStd *res);
    

}
