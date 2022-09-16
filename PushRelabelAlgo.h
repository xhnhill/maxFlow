#include "baseClass.h"
#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_unordered_set.h"
#include <atomic>
#include <vector>
#include <unordered_set>
#include <chrono>
using namespace std;

struct Graph{
  int num_vertices;
  int num_edges;
  int** capacities; // num_vertices x num_vertices array where cij = flow capacity of edge i->j
};
class MaxFlowStd{
    public:
    Graph g;
    int sc;
    int sk;
};
class MaxFlowRes{
    public:
    int val;
    int** flow;
};
class Problem{
  public:
  Graph g;
  int** compareMatrix;
  int res;
  int actualNum;
};
class Clock{
    public:
    Clock() : flg(chrono::high_resolution_clock::now()) {}
    void reset() {flg = chrono::high_resolution_clock::now();}
    double duration() const{
        return chrono::duration_cast<chrono::duration<double, std::ratio<1>>>(chrono::high_resolution_clock::now()-flg).count();
    }
    private:
    chrono::time_point<chrono::high_resolution_clock> flg;

};
class PushRelabelAlgo{
    public:
    Clock c;
     double pushRelabel(MaxFlowStd *mx,MaxFlowRes *result);
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
    

};
