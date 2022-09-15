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