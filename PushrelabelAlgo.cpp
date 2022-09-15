#include "PushRelabelAlgo.h"
using namespace std;

void RushRelabelAlgo::init(MaxFlowStd *mx){
    int nv = mx->g.num_vertices;
    ex = (int*)calloc(nv, sizeof(int));
    localEx = (int*)calloc(nv, sizeof(int));
    localLabel = (int*)calloc(numVertices, sizeof(int));
    acumEx = (atomic_int*)calloc(nv, sizeof(atomic_int));
    h = (atomic_int*)calloc(nv,sizeof(atomic_int));

    res = (int**)calloc(numVertices, sizeof(int*));
    flows = (int**)calloc(numVertices, sizeof(int*));
    reverseMap = (vector<int> *)calloc(numVertices, sizeof(vector<int>));

    #pragma omp parallel for
    for(int i =0;i<nv;i++){
        res[i] = (int*)calloc(numVertices, sizeof(int));
        flows[i] = (int*)calloc(numVertices, sizeof(int));
    }
    activeSet.clear();

}
void RushRelabelAlgo::preflow(MaxFlowStd* mx){
    init(mx);
    //preflow
    int nv = mx->g.num_vertices;
    int ne = mx->g.num_edges;
    int sc = mx->sc;
    int sk = mx->sk;
    //original graph is labeled with capacity
    int** g = mx->g.capacities;

    //change source height
    h[sc] = nv;
    localLabel[sc] = nv;

    #pragma omp parallel for 
    for(int i = 0;i<nv;i++){
        if(g[sc][i] >0 && i != sc){
            flows[sc][i] = g[sc][i];
            flows[i][sc] = -[sc][i];
            res[sc][i] = g[sc][i]-flows[sc][i];
            ex[i] = g[sc][i];
            
        }
        //initialize other residual flow
        for(int j = 0;i<nv;j++){
            res[i][j] = g[i][j] - flows[i][j] > 0?g[i][j] - flows[i][j]:res[i][j];
        }
    }
    for(int i=0;i<nv;i++){
        if(i != sc && ex[i]>0) activeSet.insert(i);
        //??
        discoveredVertices.push_back(new tbb::concurrent_vector<int>());
    }

}
void PushRelabelAlgo::globalRelabel(int nv, int sc, int sk){
    //Any need ?
    #pragma omp parallel for 
    for(int i = 0; i < nv; i++){
        h[i] = nv;
    }
    h[sk] =0;
    //bfs, using queue
    tbb::concurrent_vector<int> queue;
    queue.push_back(sink);
    //Modification here!! cache locality or single thread?
    
    for(int i=0;i<nv;i++){
        for(int j =0;j<nv;j++){
            if(res[i][j]>0){
                reverseMap[j].push_back(i);
            }
        }
    }
    while(queue.size()>0){
        #pragma omp parallel for
        for(int i=0;i<queue.size();i++){
            int cur = q[i];
            discoveredVertices[cur].clear();
            //modify here
            for(int j =0;j<reverseMap[cur].size();j++){
                
                int nxt = reverseMap[cur][j];
                if(nxt != sc){
                    //Label possible edge on res network with distance to sink
                    if(nxt != sk && h[nxt].compare_exchange_strong(nv,h[cur]+1)){
                        discoveredVertices[cur].push_back(nxt);
                    }
                }
            }
            
        }
        tbb::concurrent_vector<int> tmp;
        

    }
}