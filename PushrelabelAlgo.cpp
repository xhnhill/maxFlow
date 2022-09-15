#include "PushRelabelAlgo.h"
using namespace std;

void RushRelabelAlgo::init(MaxFlowStd *mx){
    int nv = mx->g.num_vertices;
    ex = (int*)calloc(nv, sizeof(int));
    localEx = (int*)calloc(nv, sizeof(int));
    localLabel = (int*)calloc(nv, sizeof(int));
    acumEx = (atomic_int*)calloc(nv, sizeof(atomic_int));
    h = (atomic_int*)calloc(nv,sizeof(atomic_int));

    res = (int**)calloc(nv, sizeof(int*));
    flows = (int**)calloc(nv, sizeof(int*));
    reverseMap = (vector<int> *)calloc(nv, sizeof(vector<int>));

    #pragma omp parallel for
    for(int i =0;i<nv;i++){
        res[i] = (int*)calloc(nv, sizeof(int));
        flows[i] = (int*)calloc(nv, sizeof(int));
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
            flows[i][sc] = -g[sc][i];
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
        //Difference??
        #pragma omp parallel for
        for(int i =0;i<queue.size();i++){
            int cur = queue[i];
            for(int j =0;j<discoveredVertices[cur].size();j++){
                tmp.push_back(discoveredVertices[cur]);
            }
        }
        queue.swap(tmp);

    }
}
void PushRelabelAlgo::pushRelabel(MaxFlowStd* mx,MaxFlowRes *rest){
    //begin record
    c.reset();
    //statistics which helps to decide if need bfs, introduced in prsn
    int workSinceLastGR = INT_MAX;
    //Use default para used by prsn
    float freq = 0.5; 
    int alp = 6; 
    int beta = 12;
    preflow(mx);
    int nv = mx->g.num_vertices;
    int ne = mx->g.num_edges;
    int** g = mx->g.capacities;
    int sc = mx->sc;
    int sk = mx->sk;

    while(true){
        //Modify freq * workSinceLastGR > a * numVertices + numEdges
        //original prsn algo has some mistakes,
        //one solution is run bfs after each iteration
        //the other is fix the relabel condition
        // Here for the performance, use the later one
        if(freq * workSinceLastGR > a * nv + ne){

            workSinceLastGR = 0;
            globalRelabel(nv, sc, sk);
            //Test parallel performance here
            #pragma omp parallel for
            for(int i = 0; i < nv; i++){
            localLabel[i] = h[i];
            }
            unordered_set<int> tmpSet;
            for(int v : workingSet){
            if(h[v] < nv && v != mx->sk){
                tmpSet.insert(v);
            }
            }
        activeSet.swap(tmpSet);
        }
        if(activeSet.empty()) break;
        #pragma omp parallel
        {
            #pragma omp single
            {
                for(auto itr = activeSet.begin();itr != activeSet.end();itr++){
                    
                }
            }
        }
        
    }

}