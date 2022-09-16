#include <iostream>
#include "PushRelabelAlgo.h"
using namespace std;

void PushRelabelAlgo::init(MaxFlowStd *mx){
    int nv = mx->g.num_vertices;
    ex = (int*)calloc(nv, sizeof(int));
    localEx = (int*)calloc(nv, sizeof(int));
    localLabel = (int*)calloc(nv, sizeof(int));
    acumEx = (atomic_int*)calloc(nv, sizeof(atomic_int));
    h = (atomic_int*)calloc(nv,sizeof(atomic_int));

    res = (int**)calloc(nv, sizeof(int*));
    flows = (int**)calloc(nv, sizeof(int*));
    reverseMap = (vector<int> *)calloc(nv, sizeof(vector<int>));

    work = (int*)calloc(nv, sizeof(int));
    #pragma omp parallel for
    for(int i =0;i<nv;i++){
        res[i] = (int*)calloc(nv, sizeof(int));
        flows[i] = (int*)calloc(nv, sizeof(int));
    }
    activeSet.clear();

}
void PushRelabelAlgo::preflow(MaxFlowStd* mx){
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
        for(int j = 0;j<nv;j++){
            res[i][j] = g[i][j] - flows[i][j] > 0?g[i][j] - flows[i][j]:res[i][j];
        }
    }
    tbb::concurrent_vector<int> eptVec = {};
    for(int i=0;i<nv;i++){
        discoveredVertices.push_back(eptVec);
        if(i != sc && ex[i]>0) activeSet.insert(i);
        
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
    queue.push_back(sk);
    //Modification here!! cache locality or single thread?
    
    for(int i=0;i<nv;i++){
        for(int j =0;j<nv;j++){
            if(res[i][j]>0){
                reverseMap[j].push_back(i);
            }
        }
    }
    while(!queue.empty()){
        #pragma omp parallel for
        for(int i=0;i<queue.size();i++){
            int cur = queue[i];
            discoveredVertices[cur].clear();
            //modify here
            for(int j =0;j<reverseMap[cur].size();j++){
                
                int nxt = reverseMap[cur][j];
                if(nxt != sc){
                    int rpVal = nv;
                    //Label possible edge on res network with distance to sink
                    if(nxt != sk && h[nxt].compare_exchange_strong(rpVal,h[cur]+1)){
                        discoveredVertices[cur].push_back(nxt);
                        //cout<<cur<<" "<<nxt<<"\n";
                    }
                }
            }
            
        }
        tbb::concurrent_vector<int> tmp;
        //Difference??
        #pragma omp parallel for
        for(int i=0;i<queue.size();i++){
            int cur = queue[i];
            for(int j =0;j<discoveredVertices[cur].size();j++){
                tmp.push_back(discoveredVertices[cur][j]);
                //cout<<"push "<<discoveredVertices[cur][j]<<"\n";
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
        if(freq * workSinceLastGR > alp * nv + ne){

            workSinceLastGR = 0;
            globalRelabel(nv, sc, sk);
            //Test parallel performance here
            #pragma omp parallel for
            for(int i = 0; i < nv; i++){
            localLabel[i] = h[i];
            }
            unordered_set<int> tmpSet;
            for(int v : activeSet){
            if(h[v] < nv && v != sk){
                tmpSet.insert(v);
            }
            }
        activeSet.swap(tmpSet);
        }
        if(activeSet.empty()) {
            break;
        }
        #pragma omp parallel
        {
            #pragma omp single
            {
                for(auto itr = activeSet.begin();itr != activeSet.end();itr++){
                   #pragma omp task untied
                   {
                    int v = *itr;
                    //init work
                    work[v] = 0;
                    
                    discoveredVertices[v].clear();
                    localLabel[v] = h[v];
                    //local copy of excess of cur vertex
                    int le = ex[v];
                    while(le >0){
                        int nb = nv;
                        //Use to resolve conflict between neighboring vertex
                        bool skip = false;
                        int scaned = 0;
                        bool hasPush = false;

                        for(int i =0;i<nv;i++){
                            if(res[v][i]>0){
                                if(le == 0) {
                                    break;
                                }
                                scaned++;
                                bool admissible = (localLabel[v] == h[i]+1);
                                //resolve conflicts, in this case the concurrent version always
                                // has sequential replacement
                                if(ex[i]>0){
                                   bool win = h[v] ==h[i]+1 || h[v]<h[i]-1 || (h[v] == h[i] && v<i);
                                   if(admissible && (!win)){
                                    skip = true;
                                    continue;
                                   }
                                }
                                if(admissible){
                                    int dlt = min(res[v][i],le);
                                    flows[v][i] += dlt;
                                    flows[i][v] -= dlt;
                                    le -= dlt;
                                    atomic_fetch_add(&acumEx[i], dlt);
                                    if (i != sk){
                                        discoveredVertices[v].push_back(i);
                                    }
                                    hasPush = true;
                                }
                                //Modification to the original algo
                                //The original has mistakes here
                                //Should use runtime residual network to check
                                if (le>0 && g[v][i]-flows[v][i]>0) {
                                    nb = min(nb, h[i] + 1);
                                }
                            }
                        }
                        if(le == 0 || skip){
                            break;
                        }
                        localLabel[v] = nb;
                        work[v] = work[v]+scaned+beta;
                        if(localLabel[v] == nv){
                            break;
                        }
                    }
                    acumEx[v] = acumEx[v] + (le - ex[v]);
                    if(le>0){
                        discoveredVertices[v].push_back(v);
                    }

                   } 
                }
            }
            #pragma omp taskwait
        }
        #pragma omp parallel for
        for(int i=0;i<nv;i++){
            h[i] = localLabel[i];
            ex[i] = ex[i]+acumEx[i];
            acumEx[i] = 0;
            reverseMap[i].clear();
        }
        unordered_set<int> swpSet;
        for(int i:activeSet){
            workSinceLastGR += work[i];
            for(int j =0;j<discoveredVertices[i].size();j++){
                int cur = discoveredVertices[i][j];
                if(h[cur]<nv){
                    res[i][cur] = g[i][cur]-flows[i][cur];
                    res[cur][i] = g[cur][i]-flows[cur][i];
                    swpSet.insert(cur);
                }
            }
            res[i][sk] = g[i][sk] - flows[i][sk];
            res[sk][i] = g[sk][i] - flows[sk][i];

        }
        activeSet.swap(swpSet);
        #pragma omp parallel
        {
           #pragma omp single
           {
            for (auto itr = activeSet.begin(); itr != activeSet.end(); itr++){
                #pragma omp task
                {
                    int cur = *itr;
                    ex[cur] = ex[cur]+acumEx[cur];
                    acumEx[cur] = 0;
                }
            }

           }
        }
        
        
    }
    double duration = c.duration();
        rest->val = ex[sk];
        rest->flow = flows;

}