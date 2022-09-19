#include <iostream>
#include <fstream>
#include <omp.h>
#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_unordered_set.h"
#include <atomic>
#include <vector>
#include <unordered_set>
#include <chrono>
#include <string>
using namespace std;
//Usage Example
//compile with the command : g++ -g *.cpp -o test -fopenmp -I/usr/include/tbb -ltbb
//If using spartan, please load module load gcccore/10.2.0 module load tbb/2020.3  first
//To use the code, run command like following
//  ./test single dir 1 8     
//test is the compiled code, single means just run parallelized max flow algorithm, dir is the director name, which must under the current working director
// The 1 is the name of the graph file, the name of the graph file must be interger. The following 8 means how many threads you want to set
//You could use the 1,3,5 file test, which satisfies the graph format requirement. The format specification could be seen in dir/info.txt
//The file 5 has large node numbers, please pay attention
//If you want to use the cache unfriendly version, just replace the single with cache like ./test single dir 1 8  , which is used for specific experiment
//To run the whole algorithm, find the minimum of all the max flow, run following command
//./test dir 4 whole 8 8 2
// test is the compiled code, dir is the directory name, which must under the current working directory,
// 4 is a tag number ,must be int, used in experiment, no specific meaning
//whole is also label, could be string, use to indicate test content
//The following two 8, should all be int
//The first is the total core num you plan to use
//the second is the core you want to use for the outer part parallelization. (Both value are suggested to be the same when you are using small graph)
//the reason is that in small graph, if more threads are allocated to the inner parallelization part, threads are idle
//Difference could be seen from the report.
//The required graph format is described in info.txt under dir directory
//One thing is that all the graph should have no anti-parallel edges.
//If you want to generate more graphs, please seek https://github.com/xhnhill/MaxFlowGraphGenerator.git
//I have developed the above java code for generation
//One more thing you need attention, this code's output (including system.out and files generated contains the execution time.)
//If you want test it on your own data, please define files like the required input, which is defined in dir directory with info.txt
//If you want to test the minum of all max flow,and see if your result is right, use the above required graph format, if it is wrong, the program will output wrong match on output
//If you just want to test the max flow problem, please use the required graph format, and run with the command with single parameter.
//Usually, in the test of the report, the above output code will be commented to save running time.

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
  int totalNum;
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
     double pushRelabelWrapper(MaxFlowStd *mx,MaxFlowRes *result,int coreNum);
     double cachePushRelabel(MaxFlowStd *mx,MaxFlowRes *result);
     double cachePushRelabelWrapper(MaxFlowStd *mx,MaxFlowRes *result,int coreNum);
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
    void cacheGlobalRelabel(int numV, int sc, int sk);
    void tryGlobalRelabel(MaxFlowStd *mx, MaxFlowStd *res);
    

};

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
        
        for(int i:queue){
            #pragma omp parallel for
            for(int j =0;j<discoveredVertices[i].size();j++){
                tmp.push_back(discoveredVertices[i][j]);
                //cout<<"push "<<discoveredVertices[cur][j]<<"\n";
            }
        }
        queue.swap(tmp);

    }
}
void PushRelabelAlgo::cacheGlobalRelabel(int nv, int sc, int sk){
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
            for(int j =0;j<nv;j++){
                if(i!=j && res[j][i]>0){
                    if(j==sc) continue;
                    int rpVal = nv;
                    //Label possible edge on res network with distance to sink
                    if(j != sk && h[j].compare_exchange_strong(rpVal,h[cur]+1)){
                        discoveredVertices[cur].push_back(j);
                        
                    }
                }
            }
            
        }
        tbb::concurrent_vector<int> tmp;
        //Difference??
        
        for(int i:queue){
            #pragma omp parallel for
            for(int j =0;j<discoveredVertices[i].size();j++){
                tmp.push_back(discoveredVertices[i][j]);
                //cout<<"push "<<discoveredVertices[cur][j]<<"\n";
            }
        }
        queue.swap(tmp);

    }
}
double PushRelabelAlgo::pushRelabel(MaxFlowStd* mx,MaxFlowRes *rest){
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
                    //use += is atomic on atomic variable, or using explict atomic expression!
                    acumEx[v] +=(le - ex[v]);
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
            //atomic attention
            ex[i] +=acumEx[i];
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
                    ex[cur] +=acumEx[cur];
                    acumEx[cur] = 0;
                }
            }

           }
        }
        
        
    }
    double duration = c.duration();
        rest->val = ex[sk];
        //rest->flow = flows;
        //free mem
        free(ex);
        free(acumEx);
        free(localEx);
        free(localLabel);
        free(h);
        for(int i=0;i<nv;i++){
            free(res[i]);
            free(flows[i]);
            reverseMap->clear();
        }
        free(res);
        free(flows);
        free(work);
        activeSet.clear();
        return duration;

}
double PushRelabelAlgo::cachePushRelabel(MaxFlowStd *mx,MaxFlowRes *rest){
    //begin record
    double duration =0;
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
            c.reset();
            cacheGlobalRelabel(nv, sc, sk);
             duration= duration+c.duration();
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
                    //use += is atomic on atomic variable, or using explict atomic expression!
                    acumEx[v] +=(le - ex[v]);
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
            //atomic attention
            ex[i] +=acumEx[i];
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
                    ex[cur] +=acumEx[cur];
                    acumEx[cur] = 0;
                }
            }

           }
        }
        
        
    }
    
        rest->val = ex[sk];
        //rest->flow = flows;
        //free mem
        free(ex);
        free(acumEx);
        free(localEx);
        free(localLabel);
        free(h);
        for(int i=0;i<nv;i++){
            free(res[i]);
            free(flows[i]);
            reverseMap->clear();
        }
        free(res);
        free(flows);
        free(work);
        activeSet.clear();
        return duration;
}
double PushRelabelAlgo::pushRelabelWrapper(MaxFlowStd *mx,MaxFlowRes *result,int coreNum){
    Clock c;
    c.reset();
    omp_set_num_threads(coreNum);
    pushRelabel(mx,result);
    return c.duration();
}
double PushRelabelAlgo::cachePushRelabelWrapper(MaxFlowStd *mx,MaxFlowRes *result,int coreNum){
    
    
    omp_set_num_threads(coreNum);
    return cachePushRelabel(mx,result);
    
}


//Execution part

Problem readGraph(string path){
   ifstream inStream (path);
   int fullNum; inStream>>fullNum;
   int edgePairNum; inStream>>edgePairNum;
   int actualNodeNum; inStream>>actualNodeNum;
   int res; inStream>>res;
   int** cap = new int*[fullNum];
   for(int i =0;i<fullNum;i++){
      cap[i] = new int[fullNum];
   }
   for(int i = 0;i<edgePairNum;i++){
      int sc;inStream>>sc;
      int tar;inStream>>tar;
      inStream >> cap[sc-1][tar-1];
   }
   int** compareMatrix = new int*[actualNodeNum];
   for(int i=0;i<actualNodeNum;i++){
      compareMatrix[i] = new int[actualNodeNum];
   }
   for(int i =0;i<actualNodeNum;i++){
      for(int j = 0;j<actualNodeNum;j++){
         inStream>>compareMatrix[i][j];
      }
   }
   int allMin; inStream>>allMin;
   Graph g;
   g.num_vertices = fullNum;
   g.capacities = cap;
   Problem p; p.compareMatrix=compareMatrix;p.res = allMin;p.actualNum=actualNodeNum;p.g=g;
   return p;
   
   
}
Problem readSingleProblem(string path){
   ifstream inStream (path);
   int fullNum; inStream>>fullNum;
   int edgePairNum; inStream>>edgePairNum;
   int actualNodeNum; inStream>>actualNodeNum;
   int res; inStream>>res;
   int** cap = new int*[fullNum];
   for(int i =0;i<fullNum;i++){
      cap[i] = new int[fullNum];
   }
   for(int i = 0;i<edgePairNum;i++){
      int sc;inStream>>sc;
      int tar;inStream>>tar;
      inStream >> cap[sc-1][tar-1];
   }
   
   int allMin; inStream>>allMin;
   Graph g;
   g.num_vertices = fullNum;
   g.capacities = cap;
   Problem p; p.res = allMin;p.actualNum=actualNodeNum;p.g=g;p.totalNum=fullNum;
   return p;
}
bool noWriteCompare(Problem p){
   MaxFlowStd mf;
   mf.g = p.g;
   bool flg = true;
   double avg =0;
   double count =0;
   for(int i =0;i<p.actualNum;i++){
      for(int j = 0;j<p.actualNum;j++){
         if(i!=j){
            mf.sc = i;
            mf.sk =j;
            MaxFlowRes ms;
            PushRelabelAlgo ps;
            double dt = ps.pushRelabel(&mf,&ms);
            if(ms.val !=p.compareMatrix[i][j]){
               cout<<"src-sink "<<""+to_string(i)+" "+to_string(j)<<"does not match";
               flg = false;
            }

         }
      }
   }
   
   return flg;
}
bool compare(Problem p,int corenum,string target,string dir,string idx){
   MaxFlowStd mf;
   mf.g = p.g;
   bool flg = true;
   ofstream exp;
   exp.open(dir+"/res"+to_string(corenum)+target+dir+"_"+idx);
   double avg =0;
   double count =0;
   for(int i =0;i<p.actualNum;i++){
      for(int j = 0;j<p.actualNum;j++){
         if(i!=j){
            mf.sc = i;
            mf.sk =j;
            MaxFlowRes ms;
            PushRelabelAlgo ps;
            double dt = ps.pushRelabel(&mf,&ms);
            count = count+1;
            avg=(avg+dt) / count;
            
            exp<<dt;
            exp<<" ";
            cout<<dt<<" \n";
            if(ms.val !=p.compareMatrix[i][j]){
               cout<<"src-sink "<<""+to_string(i)+" "+to_string(j)<<"does not match";
               flg = false;
            }

         }
      }
   }
   cout<<avg<<" ";
   exp.close();
   return flg;

}
bool parallelCompare(Problem p,int corenum,string target,string dir,string idx,int outNum){
   
   bool flg = true;
   ofstream exp;
   exp.open(dir+"/res"+to_string(outNum)+"_"+to_string(corenum)+target+dir+"_"+idx);
   Clock c;
   c.reset();
   #pragma omp parallel for collapse(2) num_threads(outNum)
   for(int i =0;i<p.actualNum;i++){
      for(int j = 0;j<p.actualNum;j++){
         if(i!=j){
            MaxFlowStd mf;
            mf.g = p.g;
            mf.sc = i;
            mf.sk =j;
            MaxFlowRes ms;
            PushRelabelAlgo ps;
            double dt = ps.pushRelabelWrapper(&mf,&ms,max(1,corenum/outNum));
            if(ms.val !=p.compareMatrix[i][j]){
               cout<<"src-sink "<<""+to_string(i)+" "+to_string(j)<<"does not match";
               flg = false;
            }

         }
      }
   }
   double duration = c.duration();
   exp<<duration;
   cout<<duration<<" ";
   exp.close();
   return flg;

}
void calMaxFlowParallel(Problem p,int coreNum,string dir,int outNum,int idx){
   Clock c;
   c.reset();
            MaxFlowStd mf;
            mf.g = p.g;
            mf.sc = 0;
            mf.sk =p.actualNum-1;
            MaxFlowRes ms;
            PushRelabelAlgo ps;
            double dt = ps.pushRelabelWrapper(&mf,&ms,coreNum);
            if(ms.val !=p.res){
               cout<<"does not match\n";
               
            }
            cout<<to_string(ms.val)+"\n";
            ofstream exp;
            double duration = c.duration();;
            exp.open(dir+"/res"+to_string(outNum)+"_"+to_string(coreNum)+"single"+dir+"_"+to_string(idx));
            exp<<duration;
            cout<<"dir is "+dir<<"the inner time is "<<dt;
            exp.close();
            
            

}
void calCacheMaxFlowParallel(Problem p,int coreNum,string dir,int outNum,int idx){
   Clock c;
   c.reset();
            MaxFlowStd mf;
            mf.g = p.g;
            mf.sc = 0;
            mf.sk =p.actualNum-1;
            MaxFlowRes ms;
            PushRelabelAlgo ps;
            double dt = ps.cachePushRelabelWrapper(&mf,&ms,coreNum);
            if(ms.val !=p.res){
               cout<<"does not match\n";
               
            }
            ofstream exp;
            double duration = c.duration();
            exp.open(dir+"/res"+to_string(outNum)+"_"+to_string(coreNum)+"cache"+dir+"_"+to_string(idx));
            exp<<duration;
            cout<<"dir is "+dir<<"the inner time is "<<dt;
            exp.close();
            
            

}
void testCode(){
   tbb::concurrent_vector<int> q;
   Graph g;
   g.num_vertices = 4;
   int** cap = new int*[g.num_vertices];
   for(int i =0;i<g.num_vertices;i++){
      cap[i] = new int[g.num_vertices];
   }
   cap[0][2] = 1;
   cap[0][1] = 1;
   cap[1][2] =1;
   g.capacities = cap;
   MaxFlowStd mf;
   mf.g= g;
   mf.sk=2;
   mf.sc =0;
   PushRelabelAlgo sr;
   MaxFlowRes ms;
   sr.pushRelabel(&mf,&ms);
   cout<< ms.val;
   readGraph("tmp.txt");
}

int main(int arg,char** argv)
{
   //testCode();
   if(arg>1){
      string choice(argv[1]);
      if(choice == "single"){
         string dir(argv[2]);
         int fidx = stoi(argv[3]);
         int coreNum = stoi(argv[4]);
         Problem p =readSingleProblem(dir+"/"+to_string(fidx));
            calMaxFlowParallel(p,coreNum,dir,1,fidx);
            
            for(int i=0;i<p.g.num_vertices;i++){
            free(p.g.capacities[i]);
         }
         free(p.g.capacities);
         return 0;
      }else if(choice == "cache"){
         string dir(argv[2]);
         int fidx = stoi(argv[3]);
         int coreNum = stoi(argv[4]);
         Problem p =readSingleProblem(dir+"/"+to_string(fidx));
            calCacheMaxFlowParallel(p,coreNum,dir,1,fidx);
            //cout<<p.res<<"\n";
            for(int i=0;i<p.g.num_vertices;i++){
            free(p.g.capacities[i]);
         }
         free(p.g.capacities);
         return 0;
      }
   }
   if(arg == 5){
      
      int num = stoi(argv[2]);
      int coreNum = stoi(argv[4]);
      string dir(argv[1]);
      string tar(argv[3]);
      cout<<"correct param\n";
      for(int i=1;i<=num;i++){
         omp_set_num_threads(coreNum);
         Problem p = readGraph(dir+"/"+to_string(i));
         compare(p,coreNum,tar,dir,to_string(i));
         //cout<<p.res<<"\n";
         for(int i=0;i<p.g.num_vertices;i++){
            free(p.g.capacities[i]);
         }
         free(p.g.capacities);
      }
   }else if(arg == 7){
      //This para should be set to achieve performance
      omp_set_nested(1);
      int num = stoi(argv[2]);
      int coreNum = stoi(argv[4]);
      int outNum = stoi(argv[5]);
      int fileIdx = stoi(argv[6]);
      string dir(argv[1]);
      string tar(argv[3]);
      Problem p = readGraph(dir+"/"+to_string(fileIdx));
         parallelCompare(p,coreNum,tar,dir,to_string(fileIdx),outNum);
         for(int i=0;i<p.g.num_vertices;i++){
            free(p.g.capacities[i]);
         }
         free(p.g.capacities);
   }else if(arg == 2){
      //Here we test functions of openmp
      //Only for test
      std::atomic<int> x(0);
      omp_set_nested(1);
      omp_set_num_threads(2);
      #pragma omp parallel
      {
         omp_set_num_threads(4);
         #pragma omp parallel
         {
            x++;
         }
      }
      cout<<x<<"\n";
   }
   else{
      //using default
      omp_set_num_threads(8);
      //Problem p = readGraph("default");
      //noWriteCompare(p);
      cout<<"finished\n";

   }

}
