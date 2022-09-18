#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <omp.h>
#include "tbb/concurrent_vector.h"
#include "PushRelabelAlgo.h"
using namespace std;

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
   cout<<"processor close version\n";
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

   //omp_set_num_threads(8);
   //string d = "dir";
   //Problem p = readGraph(to_string(1));
   //cout<<compare(p,8,"whole","2","1")<<" finished\n";
}



