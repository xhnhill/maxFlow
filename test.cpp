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
   if(arg == 5){
      
      int num = stoi(argv[2]);
      int coreNum = stoi(argv[4]);
      string dir(argv[1]);
      string tar(argv[3]);
      cout<<"correct param\n";
      for(int i=1;i<=num;i++){
         omp_set_num_threads(8);
         Problem p = readGraph(dir+"/"+to_string(i));
         compare(p,coreNum,tar,dir,to_string(i));
         for(int i=0;i<p.g.num_vertices;i++){
            free(p.g.capacities[i]);
         }
         free(p.g.capacities);
      }
   }else{
      //using default
      omp_set_num_threads(8);
      Problem p = readGraph("default");
      noWriteCompare(p);
      cout<<"finished\n";

   }

   //omp_set_num_threads(8);
   //string d = "dir";
   //Problem p = readGraph(to_string(1));
   //cout<<compare(p,8,"whole","2","1")<<" finished\n";
}



