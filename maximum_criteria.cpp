/*
  * File:   main.cpp
  * Author: yuihiroaki
  *maximum criteria
  * Created on 2013/07/03, 11:25
  */
#include <iostream>
#include <list>
#include <fstream>
#include <vector>
#include "mt19937ar.h"
#include <bitset>
#include <sys/time.h>
//#include <random>

using namespace std;
using namespace mt;

class Edge {
 public:
     int src;
     int dest;
     int weight;
 };
 struct subset {
     int parent;
     int rank;
 };
 class Graph {  
  private:
      static const int NUMBER_OF_VERTICES;  
      vector<int> *adj_list_edge;
      vector<int> *wgtvtx;
      vector<int> *wgtadj;
      vector<int> grapharry;
      int E; // a number of edge
      int V; // a number of vertices   
      vector<int> *adj;      
      vector<int> longestpath;
      vector<int> *mst;    
      vector<int> bucket;  
      //vector<int> rnmvtx;
   public:
      Graph() ;   
      Graph(int Vtx, int Edg):V(Vtx),E(Edg){}
      void memory_allocate(); 
      void generateEdge(int v, int w);  
      void EdgeWeight();
      void graph_read();     
      //kruskal
      Edge* edge;
      Graph* Insert();
      void KruskalMST(Graph* graph, Graph u);
      int find(struct subset subsets[], int i);
      void Union(struct subset subsets[], int x, int y);
      vector<int> FindLongestPath(); 
      vector<int> Partitioning(vector<int> longestpath, Graph G);
      vector<int> scanU(vector<int> *mst ,vector<int> longestpath, int p);    
      vector<int> scanV1(vector<int> *mst ,vector<int> longestpath, int kmin);
 };
Graph::Graph(){}
void Graph::generateEdge(int v, int w){
     adj_list_edge[v].push_back(w);
}
void Graph::memory_allocate(){
    adj_list_edge=new vector<int>[E-V];
    wgtvtx=new vector<int>[E-V];
    wgtadj=new vector<int>[E-V];
    adj=new vector<int>[V]; 
    mst=new vector<int>[V]; 
}
double second(){
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + t.tv_usec * 1e-6;
}
void Graph::EdgeWeight(){

     const double weight_begin=second();

     vector<int> rnmvtx;
     rnmvtx.resize(V);
     for(int i=0; i < V; i++) rnmvtx[i]=i;//rnmvtx.push_back(i);    
     int mtrand;
     init_genrand((unsigned)time(NULL));
     vector<int>::iterator irp;


     for(auto irp : rnmvtx){ //shufle by Fisher-Yates
         mtrand=genrand_int31()%(rnmvtx.size());
         swap(rnmvtx[irp], rnmvtx[mtrand]);
     }

     


    vector<bool> cluster(E-V, false);
    //vector<int>::iterator adj_vtx; 
    int weight;   
    int vertex_weight;
    int bucket_index=0;
     for(auto index : rnmvtx){
        vertex_weight=adj[index].size();
        for(auto adj_vtx : adj_list_edge[index]){
                 weight=vertex_weight+adj[adj_vtx].size(); 
                 cluster[weight]=true;   
                 wgtvtx[weight].push_back(adj_vtx);
                 wgtadj[weight].push_back(index);
         }
     }

   
    for(int i=0; i<E-V; ++i){
         if(cluster[i]==true){
            bucket.push_back(i);

         }
     }


    const double weight_end=second();
    //cout << "Weight : " << weight_end - weight_begin<< endl;

 }
 const int Graph::NUMBER_OF_VERTICES(289);
 void Graph::graph_read(){
     const char *graph_name="mesh3e1.txt";
     Graph vtx_edge(
289,1089
);
     vtx_edge.memory_allocate();
     
     cout << "graph : " << graph_name << endl;
     ifstream ifs(graph_name);       
     string str;
     int v1=0;
     int v2=0;
     while(getline(ifs,str)) {
        sscanf(str.data(), "%d %d", &v1, &v2);     
        if(v1!=v2){
            vtx_edge.adj[v1].push_back(v2);
            vtx_edge.adj[v2].push_back(v1); 
        
            vtx_edge.generateEdge(v1,v2);

        }
     }  


     const double begin_time=second();
     vtx_edge.EdgeWeight();
     Graph* p_graph=vtx_edge.Insert();  
     p_graph->KruskalMST(p_graph,vtx_edge);    
     vector<int> LP=p_graph->FindLongestPath(); 
     p_graph->Partitioning(LP,vtx_edge);
     const double end_time = second();
     cout << "maximum criteria: " <<end_time- begin_time << endl;

     delete p_graph;

 }
 Graph* Graph::Insert(){

      const clock_t Ins_begin=clock();
      Graph* graph=new Graph ; 
      graph->V=V;
      graph->E=E;
      graph->edge=new Edge[graph->E] ;
      vector<int>::reverse_iterator brt;
      int weightCount=1;
      int wgtvtx_size;
      for(brt=bucket.rbegin(); brt!=bucket.rend(); ++brt){    
      //for(auto brt : bucket){
          wgtvtx_size=wgtvtx[*brt].size();
          for(int i=0; i<wgtvtx_size; ++i){
            graph->edge[weightCount-1].src=wgtvtx[*brt][i];
            graph->edge[weightCount-1].dest=wgtadj[*brt][i];
            graph->edge[weightCount-1].weight=weightCount;
            ++weightCount;
            //cout << weightCount << endl;
         }
      }
     const clock_t Ins_end=clock();
     //cout << "Insert : " << float(Ins_end - Ins_begin)/CLOCKS_PER_SEC << endl;   
     //delete graph;
     return graph;
  }
 int Graph::find(subset subsets[], int i){  
     if(subsets[i].parent!=i) subsets[i].parent=find(subsets, subsets[i].parent);
     return subsets[i].parent;
}
void Graph::Union(subset subsets[], int x, int y){
     int xroot=find(subsets, x);
     int yroot=find(subsets, y);
     if(subsets[xroot].rank<subsets[yroot].rank){
         subsets[xroot].parent=yroot;
     }else if(subsets[xroot].rank>subsets[yroot].rank){
         subsets[yroot].parent=xroot;
     }else{
         subsets[yroot].parent=xroot;
         subsets[xroot].rank++;
     }
}
void Graph::KruskalMST(Graph* graph, Graph u){    
     const clock_t mst_begin=clock();   
     int V=graph->V;
     Edge result[V];  
     int e=0; 
     int i=0;  
     subset *subsets=(subset*)malloc(V * sizeof(subset));
     for(int v=0; v<V; ++v){
         subsets[v].parent=v;
         subsets[v].rank=0;
     }

     while (e<V-1){
         Edge next_edge = graph->edge[++i];
         int x=find(subsets, next_edge.src);
         int y=find(subsets, next_edge.dest);

         if (x!=y){
             result[++e]=next_edge;
             Union(subsets, x, y);
         }
     }

     mst=new vector<int>[E];
     //cout << "Following are the edges in the constructed MST " << endl;
     for (i= 0; i<e; ++i){
       //cout << result[i].src +1 << "<->" << result[i].dest +1<< ", ";
       if(result[i].src!=result[i].dest){
           mst[result[i].src].push_back(result[i].dest);
           mst[result[i].dest].push_back(result[i].src);
       }
     }

     const clock_t mst_end = clock();
    // cout << "Kruskal : " << float(mst_end- mst_begin )/CLOCKS_PER_SEC << endl;
     return;
 }
vector<int> Graph::FindLongestPath(){   
    const clock_t lp_begin = clock();
    init_genrand((unsigned)time(NULL));
    int first_vtx=genrand_int31() % V;
    vector<bool> visited(V, false);
    
    list<int> first_path;
    visited[first_vtx]=true;
    first_path.push_back(first_vtx);
    vector<int>::iterator vtx;
    while(!first_path.empty()){
        first_vtx=first_path.front();      
        first_path.pop_front();
       for(vtx=mst[first_vtx].begin(); vtx!=mst[first_vtx].end(); ++vtx){   
            if(!visited[*vtx]){
                visited[*vtx]=true;
                first_path.push_back(*vtx);
            }
        }      
    }

    visited.assign(visited.size(), false);
    //second find the longest path from the last, 
    //vertex which we found in the first operation.
    int far_vtx=first_vtx;// the most far vertex 
    //Create a queue for BFS
    list<int> queue;
    list<int> adj_size;
    vector<int> add_adj_size;
    visited[far_vtx]=true;
    queue.push_back(far_vtx);
    vector<int> node;
    //vector<int>::iterator vtx2;
    while(!queue.empty()){
        far_vtx=queue.front();
        node.push_back(far_vtx);
        adj_size.push_back(mst[far_vtx].size()-1);
        queue.pop_front();
        for(auto vtx2 : mst[far_vtx]){//.begin(); vtx2!=mst[far_vtx].end(); ++vtx2){
            if(!visited[vtx2]){                
                visited[vtx2]=true;
                queue.push_back(vtx2);
            }
        }
    }
  //cout << endl << "show the adjacent sizes" << endl;
 // list<int>::iterator itx;
  adj_size.pop_front();
  adj_size.push_front(1);
  int sum=0;
  for(auto itx : adj_size){//.begin(); itx!=adj_size.end(); ++itx){     
        sum=sum+(itx);
        add_adj_size.push_back(sum);
  }
  int size_equence=0;
  int tmp_num=0;
  add_adj_size[0]=1;
  list<int> level;
  vector<int>::iterator itr=add_adj_size.begin(); 
  while(itr!=add_adj_size.end()){
         size_equence=size_equence+add_adj_size[size_equence]-tmp_num; 
         level.push_back(size_equence);
         tmp_num=size_equence;
         if(size_equence>=node.size()-1)break;     
         ++itr;                
  }
  visited.assign(visited.size(), false);
  level.push_front(0);
  list<int>::reverse_iterator level_itr = level.rbegin();
  vector<int>::reverse_iterator node_begin= node.rbegin(); 
  vector<int>::iterator trace_itr;
  bool bFlag=false;
  vector<int> longestpath;
  int node_back=node.back(); // making short variable for print 
  //int lp_back=longestpath.back();// making short variable for print 
  longestpath.push_back(node_back);// node is vector
  while(level_itr!=level.rend()){
     level.pop_back();
     while(node_begin!=node.rend()){ 
         if(*node_begin==node[level.back()]){
             break;
         }else{
             visited[*node_begin]=true;  
         }
         ++node_begin;      
     }
    ++level_itr;
    int lp_back=longestpath.back();// making short variable for print 
    if(bFlag==false){
       for(trace_itr = mst[node_back].begin(); trace_itr != mst[node_back].end(); ++trace_itr){      
          if(!visited[*trace_itr])longestpath.push_back(*trace_itr);
       } 
       bFlag=true;
    }else{
       for(trace_itr=mst[lp_back].begin(); trace_itr!= mst[lp_back].end(); ++trace_itr){    
          if(!visited[*trace_itr]){
              longestpath.push_back(*trace_itr);
              break;
          }              
       }
    }  
  }
  //cout << endl<< "Longest Path for mathematica vertices" << endl;
  //vector<int>::iterator pi;
  //for(pi=longestpath.begin();pi!=longestpath.end();pi++)cout<< *pi+1 << ", ";
  const clock_t lp_end = clock();
  //cout << "Longest Path  :  " << float(lp_end- lp_begin)/CLOCKS_PER_SEC << endl;
  //cout << "Longest Path " << longestpath.size() << endl;

  return longestpath;
}
vector<int> Graph::Partitioning(vector<int> longestpath, Graph G){
    const clock_t Pa_begin = clock();
    bitset<NUMBER_OF_VERTICES> bvU,bvAdjU,bvNewseparator,bvSeparator,bvAdjV1,
    bvV1,bvV2,bvV,bvW, S;
    int p1=longestpath.front();
    
    vector<int>::iterator Adjp1;
    for(Adjp1=G.adj[p1].begin(); Adjp1!=G.adj[p1].end(); Adjp1++)bvAdjU.set(*Adjp1 ,1);//Adj(p1)
    
    bvSeparator=bvAdjU; // Separator <- Adj(p1)
    //bvW.set(p1,1);      //W1<-{p1}
    bvU.set(p1,1);      //U1<- W1
    bvAdjU.reset();
    int max=0;
    int kmin;
    
    vector<int> V1;
    vector<int> U;
    //vector<int>::iterator Uvtx;
    //vector<int>::iterator AdjU;
    int SizeOfV1,SizeOfV2,square;
    int longestpath_size=longestpath.size();
    for(int p=1; p<longestpath_size; ++p){
          U =scanU(mst,longestpath,p);  
          for(auto Uvtx : U){                    
              for(auto AdjU : G.adj[Uvtx]){
                        bvAdjU.set(AdjU,1);//Adj(U)
              } 
          }  

          for( auto Uvtx : U ) bvW.set(Uvtx,1);

          bvNewseparator=(bvAdjU & ~bvW)|(bvSeparator & ~bvW);
          //bvNewseparator=(bvAdjU|bvSeparator) & ~bvW;
          //bvW=bvW | bvU;//Wk<-Wk-1+Uk
          SizeOfV1=bvW.count(); // + => |,- => & ~
          SizeOfV2=V-bvNewseparator.count()-SizeOfV1; 
          square=SizeOfV1*SizeOfV2; 
          
          if(max<square){
               max=square;              
               kmin=p;  
           }
           bvAdjU.reset();
           //bvU.reset();
           bvSeparator=bvNewseparator;
           bvNewseparator.reset(); 
           U.clear();
    }
    //cout << "kmin " << kmin << endl;
    V1 = scanV1(mst,longestpath, kmin);   
    vector<int>::iterator scanV1;  
    vector<int>::iterator adjV1;
    for(scanV1=V1.begin(); scanV1!=V1.end(); scanV1++){           
       bvV1.set(*scanV1,1); 
       for(adjV1=G.adj[*scanV1].begin(); adjV1!=G.adj[*scanV1].end(); ++adjV1){            
                  bvAdjV1.set(*adjV1,1);           
        }                     
   } 
    bvV.set();    
    S =(bvAdjV1) & (bvV & ~bvV1);
   // cout << bvNewseparator.count() << endl;
    bvV.set();
    bvV2=bvV & ~bvV1 & ~S;   
    const clock_t Pa_end = clock();
   // cout << "Partitioning : " << float(Pa_end- Pa_begin)/CLOCKS_PER_SEC << endl;
    cout << endl;
/*
    cout << "|V| " << V << ": |E| " << E << endl;
    cout <<endl << "*none-refinement process " << endl; 
    cout << "|S|  " << S.count() <<endl;
    cout << "|V1| " << bvV1.count() <<endl;
    cout << "|V2| " << bvV2.count() <<endl; 
*/
    double FullmatrixCost=bvV.count()*bvV.count();
/*
    double CostDP = FullmatrixCost-2*bvV1.count()*bvV2.count();
    cout << "ρ for none-refiment : " << CostDP/FullmatrixCost << endl;
    */
    vector<int> Separator;
    for(int i=V-1;i >= 0;i--) if(S[i]==1) Separator.push_back(i);  
    
    //cout <<endl << "*refinement process " <<endl;
    vector<int>::iterator vtxOfSep; 
    vector<int>::iterator adj_Sep;
    //vector<int> Refinement;
    bitset<NUMBER_OF_VERTICES> bvAdjSep, bvref, bvRefinement;
    
    for(vtxOfSep=Separator.begin(); vtxOfSep!=Separator.end(); vtxOfSep++){
        for(adj_Sep=G.adj[*vtxOfSep].begin(); adj_Sep!=G.adj[*vtxOfSep].end(); ++adj_Sep){
            bvAdjSep.set(*adj_Sep,1);           
        }
        bvref=(bvAdjSep & bvV2);
        if(bvref.any()){
            bvRefinement.set(*vtxOfSep);
           // Refinement.push_back(*vtxOfSep);
        }
        bvref.reset();
        bvAdjSep.reset();
    }
    bitset<NUMBER_OF_VERTICES> bvRefV1;
    bvRefV1=bvV1 | (S & ~bvRefinement);
/*
    cout << "|S'|  " << bvRefinement.count() <<endl;
    cout << "|V'1| " << bvRefV1.count() << endl;
    cout << "|V'2| " << bvV2.count() <<endl;    
*/
    double CostDP_Ref=FullmatrixCost-2*bvRefV1.count()*bvV2.count();
  //  cout <<"ρ for refinement : "  << CostDP_Ref/FullmatrixCost << endl;
    cout << endl;
  
    return grapharry;//CostDP_Ref/FullmatrixCost;
} 
vector<int> Graph::scanU(vector<int> *mst, vector<int> longestpath, int p){
    int s=longestpath[p];
    vector<bool> visited(V, false);
    if(p==longestpath.size()-1){
      visited[longestpath[p-1]]=true;
    }else{    
      visited[longestpath[p-1]]=true;
      visited[longestpath[p+1]]=true;
    }
    list<int> queue;
    visited[s]=true;
    queue.push_back(s);
    vector<int>::iterator vtx;
    vector<int> U;
    while(!queue.empty()){
        s=queue.front();
        U.push_back(s);
        queue.pop_front();
        for(vtx=mst[s].begin(); vtx!=mst[s].end(); ++vtx){  
            if(!visited[*vtx]){
                visited[*vtx]=true;
                queue.push_back(*vtx);
            }
        }
    }
    if(p==longestpath.size()-1){
      visited[longestpath[p-1]]=false;
    }else{    
      visited[longestpath[p-1]]=false;
      visited[longestpath[p+1]]=false;
    }
    return U;
}
vector<int>Graph::scanV1(vector<int> *mst, vector<int> longestpath, int kmin){

    int s=longestpath.front();
    vector<bool> visited(V, false);
    visited[longestpath[kmin+1]]=true;
    list<int> queue;
    visited[s]=true;
    queue.push_back(s);
    vector<int>::iterator vtx;
    vector<int> V1;
    while(!queue.empty()){
        s = queue.front();
        V1.push_back(s);
        queue.pop_front();
        for(vtx=mst[s].begin(); vtx!=mst[s].end(); ++vtx){
            if(!visited[*vtx]){
                visited[*vtx]=true;
                queue.push_back(*vtx);
            }
        }
    }
    visited[longestpath[kmin+1]]=false;
    return V1;
}
int main(){
     Graph g;
     g.graph_read();
     return 0;
}
