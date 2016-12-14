#include <Rcpp.h>
using namespace Rcpp;

// This class represents an undirected graph using adjacency list representation
class Graph
{
  int V;    // No. of vertices
  std::list<int> *adj;    // Pointer to an array containing adjacency lists
public:
  Graph(int V);  // Constructor
  void addEdge(int v, int w); // function to add an edge to graph
  void init_w_adjlist(const List & adjlist, int n); // add a whole adjacency list
  IntegerVector getEdges(int v); // function to add an edge to graph
  void BFS(int s);  // prints BFS traversal from a given source s
};

Graph::Graph(int V)
{
  this->V = V;
  adj = new std::list<int>[V];
}

void Graph::addEdge(int v, int w)
{
  if (std::find(adj[v].begin(), adj[v].end(), w) == adj[v].end()) // w not in adj[v]
    adj[v].push_back(w); // Add w to v’s list.
  if (std::find(adj[w].begin(), adj[w].end(), v) == adj[w].end()) // v not in adj[w]
    adj[w].push_back(v); // Add v to w’s list.
}
void Graph::init_w_adjlist(const List & adjlist, int n)
{
  for (int i=0; i<n; i++) {
    int ni = as<IntegerVector>(adjlist[i]).size();
    for (int j=0; j<ni; j++) 
      addEdge(i, as<IntegerVector>(adjlist[i])[j] - 1);
  }
}

IntegerVector Graph::getEdges(int v)
{
  IntegerVector rval = wrap(adj[v]);
  return rval;
}


void Graph::BFS(int s)
{
  // Mark all the vertices as not visited
  bool *visited = new bool[V];
  for(int i = 0; i < V; i++)
    visited[i] = false;
  
  // Create a queue for BFS
  std::list<int> queue;
  
  // Mark the current node as visited and enqueue it
  visited[s] = true;
  queue.push_back(s);
  
  // 'i' will be used to get all adjacent vertices of a vertex
  std::list<int>::iterator i;
  
  while(!queue.empty())
  {
    // Dequeue a vertex from queue and print it
    s = queue.front();
    std::cout << s << " ";
    queue.pop_front();
    
    // Get all adjacent vertices of the dequeued vertex s
    // If a adjacent has not been visited, then mark it visited
    // and enqueue it
    for(i = adj[s].begin(); i != adj[s].end(); ++i)
    {
      if(!visited[*i])
      {
        visited[*i] = true;
        queue.push_back(*i);
      }
    }
  }
}

// program to test methods of graph class
// [[Rcpp::export]]
void test(const List & adjlist)
{
  // Create a graph 
  int n = adjlist.size();
  Graph g(n);
  g.init_w_adjlist(adjlist, n);
  
  IntegerVector ll = g.getEdges(2);
  print(ll);
}

/*** R
testthat::test_that("adjgraph gives no error", {
  labelmat = matrix(1:100, 10, 10)
  B = Matrix::bandSparse(10, k = c(-1, 1, -1)); I = diag(10)
  P = kronecker(B, I) + kronecker(I, B)
  adjlist = sapply(1:nrow(P), function(i) which(P[i,] == 1))
  test(adjlist)
})
*/