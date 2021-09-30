# Graph Quick Tutorial

## 1. DFS
   - DFS starts traversing nodes from the root node to the farthest of the nodes (i.e. depth wise).
   - DFS generally uses stack for implementation.
   - Time Complexity - O(V + E) (in general graphs) or O(max(V, E) = O(E) in dense graph.
   - Space Complexity - O(V) or O(b.d) where b is the branching factor and d is distance from root node i.e. as we cover depth first, hence space is consumed less as compared to BFS.

####Applications

- Topological Sorting
- To test if graph is bipartite
- Solving puzzles with only one solution such as mazes


####Code

- ####DFS With Recursion
```
private List<Integer> adjVertices = new ArrayList<>();

public void dfs(int start) {
     boolean[] isVisited = new boolean[adjVertices.size()];
     dfsRecursive(start, isVisited);
}

private void dfsRecursive(int current, boolean[] isVisited) {
     isVisited[current] = true;
     visit(current);
     for (int dest : adjVertices.get(current)) {
     if (!isVisited[dest])
     dfsRecursive(dest, isVisited);
     }
}
```

- ####DFS Without Recursion
```
private List<Integer> adjVertices = new ArrayList<>();

public void dfsWithoutRecursion(int start) {
    Stack<Integer> stack = new Stack<Integer>();
    boolean[] isVisited = new boolean[adjVertices.size()];
    stack.push(start);
    while (!stack.isEmpty()) {
        int current = stack.pop();
        if(!isVisited[current]){
            isVisited[current] = true;
            visit(current);
            for (int dest : adjVertices.get(current)) {
                if (!isVisited[dest])
                    stack.push(dest);
            }
    }
}
```

## 2. BFS
   - BFS starts traversing nodes from the root node and visits all other nodes in a level by level manner (i.e. visiting the ones closest to the root first).
   - BFS generally uses queue for implementation.
   - Time Complexity - O(V + E) (in general graphs) or O(max(V, E) = O(E) in dense graph.
   - Space Complexity - O(V) or O(b^d) where b is the branching factor and d is distance from root node i.e. as we move level
   by level deeper, space grows exponentially.
   - It is optimal algorithm to generally find the shortest path problems in unweighted graph.

####Applications

- Shortest Path problems for unweighted graph (Dijkstra's Shortest Path)
- Minimum Spanning Tree for unweighted graph (Prim's MST)
- Peer to Peer networks
- Garbage Collection
- GPS Navigation Syste,

####Code

```
private List<Integer> adjVertices = new ArrayList<>();

public void bfs(int start) {
   // ArrayDeque<Integer> queue = new ArrayDeque<>();
   LinkedList<Integer> queue = new LinkedList<>();
   queue.addLast(src);
   boolean[] visited = new boolean[adjVertices.size()];
   while(queue.size() > 0){
      int rem = queue.removeFirst();

      if(visited[rem] == true){
         continue;
      }
      visited[rem] = true;
      
      for (int dest : graph.get(rem)) {
         if (visited[dest] == false) {
            queue.add(dest);
         }
      }
   }
}
```

## 3. Topological Sorting

- It is always applicable on directed acyclic graph(DAG)
- It is a linear ordering of vertices such that every directed graph edge u -> v, vertex u comes before v in the ordering.
- Time Complexity - O(V + E)
- Space Complexity - O(V)
####Code

```
public static void main(String[] args) throws Exception {
      BufferedReader br = new BufferedReader(new InputStreamReader(System.in));

      int vtces = Integer.parseInt(br.readLine());
      ArrayList<Edge>[] graph = new ArrayList[vtces];
      for (int i = 0; i < vtces; i++) {
         graph[i] = new ArrayList<>();
      }

      int edges = Integer.parseInt(br.readLine());
      for (int i = 0; i < edges; i++) {
         String[] parts = br.readLine().split(" ");
         int v1 = Integer.parseInt(parts[0]);
         int v2 = Integer.parseInt(parts[1]);
         graph[v1].add(new Edge(v1, v2));
      }

      boolean[] visited = new boolean[vtces];
      Stack<Integer> st = new Stack<>();
      for(int v = 0; v < vtces; v++){
         if(visited[v] == false){
            topological(graph, v, visited, st);
         }
      }

      while(st.size() > 0){
         System.out.println(st.pop());
      }
   }

   public static void topological(ArrayList<Edge>[] graph, int src, boolean[] visited, Stack<Integer> st) {
      visited[src] = true;
      for (Edge e : graph[src]) {
         if (!visited[e.nbr]) {
            topological(graph, e.nbr, visited, st);
         }
      }
      st.push(src);
   }
```

## 4. Kahn's Algorithm
- This is another way to do topological sorting based on the fact that a DAG graph will have at least one vertex with in-degree 0 and one vertex with out-degree 0.
- Time Complexity - O(V + E)
- Space Complexity - O(V)

####Algorithm

1. Compute in-degree for each of the vertex and initialize the count of visited nodes to 0.
2. Pick all the vertices with in-degree 0 and add them to queue.
3. Remove a vertex from queue and start processing.
   - Increment count of visited nodes by 1.
   - Decrease in-degree by 1 for all neighbouring nodes.
   - If in-degree of a neighbouring nodes is reduced to zero, then add it to the queue.
4. Repeat Step 3 until the queue is empty.
5. If count of visited nodes is not equal to the number of nodes in the graph then the topological sort is not possible for the given graph.

####Applicaitons

- It can quickly detect cycles in a graph.
   
####Code

```
public static void main(String[] args) throws NumberFormatException, IOException {
     BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
     String[] st = br.readLine().split(" ");
     int n = Integer.parseInt(st[0]);
     int m = Integer.parseInt(st[1]);

     int[][] prerequisites = new int[m][2];
     for (int i = 0; i < m; i++) {
         st = br.readLine().split(" ");
         prerequisites[i][0] = Integer.parseInt(st[0]);
         prerequisites[i][1] = Integer.parseInt(st[1]);
     }
     ArrayList<ArrayList<Integer>> graph = new ArrayList<>();
     for (int i = 0; i < n; i++) {
         graph.add(new ArrayList<>());
     }

     for (int i = 0; i < prerequisites.length; i++) {
         int u = prerequisites[i][0];
         int v = prerequisites[i][1];

         graph.get(v).add(u);
     }

     int[] ans = findOrder(n, graph);

     for (int val : ans) {
         System.out.print(val + " ");
     }
 }

 public static int[] findOrder(int numCourses, ArrayList<ArrayList<Integer>> graph) {

     int[] ans = new int[numCourses];
     int[] indegree = new int[numCourses];

     for (int i = 0; i < numCourses; i++) {
         for (int nbrs : graph.get(i)) {
             indegree[nbrs]++;
         }
     }

     LinkedList<Integer> queue = new LinkedList<>();
     for (int i = 0; i < indegree.length; i++) {
         if (indegree[i] == 0) {
             queue.addLast(i);
         }
     }

     int idx = 0;
     while (queue.size() > 0) {
         int rem = queue.removeFirst();
         ans[idx] = rem;
         idx++;

         for (int nbrs : graph.get(rem)) {
             indegree[nbrs]--;
             if (indegree[nbrs] == 0) {
                 queue.addLast(nbrs);
             }
         }
     }

     if (idx == numCourses) {
         return ans;
     } else {
         return new int[] { -1 };
     }
 }
```

## 5. Dijkstra Algorith
- To find shortest path in weighted graph having positive weights.
- It can only work with graphs having positive weights on the edges. For negative weights, this could give wrong value as it could be a possibility to traverse the same vertex again and again due to negative weight on the edge.
- Time Complexity - O(E.log(V))

```
static class Edge {
   int src;
   int nbr;
   int wt;

   Edge(int src, int nbr, int wt) {
      this.src = src;
      this.nbr = nbr;
      this.wt = wt;
   }
}

public static void main(String[] args) throws Exception {
   BufferedReader br = new BufferedReader(new InputStreamReader(System.in));

   int vtces = Integer.parseInt(br.readLine());
   ArrayList<Edge>[] graph = new ArrayList[vtces];
   for (int i = 0; i < vtces; i++) {
      graph[i] = new ArrayList<>();
   }

   int edges = Integer.parseInt(br.readLine());
   for (int i = 0; i < edges; i++) {
      String[] parts = br.readLine().split(" ");
      int v1 = Integer.parseInt(parts[0]);
      int v2 = Integer.parseInt(parts[1]);
      int wt = Integer.parseInt(parts[2]);
      graph[v1].add(new Edge(v1, v2, wt));
      graph[v2].add(new Edge(v2, v1, wt));
   }

   int src = Integer.parseInt(br.readLine());

   PriorityQueue<Pair> queue = new PriorityQueue<>();
   queue.add(new Pair(src, src + "", 0));
   boolean[] visited = new boolean[vtces];
   while(queue.size() > 0){
      Pair rem = queue.remove();

      if(visited[rem.v] == true){
         continue;
      }
      visited[rem.v] = true;
      System.out.println(rem.v + " via " + rem.psf + " @ " + rem.wsf);
      
      for (Edge e : graph[rem.v]) {
         if (visited[e.nbr] == false) {
            queue.add(new Pair(e.nbr, rem.psf + e.nbr, rem.wsf + e.wt));
         }
      }
   }
}

static class Pair implements Comparable<Pair> {
   int v;
   String psf;
   int wsf;

   Pair(int v, String psf, int wsf){
      this.v = v;
      this.psf = psf;
      this.wsf = wsf;
   }

   public int compareTo(Pair o){
      return this.wsf - o.wsf;
   }
}

```


## 6. Prim's Algorithm
- To find Minimum Spanning Tree from a graph.
- It starts with a single node and explore all adjacent nodes with all the connecting edges at every step. 
- The edges with the minimal weight causing no cycle in the graph got selected.
- It is different from Dijkstra's in implementation in the fact that while adding neighbours, the weight of vertex is not added to the neighbouring vertices.

```
static class Edge {
   int src;
   int nbr;
   int wt;

   Edge(int src, int nbr, int wt) {
      this.src = src;
      this.nbr = nbr;
      this.wt = wt;
   }
}

public static void main(String[] args) throws Exception {
   BufferedReader br = new BufferedReader(new InputStreamReader(System.in));

   int vtces = Integer.parseInt(br.readLine());
   ArrayList<Edge>[] graph = new ArrayList[vtces];
   for (int i = 0; i < vtces; i++) {
      graph[i] = new ArrayList<>();
   }

   int edges = Integer.parseInt(br.readLine());
   for (int i = 0; i < edges; i++) {
      String[] parts = br.readLine().split(" ");
      int v1 = Integer.parseInt(parts[0]);
      int v2 = Integer.parseInt(parts[1]);
      int wt = Integer.parseInt(parts[2]);
      graph[v1].add(new Edge(v1, v2, wt));
      graph[v2].add(new Edge(v2, v1, wt));
   }

   int src = 0;
   PriorityQueue<Pair> queue = new PriorityQueue<>();
   queue.add(new Pair(src, -1, 0));
   Integer[] visited = new Integer[vtces];
   while(queue.size() > 0){
      Pair rem = queue.remove();

      if(visited[rem.v] != null){
         continue;
      }
      visited[rem.v] = rem.p;
      if(rem.p != -1){
         System.out.println("[" + rem.v + "-" + rem.p + "@" + rem.wt + "]");
      }
      
      for (Edge e : graph[rem.v]) {
         if (visited[e.nbr] == null) {
            queue.add(new Pair(e.nbr, rem.v, e.wt));
         }
      }
   }
}

static class Pair implements Comparable<Pair> {
   // vertex, path, weigth
   Integer v;
   Integer p;
   int wt;

   Pair(Integer v, Integer p, int wt){
      this.v = v;
      this.p = p;
      this.wt = wt;
   }

   public int compareTo(Pair o){
      return this.wt - o.wt;
   }
}
```

## 7. Bellman Ford
- To find shortest path in directed weighted graph having negative weights.
- It doesn't work with undirected graphs with negative edges as it wil declared a negative cycle.
- The shortest path between 2 vertices will get calculated  at max in (number of max edges - 1) iteration.
- Time Complexity - O(E.log(V))

####Algorithm
1. Initialize the shortest path array with Integer.MAX_VALUE (infinity)
2. Initialize the first vertex with zero in this array
3. Iterate vertices - 1 times.
4. For each edge, check if path[u] + wt < path[v]. If it holds true, update path[v], else continue.
5. After iteration is completed, array holds the shortest path from starting vertex (could be 0) to each vertex.
6. There could be case where in path Integer.MAX_VALUE is still set. This could happen in case of disconnected graph.

####Code
```
int[] path = new int[vtx];
Arrays.fill(path, Integer.MAX_VALUE);
path[0] = 0;

for (int i = 0; i < vtx - 1; i++) {
   for (int j = 0; j < edges.length; j++) {
      int u = edges[j][0];
      int v = edges[j][1];
      int wt = edges[j][2];
      
      if (path[u] == Integer.MAX_VALUE) {
         continue;
      }
      
      if (path[u] + wt < path[v]) {
         path[v] = path[u] + wt;
      }
   }
}
```

## 8. Kosaraju Algorithm
- To find strongly connected components (SCC)

####Algorithm
1. Do random order DFS and store it in stack while backtracking.
2. Reverse the edges of all vertices
3. Do dfs on the vertices mentioned in the stack
4. It will give out the number of strongly connected components

####Code

```
public static int kosaraju(int vtx, ArrayList<ArrayList<Integer>> graph) {
   boolean[] visited = new boolean[vtx];
   LinkedList<Integer> stack = new LinkedList<>();
   
   // Step 1
   for (int i = 0; i < vtx; i++) {
      if (visited[i] == false) {
         dfs(i, graph, stack, visited);
      }
   }
   
   // Step 2
   ArrayList<ArrayList<Integer>> newGraph = new ArrayList<>();
   for (int i = 0; i < vtx; i++) {
      newGraph.add(new ArrayList<>());
   }
   for (int i = 0; i < vtx; i++) {
      ArrayList<Integer> nbrs = graph.get(i);
      for (int nbr: nbrs) {
         newGraph.get(nbr).add(i);
      }
   }
   
   // Step 3
   int ans = 0;
   visited = new boolean[vtx];
   while (stack.size() > 0) {
      int rem = stack.removeFirst();
      if (visited[rem] == false) {
         dfs2(rem, newGraph, vis);
         ans++;
      }
   }
   
   return ans;
}

public static void dfs(int src, ArrayList<ArrayList<Integer>> graph, LinkedList<Integer> stack, boolean[] visited) {
   visited[src] = true;
   
   ArrayList<Integer> nbrs = graph.get(src);
   
   for (int nbr: nbrs) {
      if (visited[nbr] == false) {
         dfs(nbr, graph, stack, visited);
      }
   }
   stack.addFirst(src);
}

public static void dfs2(int src, ArrayList<ArrayList<Integer>> graph, boolean[] visited) {
   visited[src] = true;
   
   ArrayList<Integer> nbrs = graph.get(src);
   
   for (int nbr: nbrs) {
      if (visited[nbr] == false) {
         dfs2(nbr, graph, visited);
      }
   }
}

```

## 9. Disjoint Set Union (Union Find)
- It is used in applications where we need to group components based on similarity and there is a transitive dependency as well. If there is no similarity in transitive, then union find is of no use.
- In union find algorithm, find will return the leader and union will merge the leader one into another.
- It is usually used in graph problems but this can be used in other applications too.
- Time Complexity - O(n) & Space Complexity - O(n)
- To optimize this, we need to apply following algorithms
   - Path Compression
   - Union by Rank
- Post optimization, Time Complexity - O(1) & Space Complexity - O(1)

####Code

```
public int find(int x, int[] parent) {
   if (parent[x] == x) {
      return x;
   }
   int temp = find(x, parent);
   
   parent[x] = temp; // Path compression
   
   return temp;
}

public void union(int x, int y, int[] parent, int[] rank) {
   int leaderX = find(x, parent);
   int leaderY = find(y, parent);
   if (leaderX != leaderY) {
      
      // Rank by Union
      
      if (rank[leaderX] > rank[leaderY]) {
         parent[leaderY] = leaderX;
      } else if (rank[leaderX] < rank[leaderY]) {
         parent[leaderX] = leaderY;
      } else {
         parent[leaderY] = leaderX;
         rank[leaderX]++;
      }
   }
}
```

## 10. Kruskal's Algorithm
   To find Minimum Spanning Tree using Disjoint Set Union (Union Find)

## 11. Eulerian Path