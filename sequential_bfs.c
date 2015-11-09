/* ==========  ========== ========== ========= */
//            Author - Anuj & Aniruth             //
//                 CSE 603                     //
/* ========== ========== ========== ========== */
  
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
//#include <mpi.h>
  
struct Edge {
    int vertex;
    struct Edge * next;
};
  
// Inserts Node to the Linked List by Head Insertion - O(1)
// Returns address of head which is the newly created node.
struct Edge * AddEdge(struct Edge * currentHead, int newVertex)
{
    struct Edge * newHead
                 = (struct Edge *) malloc(sizeof(struct Edge));
  
    newHead->vertex = newVertex;
    newHead->next = currentHead;
  
    return newHead;
}
  
void BreadthFirstSearch(
                        struct Edge * adjacencyList[],
                        int vertices,
                        int parent[],
                        int level[],
                        int startVertex
                       )
{
    struct Edge * traverse;
    int i, par, lev, flag = 1;
    // 'lev' represents the level to be assigned
    // 'par' represents the parent to be assigned
    // 'flag' used to indicate if graph is exhausted
  
    lev = 0;
    level[startVertex] = lev;
    // We start at startVertex
  
    while (flag) {
        flag = 0;
        for (i = 1; i <= vertices; ++i) {
            if (level[i] == lev) {
                flag = 1;
                traverse = adjacencyList[i];
                par = i;
  
                while (traverse != NULL) {
                    if (level[traverse->vertex] != -1) {
                        traverse = traverse->next;
                        continue;
                    }
  
                    level[traverse->vertex] = lev + 1;
                    parent[traverse->vertex] = par;
                    traverse = traverse->next;
                }
            }
        }
  
        ++lev;
    }
}
  
int main(int argc, char **argv)
{
    int vertices, edges, i, v1, v2;
    char ch, file_name[25];
    FILE *fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    int startingVertex = 0;

    // ------------ opening and reading input from file graphData.txt ------------
    fp = fopen("graphData.txt","r"); // read mode
 
    if( fp == NULL )
    {
       perror("Error while opening the input data file file.\n");
       exit(EXIT_FAILURE);
    }

    // Reading from file and setting values of all attributes
    int counter = 0;
    int adjCounter = 0;
    while ((read = getline(&line, &len, fp)) != -1) {
        
        int val = atoi(line);
        if(counter == 0){
          // saving number of vertices
          vertices = val;
        }
        else if (counter == 1){
          // saving number of edges
          edges = val;
        }
        else{
          // closing file and breaking loop
          startingVertex = val;
          fclose(fp);
          break;
        }

        // incrementing counter
        counter++;

    }

    counter = 0;
    adjCounter = 0;
  
    struct Edge * adjacencyList[vertices + 1];
    // Size is made (vertices + 1) to use the
    // array as 1-indexed, for simplicity
  
    int parent[vertices + 1];
    // Each element holds the Node value of its parent
    int level[vertices + 1];
    // Each element holds the Level value of that node
  
    // Initializing array
    for (i = 0; i <= vertices; ++i) {
        adjacencyList[i] = NULL;
        parent[i] = 0;
        level[i] = -1;
    }
  
    // again reading from file
    fp = fopen("graphData.txt","r"); // read mode
 
    if( fp == NULL )
    {
       perror("Error while opening the input data file file.\n");
       exit(EXIT_FAILURE);
    }

    int newEdge = 1, prevEdge = 0;
    while ((read = getline(&line, &len, fp)) != -1) {
        
        int val = atoi(line);
        if(counter == 0){
          // do nothing as number of vertices already saved
        }
        else if (counter == 1){
          // do nothing as number of edges already saved
        }
        else{
          // saving all edges
          if(newEdge == 1){
            prevEdge = val;
            newEdge = 0;
          }
          else{
            v1 = prevEdge;
            v2 = val;

            // Adding edge v1 --> v2
            adjacencyList[v1] = AddEdge(adjacencyList[v1], v2);
  
            // Adding edge v2 --> v1
            // Remove this if you want a Directed Graph
            adjacencyList[v2] = AddEdge(adjacencyList[v2], v1);

            newEdge = 1;
          }
        }

        // incrementing counter
        counter++;

    }
  
    // Printing Adjacency List
    printf("\nAdjacency List -\n\n");
    for (i = 1; i <= vertices; ++i) {
        printf("adjacencyList[%d] -> ", i);
  
        struct Edge * traverse = adjacencyList[i];
  
        while (traverse != NULL) {
            printf("%d -> ", traverse->vertex);
            traverse = traverse->next;
        }
  
        printf("NULL\n");
    }
  
    clock_t start = time(NULL), diff;
    printf("\n\nStarting vertex is: %d\n\n", startingVertex);
     
    // ---------------------Breadth first search start------------------------
    // starting mpi here (above this section is sequential code)
    //MPI_Init(&argc, &argv);

    struct Edge * traverse;
    int j, par, lev, flag = 1;
    // 'lev' represents the level to be assigned
    // 'par' represents the parent to be assigned
    // 'flag' used to indicate if graph is exhausted
  
    lev = 0;
    level[startingVertex] = lev;
    // We start at startVertex
  
    while (flag) {
        flag = 0;
        for (j = 1; j <= vertices; ++j) {
            if (level[j] == lev) {
                flag = 1;
                traverse = adjacencyList[j];
                par = j;
  
                while (traverse != NULL) {
                    if (level[traverse->vertex] != -1) {
                        traverse = traverse->next;
                        continue;
                    }
  
                    level[traverse->vertex] = lev + 1;
                    parent[traverse->vertex] = par;
                    traverse = traverse->next;
                }
            }
        }
  
        ++lev;
    }

    //MPI_Finalize();
    // finalizing and closing mpi (below this section is sequential output code
    // ---------------------Breadth first search ends-------------------------

    diff = time(NULL) - start;
  
    // Printing Level and Parent Arrays
    printf("\nLevel and Parent Arrays -\n");
    for (i = 1; i <= vertices; ++i) {
        printf("Level of Vertex %d is %d\n", i, level[i]);
    }

    /* Print out the difference */
    printf ( "\nTime consumed was: %.2f \n\n", (double)diff);
    
    return 0;
}
