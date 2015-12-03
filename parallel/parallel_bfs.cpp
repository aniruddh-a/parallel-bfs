#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/adjacency_iterator.hpp>
#include <boost/graph/graph_traits.hpp>
#include <vector>
#include <mpi.h>
#include <boost/lockfree/queue.hpp>
#include <boost/graph/graphviz.hpp>
#include <math.h>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/small_world_generator.hpp>
#include <boost/graph/property_maps/null_property_map.hpp>

using namespace std;
using namespace boost;

// define globally accessible graph
typedef adjacency_matrix<undirectedS> Graph;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef property_map<Graph, vertex_index_t>::type IndexMap;

// graph generator
typedef boost::small_world_iterator<boost::minstd_rand, Graph> SWGen;

Vertex v;

int processor_rank;
int randNum = 0;

typedef boost::unordered_map<int, int> boost_map;
boost_map vertex_to_position_map;

int total_no_of_processors;
int total_no_of_vertices;
int source_vertex;
int level = 1;
// std::vector<vector<int> > adjacencyMatrix(1000000);

class CustomGraph {
private:
      std::vector<vector<int> > adjacencyMatrix;
      int vertexCount;
public:
      //bool** adjacencyMatrix;
      //int vertexC = 1000000;
      
      //bool *adjacencyMatrix = (bool *)malloc(vertexCount * vertexCount * sizeof(bool));
      CustomGraph(int vertexCount) {
            this->vertexCount = vertexCount;
            adjacencyMatrix = vector<vector<int> >(total_no_of_vertices);
            //adjacencyMatrix = new bool*[vertexCount];
            /*for (int i = 0; i < vertexCount; i++) {
                  //adjacencyMatrix[i] = new bool[vertexCount];
                  for (int j = 0; j < vertexCount; j++)
                        *(adjacencyMatrix + i*vertexCount + j) = false;
                        //adjacencyMatrix[i][j] = false;
            }*/
      }
 
      void addEdge(int i, int j) {
            if (i >= 0 && i < vertexCount && j > 0 && j < vertexCount) {
                  //adjacencyMatrix[i][j] = true;
                  //adjacencyMatrix[j][i] = true;
                  //*(adjacencyMatrix + i*vertexCount + j) = true;
                  //*(adjacencyMatrix + j*vertexCount + i) = true;
                  adjacencyMatrix[i].push_back(j);
            }
      }
 
      /*void removeEdge(int i, int j) {
            if (i >= 0 && i < vertexCount && j > 0 && j < vertexCount) {
                  //adjacencyMatrix[i][j] = false;
                  //adjacencyMatrix[j][i] = false;
                  //*(adjacencyMatrix + i*vertexCount + j) = false;
                  //*(adjacencyMatrix + j*vertexCount + i) = false;
                  adjacencyMatrix[i].push_back(j);
            }
      }*/
 
      bool isEdge(int i, int j) {
            if(std::find(adjacencyMatrix[i].begin(), adjacencyMatrix[i].end(), j) != adjacencyMatrix[i].end())
                  return true;
            else
                  return false;
      }

      void printGraph(){
        for(int i = 0; i < vertexCount; i++){
           for(int j = 0; j < vertexCount; j++){
             if(isEdge(i, j)){
               std::cout << "Edge " << i << " " << j << std::endl;
             }
           }
        } 
      }
 
      /*~CustomGraph() {
            for (int i = 0; i < vertexCount; i++)
                  delete[] adjacencyMatrix[i];
            delete[] adjacencyMatrix;
      }*/
};

void create_graph(void);
void BFS(CustomGraph);
//void create_positions_map(Graph g);

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &processor_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &total_no_of_processors);

  if (argc < 3) {
    if (processor_rank == 0)
      std::cout << "Usage: <program_name> no-of-vertices source-vertex"
          << std::endl;
    return MPI_Finalize();
  }

  // argv[1] is no of vertices
  total_no_of_vertices = atoi(argv[1]);
  source_vertex = atoi(argv[2]);
  if (processor_rank == 0) {
    std::cout << " total no of vertices from command line : "
        << total_no_of_vertices << std::endl;
    std::cout << " source vertex from command line : " << source_vertex
        << std::endl;
  }

  if (processor_rank == 0) {
    std::cout << " total no of vertices : " << total_no_of_vertices
        << " total no of processors : " << total_no_of_processors
        << std::endl;
  }

  //MPI_Barrier (MPI_COMM_WORLD);

  create_graph();

  return MPI_Finalize();
}

/*void create_positions_map(Graph g) {
  IndexMap index = get(vertex_index, g);

  int position = 0;
  typedef graph_traits<Graph>::vertex_iterator vertex_iter;
  std::pair<vertex_iter, vertex_iter> vp;
  for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
    Vertex v = *vp.first;
    int key = index[v];
    vertex_to_position_map.insert(std::pair<int, int>(key, position));
    ++position;
  }

  if (processor_rank == 0)
    std::cout << "total size of map " << vertex_to_position_map.size()
        << std::endl;
}*/

void create_graph(void) {
  enum {
    A, B, C, D, E, F, G, H, I, J, K, L, M, N
  };
  //const char* name = "0123456789ABCD";

  //Graph g(total_no_of_vertices);
  /*for (int count = 1; count < total_no_of_vertices; ++count) {
    add_edge(0, count, g);
  }*/

  

  /*

   add_edge(F, A, g);
   add_edge(F, B, g);
   add_edge(F, C, g);
   add_edge(F, D, g);
   add_edge(F, E, g);

   add_edge(F, G, g);

   add_edge(E, L, g);

   add_edge(G, H, g);
   add_edge(G, I, g);
   add_edge(G, J, g);
   add_edge(G, K, g);
   add_edge(G, L, g);

   add_edge(L, M, g);
   */

  //add_edge(B, C, g);
  //add_edge(B, F, g);
  //add_edge(C, A, g);
  //add_edge(D, E, g);
  //add_edge(A, E, g);
  //add_edge(F, D, g);
  //std::cout << "vertex set: ";
  //boost::print_vertices(g, name);
  //std::cout << std::endl;
  //std::cout << "edge set: ";
  //boost::print_edges(g, name);
  //std::cout << std::endl;
  // std::cout << "out-edges: " << std::endl;
  // boost::print_graph(g, name);
  // std::cout << std::endl;
  // printing graph
  std::string name2;
  //boost::dynamic_properties dp(ignore_other_properties);

  //boost::dynamic_properties dp(boost::ignore_other_properties);

  //dp.property("node_id",     boost::get(&DotVertex::name, g));

  //dp.property("node_id",     boost::get(&DotVertex::name, graphviz));
  //auto index = get(vertex_color, g);
  //dp.property("node_id", index);
  //dp.property("node_id",     boost::get(&DotVertex::name, g));

  //if(processor_rank == 0){
  boost::minstd_rand gen;

  /*/ Create graph with 100 nodes
   Graph g2(SWGen(gen, total_no_of_vertices, 6, 1), SWGen(), total_no_of_vertices);
   std::ofstream outf("mil_graph.dot");
   write_graphviz(outf, g2);
  }*/

  // if (true)
    //std::cout << " star graph " << std::endl;

  //std::ifstream dot("test.dot");
  //boost::read_graphviz(dot, g, dp);
  /*
   std::ofstream outf("test2.dot");
   if (processor_rank == 0) {
   write_graphviz(outf, g);
   }

   */
   //std::cout << "before graph" << std::endl;
   CustomGraph g(total_no_of_vertices);
  //std::cout << "after graph 1" << std::endl;
  //if(processor_rank == 0){
  for(int i = 1; i < total_no_of_vertices; i++){
     g.addEdge(0,i);
   }
  //}

  MPI_Barrier (MPI_COMM_WORLD);

  if(processor_rank == 0)
    std::cout << "graph generation complete" << std::endl;


  //return;
  BFS(g);

}

int find_owner(int s_vertex) {

  //int position = vertex_to_position_map.at(s_vertex);

  if (total_no_of_processors == 1)
    return 0;

  return s_vertex % total_no_of_processors;
}

void BFS(CustomGraph g) {
  // get the property map for vertex indices
  //IndexMap index = get(vertex_index, g);

  // std::cout << "vertices(g) = " << endl;
  // typedef graph_traits<Graph>::vertex_iterator vertex_iter;
  // std::pair<vertex_iter, vertex_iter> vp;
  // for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
  //     //std::cout << index[v] <<  " ";
  //   Vertex v = *vp.first;
  //   distance_from_source.push_back(-1);

  // }

  vector<int> distance_from_source(total_no_of_vertices, -1);

  std::vector<int> FS;
  std::vector<int> NS;

  Vertex v;

  //find owner for source vertex
  int owner = find_owner(source_vertex);

  //MPI_Barrier (MPI_COMM_WORLD);
  int source_vertex_bcast_size = 1;

  if (owner == processor_rank) {
    //std::cout << " owner for vertex" << source_vertex << "  : " << owner
    //  << std::endl;
    //cout << " My processor_rank is " << processor_rank << endl;
    //cout << " pushed to queue: " << source_vertex << endl;
    FS.push_back(source_vertex);
    distance_from_source[source_vertex] = 0;
  }

  // Resize FS with bcast size and then bcast to all ranks
  // Resize for all ohers
  if (processor_rank != owner) {
    FS.resize(source_vertex_bcast_size);
    distance_from_source.resize(total_no_of_vertices);
  }

  // if(processor_rank == owner)
  MPI_Bcast(&FS[0], FS.size(), MPI_INT, owner, MPI_COMM_WORLD);
  MPI_Bcast(&distance_from_source[0], distance_from_source.size(), MPI_INT,
      owner, MPI_COMM_WORLD);

  //Resize distance_from_source with size as 1 and then bcast to all ranks
  //Resize for all others

  MPI_Barrier (MPI_COMM_WORLD);

  //if (processor_rank == 0)
  //std::cout << "FS[source_vertex] :" << FS[0] << std::endl;
  //MPI_Barrier(MPI_COMM_WORLD);

  int loopCount = 1;

  // printing newline
  if (processor_rank == 0)
    std::cout << endl;

  const clock_t begin_time = clock();
  double starttime, endtime, startTimeFindOwner, endTimeFindOwner,
      totalTimeFindOwner;

  if (processor_rank == 0)
    starttime = MPI_Wtime();

  MPI_Barrier(MPI_COMM_WORLD);

  while (!FS.empty()) {

    // Printing current FS
    if (processor_rank == 0) {
      std::cout << "------------------ FS ------------------"
          << std::endl;
      for (int i = 0; i < FS.size(); i++) {
        std::cout << "FS[" << i << "] = " << FS[i] << std::endl;
      }
    }

    // -------------- defining all local variables ----------------- //

    std::vector<vector<int> > recv_buffers(total_no_of_processors,
        vector<int>(total_no_of_vertices, -1));

    std::vector<vector<int> > local_buffers(total_no_of_processors);

    std::vector<int> recieve_count(total_no_of_processors);

    std::vector<int> send_count(total_no_of_processors);

    std::vector<int> recieve_ns_count(total_no_of_processors);

    std::vector<int> send_ns_count(total_no_of_processors);

    std::vector<int> updated_count(total_no_of_processors);

    std::vector<vector<int> > updated_vector(total_no_of_processors);

    std::vector<int> combinedNS(total_no_of_vertices, -1);

    int newFSSize = 0;

    //divide loop in parallel
    int for_start, for_end, queue_size;
    queue_size = FS.size();

    for_start = (queue_size / total_no_of_processors) * processor_rank;
    if (queue_size % total_no_of_processors > processor_rank) {
      for_start += processor_rank;
      for_end = for_start + (queue_size / total_no_of_processors) + 1;
    } else {
      for_start += queue_size % total_no_of_processors;
      for_end = for_start + (queue_size / total_no_of_processors);
    }
    // std::cout << "processor rank : " << processor_rank << " for_start " << for_start << " for_end " << for_end << endl;
    for (int count = for_start; count < for_end; ++count) {
      //std::cout << "inside parallel for loop " << endl;
      //take each u in FS
      int u = FS.at(count);

      //get adjacent vertices of u
      //typedef graph_traits<Graph>::adjacency_iterator adjacency_iterator;
      //std::pair<adjacency_iterator, adjacency_iterator> neighbours =
        //  adjacent_vertices(vertex(u, g), g);

      for (int i = 0; i < total_no_of_vertices; i++) {
        if(!g.isEdge(u, i)) continue;

        int v = i;
        if (distance_from_source[v] != -1) {
          continue;
        }

        //startTimeFindOwner = MPI_Wtime();

        owner = find_owner(v);

        //endTimeFindOwner = MPI_Wtime();

        //totalTimeFindOwner += endTimeFindOwner - startTimeFindOwner;

        if (owner != processor_rank) {
          local_buffers[owner].push_back(v);
          send_count[owner]++;
          // std::cout << "inside if vertex: " << v << " owner: "
          //     << owner << " rank: " << processor_rank
          //     << std::endl;
        } else {
          if (distance_from_source[v] == -1) {
            // means vertex has not yet been traversed then giving it value of level
            // and pushing it in NS list
            distance_from_source[v] = level;
            NS.push_back(v);
            updated_count[processor_rank] =
                updated_count[processor_rank] + 1;
            updated_vector[processor_rank].push_back(v);
            // std::cout << "inside else vertex: " << v << " owner: "
            //     << owner << " rank: " << processor_rank
            //     << std::endl;
          }
        }

        // ++local_buffers_size_metadata[owner];.\
        // MPI_Bcast(&local_buffers_size_metadata[0], local_buffers_size_metadata.size(), MPI_INT, count, MPI_COMM_WORLD);
        // MPI_Barrier(MPI_COMM_WORLD);
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < total_no_of_processors; i++) {
      // here gathering send counts into recv counts for each rank
      MPI_Gather(&send_count[i], 1, MPI_INT, &recieve_count[0], 1,
          MPI_INT, i, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    // std::cout << " line 306 " << std::endl;

    // --------------------------------------- from here --------------------------------------//

    //calculate displacement vector now
    std::vector<int> displacement_vector3(total_no_of_processors);
    int total_ns_size3 = 0;
    int temp_displ3 = 0;
    for (int count = 1; count < recieve_count.size(); ++count) {
        displacement_vector3[count] = displacement_vector3[count - 1]
            + recieve_count[count - 1];
    }

    // if(processor_rank == 0){
    //   std::cout << "------------------ Printing displ vector buffer for rank 0 ---------------------------" << std::endl;
    //   for(int count =0; count < recieve_ns_count.size(); ++count){
    //     std::cout << "displ "   << displacement_vector[count] << std::endl;
    //   }
    // }

    MPI_Barrier(MPI_COMM_WORLD);

    for (int count = 0; count < total_no_of_processors; ++count) {
    MPI_Gatherv(local_buffers[count].data(), local_buffers[count].size(), MPI_INT, recv_buffers[count].data(),
        recieve_count.data(), displacement_vector3.data(), MPI_INT, count,
        MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // ------------------------------- till here -------------------------------------//

    //MPI_Request request;

    /*/ for verifying printing local buffers of rank 1
     if(processor_rank == 1){
     std::cout << "--------------- Printing local buffers for process 1 --------------------" << std::endl;
     for(int j = 0; j < local_buffers.size(); j++){
     for(int i = 0; i < local_buffers[j].size(); i++){
     std::cout << "rank 1 local_buffers" << j << "[" << i << "] is " << local_buffers[j][i] << std::endl;
     }
     }
     }*/

    /*/ Now we will send data in local buffers to all respective processors
    for (int i = 0; i < total_no_of_processors; i++) {
      if (send_count[i]) {
        MPI_Isend(&local_buffers[i][0], send_count[i], MPI_INT, i, 1,
            MPI_COMM_WORLD, &request);
        MPI_Request_free(&request);
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
     */
    // std::cout << " line 330 " << std::endl;
    //  for verifying printing recv count buffer
    //    if(processor_rank == 0){
    //      std::cout << "------------------ Printing recieve count buffer for rank 0 ---------------------------" << std::endl;
    //      for(int i = 0; i < recieve_count.size(); i++){
    //        std::cout << "recieve_count[" << i << "] is " << recieve_count[i] << std::endl;
    //      }
    //    }
    // }

    std::vector<vector<int> > tempVector(total_no_of_processors,
        vector<int>(total_no_of_processors));

    /*/ Now recieving data
    for (int i = 0; i < total_no_of_processors; i++) {
      if (recieve_count[i]) {
        MPI_Recv(&recv_buffers[i][0], recieve_count[i], MPI_INT, i, 1,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }*/

    // for verification printing recv buffer of processor 3
    // if (processor_rank == 7) {
    //   std::cout
    //       << "--------------- Printing recv buffer for process 3 --------------------"
    //       << std::endl;
    //   for (int m = 0; m < recv_buffers.size(); m++) {
    //     for (int n = 0; n < recv_buffers[m].size(); n++) {
    //       std::cout << "rank 7 recv_buffer" << m << "[" << n
    //           << "] is " << recv_buffers[m][n] << std::endl;
    //     }
    //   }
    // }

    /*/ now our local buffers are synchronized in recv buffers we will loop in processor's recv
    for (int j = 0; j < recv_buffers.size(); j++) {
      if (recv_buffers[j].size() != 0) {
        for (int i = 0; i < recv_buffers[j].size(); i++) {
          if (recv_buffers[j][i] != -1) {
            if (distance_from_source[recv_buffers[j][i]] == -1) {
              // means vertex has not yet been traversed then giving it value of level
              // and pushing it in NS list
              distance_from_source[recv_buffers[j][i]] = level;
              NS.push_back(recv_buffers[j][i]);
              //std::cout << "processor rank inside updating loop " << processor_rank << " and level is " << level << std::endl;
              updated_count[processor_rank] =
                  updated_count[processor_rank] + 1;
              updated_vector[processor_rank].push_back(
                  recv_buffers[j][i]);
            }
          }
        }
      }
    }*/

    for (int j = 0; j < recv_buffers[processor_rank].size(); j++) {
      if (recv_buffers[processor_rank].size() != 0) {
          if (recv_buffers[processor_rank][j] != -1) {
            if (distance_from_source[recv_buffers[processor_rank][j]] == -1) {
              // means vertex has not yet been traversed then giving it value of level
              // and pushing it in NS list
              distance_from_source[recv_buffers[processor_rank][j]] = level;
              NS.push_back(recv_buffers[processor_rank][j]);
              //std::cout << "processor rank inside updating loop " << processor_rank << " and level is " << level << std::endl;
              updated_count[processor_rank] =
                  updated_count[processor_rank] + 1;
              updated_vector[processor_rank].push_back(
                  recv_buffers[processor_rank][j]);
            }
          }
      }
    }

    /*
     for(int i = 0; i < updated_vector[processor_rank].size(); i++){
     std:: cout << "updated vector[" << i << "] = " << updated_vector[processor_rank][i] << std::endl;
     }*/

    // here we will combile all NSs and put it in FS
    send_ns_count[processor_rank] = NS.size();

    //MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < total_no_of_processors; i++) {
      // here gathering send counts into recv counts for each rank
      MPI_Gather(&send_ns_count[processor_rank], 1, MPI_INT,
          &recieve_ns_count[0], 1, MPI_INT, i, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // / for verifying printing recieve ns count buffer
    // if (processor_rank == 0) {
    //   std::cout
    //       << "------------------ Printing ns counts recv buffer for rank 0 ---------------------------"
    //       << std::endl;
    //   for (int i = 0; i < recieve_ns_count.size(); i++) {
    //     std::cout << "recieve_ns_count[" << i << "] is "
    //         << recieve_ns_count[i] << std::endl;
    //   }
    // }

    int totalNS = 0;

    //calculate displacement vector now
    std::vector<int> displacement_vector(total_no_of_processors);
    int total_ns_size = 0;
    int temp_displ = 0;
    if (processor_rank == 0) {
      for (int count = 1; count < recieve_ns_count.size(); ++count) {
        displacement_vector[count] = displacement_vector[count - 1]
            + recieve_ns_count[count - 1];
      }
    }

    // if(processor_rank == 0){
    //   std::cout << "------------------ Printing displ vector buffer for rank 0 ---------------------------" << std::endl;
    //   for(int count =0; count < recieve_ns_count.size(); ++count){
    //     std::cout << "displ "   << displacement_vector[count] << std::endl;
    //   }
    // }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Gatherv(NS.data(), NS.size(), MPI_INT, combinedNS.data(),
        recieve_ns_count.data(), displacement_vector.data(), MPI_INT, 0,
        MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    // / for verification printing combined NS
    // if(processor_rank == 0){
    // std::cout << "Printing combined NS in rank 0" << std::endl;
    // for(int i = 0; i < combinedNS.size(); i++){
    // std::cout << "combinedNS[" << i << "] = " << combinedNS[i] << std::endl;
    // }
    // }

    // Push data from combine NS to FS by rank 0
    FS.clear();

    // refilling FS from combined NS vector by rank 0
    if (processor_rank == 0) {

      for (int i = 0; i < combinedNS.size(); i++) {
        if (combinedNS[i] != -1) {
          FS.push_back(combinedNS[i]);
          // std::cout << "new FS[" << i << "]" << combinedNS[i] << std::endl;
        }
      }

      // std::cout << " FS new size is " << FS.size() << std::endl;
    }

    // now broadcasting new FS size

    newFSSize = FS.size();

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(&newFSSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    if (processor_rank != 0) {
      FS.resize(newFSSize);
    }

    MPI_Bcast(&FS[0], FS.size(), MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    /*/ printing distance vector
     std::cout << "Distance vector of rank " << processor_rank << std::endl;
     for(int i = 0; i < distance_from_source.size(); i++){
     std::cout << "distance_from_source" << processor_rank << "[" << i << "] = " <<  distance_from_source[i] << std::endl;
     }*/

    // --------------------------- udpating distance vector for all ranks ---------------------------
    /*
     for (int i = 0; i < total_no_of_processors; i++) {
     // here gathering send counts into recv counts for each rank
     MPI_Gather(&updated_count[processor_rank], 1, MPI_INT,
     &updated_count[0], 1, MPI_INT, i, MPI_COMM_WORLD);
     }

     MPI_Barrier(MPI_COMM_WORLD);
     */
    /*/ for verifying printing recieve ns count buffer
     if(processor_rank == 0){
     std::cout << "------------------ Printing updated counts buffer for rank 0 ---------------------------" << std::endl;
     for(int i = 0; i < updated_count.size(); i++){
     std::cout << "updated_count[" << i << "] is " << updated_count[i] << std::endl;
     }
     }*/


    // Now after getting combined updated vectors we will udpate distance vector of rank 0
    if (processor_rank == 0) {
      for (int i = 0; i < combinedNS.size(); i++) {
        if (combinedNS[i] != -1)
          distance_from_source[combinedNS[i]] = level;
      }
    }

    // Now we will broadcast distance_from_source vector to all ranks
    MPI_Bcast(&distance_from_source[0], distance_from_source.size(),
        MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    // Now verifying result of broadcast by printing values of distance vector for rank 2
    // if (FS.empty()) {
    if (FS.empty()) {
      if (processor_rank == 0) {
        std::cout << endl;
        std::cout
            << "------------------- resulting distance_from_source ---------------------"
            << std::endl;
        for (int i = 0; i < distance_from_source.size(); i++) {
          std::cout << "distance_from_source[" << i << "] = "
              << distance_from_source[i] << std::endl;
        }
        std::cout << endl;
        // printing elapsed time
        std::cout << "----------- total running time --------------" << std::endl;
        std::cout << "total no of processors :" << total_no_of_processors << " total no of vertices :" << total_no_of_vertices  << std::endl;
        std::cout << float(clock() - begin_time) / CLOCKS_PER_SEC;
        std::cout << endl;
      }
    }

    // before loop continution incrementing level value
    level++;

    // clearing NS for all ranks
    NS.clear();

    MPI_Barrier(MPI_COMM_WORLD);

    /*/ printing time taken
     if(FS.empty()){
     std::cout << std::endl;
     //std::cout << "totalTimeFindOwner for rank " << processor_rank << " = " << totalTimeFindOwner << std::endl;
     std::cout << endl;

     if (processor_rank == 0) {
     // printing elapsed time
     endtime   = MPI_Wtime();
     std::cout << (endtime-starttime);
     std::cout << endl;
     }
     }*/
  }

  //  std::cout << "adjacent vertices: " << endl;
  //   typename graph_traits<Graph>::adjacency_iterator ai;
  //   typename graph_traits<Graph>::adjacency_iterator ai_end;
  //   Vertex v = *vp.first;

  //   for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
  //    Vertex v = *vp.first;
  //    std::cout << index[v] <<  " -->" ;
  //     for (boost::tie(ai, ai_end) = adjacent_vertices(v, g);
  //         ai != ai_end; ++ai)
  //      std::cout << index[*ai] <<  " ";
  //  std::cout << std::endl;

  // }
}


