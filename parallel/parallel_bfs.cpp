#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <mpi.h>


// boost headers
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/adjacency_iterator.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/small_world_generator.hpp>
#include <boost/graph/rmat_graph_generator.hpp>
#include <boost/graph/property_maps/null_property_map.hpp>

using namespace std;
using namespace boost;

// define globally accessible graph
typedef adjacency_matrix<undirectedS> Graph;

// typedef adjacency_list<> Graph;
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
int no_of_vertices_per_graph;
int multiplier;

int level = 1;
vector<vector<long int> > adjacencyMatrix;
void create_graph(void);
void BFS();
//void create_positions_map(Graph g);

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &processor_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &total_no_of_processors);

	if (argc < 3) {
		if (processor_rank == 0)
			cout
					<< "Usage: <program_name> no-of-vertices no-of-vertices-per-graph source-vertex "
					<< endl;
		return MPI_Finalize();
	}

	// argv[1] is no of vertices
	total_no_of_vertices = atoi(argv[1]);
	source_vertex = atoi(argv[3]);
	no_of_vertices_per_graph = atoi(argv[2]);
	multiplier = total_no_of_vertices / no_of_vertices_per_graph;

	if (processor_rank == 0) {
		cout << " total no of vertices from command line : "
				<< total_no_of_vertices << endl;
		cout << " source vertex from command line : " << source_vertex
				<< endl;

		cout << " no of vertices per graph from command line : "
				<< no_of_vertices_per_graph << endl;
		cout << " no of graphs : " << multiplier << endl;
	}

	if (processor_rank == 0) {
		cout << " total no of vertices : " << total_no_of_vertices
				<< " total no of processors : " << total_no_of_processors
				<< endl;
	}

	create_graph();

	return MPI_Finalize();
}

void create_graph(void) {
	enum {
		A, B, C, D, E, F, G, H, I, J, K, L, M, N
	};

	boost::minstd_rand gen;

	adjacencyMatrix.resize(total_no_of_vertices);

	for (int i = 0; i < multiplier; ++i) {

		// generating new random graph
		Graph g(SWGen(gen, no_of_vertices_per_graph, 10, 1), SWGen(),
				no_of_vertices_per_graph);
		// adding edge between this graph and previous graph
		if (i > 0)
			adjacencyMatrix[((i * no_of_vertices_per_graph)) - 1].push_back(
					((i * no_of_vertices_per_graph)));

		for (long int count = 0; count < no_of_vertices_per_graph; ++count) {

			IndexMap index = get(vertex_index, g);
			typedef graph_traits<Graph>::adjacency_iterator adjacency_iterator;
			pair<adjacency_iterator, adjacency_iterator> neighbours =
					adjacent_vertices(vertex(count, g), g);

			for (; neighbours.first != neighbours.second; ++neighbours.first) {
				long int v = index[*neighbours.first];
				adjacencyMatrix[((i * no_of_vertices_per_graph) + count)].push_back(
						((i * no_of_vertices_per_graph) + v));
			}
		}
	} //for ends here

	//write graph to file
	if (processor_rank == 0) {

		ofstream outff;
		string file_name;
		file_name = "graph_" + to_string(total_no_of_vertices) + ".graph";

		outff.open(file_name, ios_base::app);

		// writing adjacency matrix to file
		for (int i = 0; i < adjacencyMatrix.size(); ++i) {
			for (int j = 0; j < adjacencyMatrix[i].size(); j++) {
				outff << i << " " << adjacencyMatrix[i][j] << endl;
			}
		}

	}

	cout << "Total size of resultant graph is " << adjacencyMatrix.size()
			<< endl;

	MPI_Barrier (MPI_COMM_WORLD);

	if (processor_rank == 0)
		cout << "graph generation complete" << endl;

	// return;
	BFS();

}

int find_owner(int s_vertex) {

	//int position = vertex_to_position_map.at(s_vertex);

	if (total_no_of_processors == 1)
		return 0;

	return s_vertex % total_no_of_processors;
}

void BFS() {

	vector<int> distance_from_source(total_no_of_vertices, -1);

	vector<int> FS;
	vector<int> NS;

	Vertex v;

	//find owner for source vertex
	int owner = find_owner(source_vertex);

	int source_vertex_bcast_size = 1;

	if (owner == processor_rank) {
		//cout << " owner for vertex" << source_vertex << "  : " << owner
		//  << endl;
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

	MPI_Bcast(&FS[0], FS.size(), MPI_INT, owner, MPI_COMM_WORLD);
	MPI_Bcast(&distance_from_source[0], distance_from_source.size(), MPI_INT,
			owner, MPI_COMM_WORLD);

	int loopCount = 1;

	if (processor_rank == 0)
		cout << endl;

	double starttime, endtime, startTimeFindOwner, endTimeFindOwner,
			totalTimeFindOwner;

	MPI_Barrier (MPI_COMM_WORLD);

	if (processor_rank == 0)
		starttime = MPI_Wtime();

	while (!FS.empty()) {

		/*/ Printing current FS
		 if (processor_rank == 0) {
		 cout << "------------------ FS ------------------"
		 << endl;
		 for (int i = 0; i < FS.size(); i++) {
		 cout << "FS[" << i << "] = " << FS[i] << endl;
		 }
		 }*/

		// -------------- defining all local variables ----------------- //
		vector<vector<int> > recv_buffers(total_no_of_processors,
				vector<int>(total_no_of_vertices, -1));

		vector<vector<int> > local_buffers(total_no_of_processors);

		vector<int> recieve_count(total_no_of_processors);

		vector<int> send_count(total_no_of_processors);

		vector<int> recieve_ns_count(total_no_of_processors);

		vector<int> send_ns_count(total_no_of_processors);

		vector<int> updated_count(total_no_of_processors);

		vector<vector<int> > updated_vector(total_no_of_processors);

		vector<int> combinedNS(total_no_of_vertices, -1);

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

		for (int count = for_start; count < for_end; ++count) {
			//take each u in FS
			int u = FS.at(count);

			for (int i = 0; i < adjacencyMatrix[u].size(); i++) {

				int v = adjacencyMatrix[u][i];
				// cout << " u-v" << i <<" " << v << endl;
				if (distance_from_source[v] != -1) {
					continue;
				}

				owner = find_owner(v);

				if (owner != processor_rank) {
					local_buffers[owner].push_back(v);
					send_count[owner]++;
					// cout << "inside if vertex: " << v << " owner: "
					//     << owner << " rank: " << processor_rank
					//     << endl;
				} else {
					if (distance_from_source[v] == -1) {
						// means vertex has not yet been traversed then giving it value of level
						// and pushing it in NS list
						distance_from_source[v] = level;
						NS.push_back(v);
						updated_count[processor_rank] =
								updated_count[processor_rank] + 1;
						updated_vector[processor_rank].push_back(v);
						// cout << "inside else vertex: " << v << " owner: "
						//     << owner << " rank: " << processor_rank
						//     << endl;
					}
				}

			}
		}

		for (int i = 0; i < total_no_of_processors; i++) {
			// here gathering send counts into recv counts for each rank
			MPI_Gather(&send_count[i], 1, MPI_INT, &recieve_count[0], 1,
					MPI_INT, i, MPI_COMM_WORLD);
		}

		//calculate displacement vector now
		vector<int> displacement_vector3(total_no_of_processors);
		int total_ns_size3 = 0;
		int temp_displ3 = 0;
		for (int count = 1; count < recieve_count.size(); ++count) {
			displacement_vector3[count] = displacement_vector3[count - 1]
					+ recieve_count[count - 1];
		}

		// if(processor_rank == 0){
		//   cout << "------------------ Printing displ vector buffer for rank 0 ---------------------------" << endl;
		//   for(int count =0; count < recieve_ns_count.size(); ++count){
		//     cout << "displ "   << displacement_vector[count] << endl;
		//   }
		// }

		for (int count = 0; count < total_no_of_processors; ++count) {
			MPI_Gatherv(local_buffers[count].data(),
					local_buffers[count].size(), MPI_INT,
					recv_buffers[count].data(), recieve_count.data(),
					displacement_vector3.data(), MPI_INT, count,
					MPI_COMM_WORLD);
		}

		//MPI_Request request;

		/*/ for verifying printing local buffers of rank 1
		 if(processor_rank == 1){
		 cout << "--------------- Printing local buffers for process 1 --------------------" << endl;
		 for(int j = 0; j < local_buffers.size(); j++){
		 for(int i = 0; i < local_buffers[j].size(); i++){
		 cout << "rank 1 local_buffers" << j << "[" << i << "] is " << local_buffers[j][i] << endl;
		 }
		 }
		 }*/

		// cout << " line 330 " << endl;
		//  for verifying printing recv count buffer
		//    if(processor_rank == 0){
		//      cout << "------------------ Printing recieve count buffer for rank 0 ---------------------------" << endl;
		//      for(int i = 0; i < recieve_count.size(); i++){
		//        cout << "recieve_count[" << i << "] is " << recieve_count[i] << endl;
		//      }
		//    }
		// }
		vector<vector<int> > tempVector(total_no_of_processors,
				vector<int>(total_no_of_processors));

		// for verification printing recv buffer of processor 7
		// if (processor_rank == 7) {
		//   cout
		//       << "--------------- Printing recv buffer for processor 7 --------------------"
		//       << endl;
		//   for (int m = 0; m < recv_buffers.size(); m++) {
		//     for (int n = 0; n < recv_buffers[m].size(); n++) {
		//       cout << "rank 7 recv_buffer" << m << "[" << n
		//           << "] is " << recv_buffers[m][n] << endl;
		//     }
		//   }
		// }

		for (int j = 0; j < recv_buffers[processor_rank].size(); j++) {
			if (recv_buffers[processor_rank].size() != 0) {
				if (recv_buffers[processor_rank][j] != -1) {
					if (distance_from_source[recv_buffers[processor_rank][j]]
							== -1) {
						// means vertex has not yet been traversed then giving it value of level
						// and pushing it in NS list
						distance_from_source[recv_buffers[processor_rank][j]] =
								level;
						NS.push_back(recv_buffers[processor_rank][j]);
						//cout << "processor rank inside updating loop " << processor_rank << " and level is " << level << endl;
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
		  cout << "updated vector[" << i << "] = " << updated_vector[processor_rank][i] << endl;
		 }*/

		// here we will combile all NSs and put it in FS
		send_ns_count[processor_rank] = NS.size();

		for (int i = 0; i < total_no_of_processors; i++) {
			// here gathering send counts into recv counts for each rank
			MPI_Gather(&send_ns_count[processor_rank], 1, MPI_INT,
					&recieve_ns_count[0], 1, MPI_INT, i, MPI_COMM_WORLD);
		}

		// / for verifying printing recieve ns count buffer
		// if (processor_rank == 0) {
		//   cout
		//       << "------------------ Printing ns counts recv buffer for rank 0 ---------------------------"
		//       << endl;
		//   for (int i = 0; i < recieve_ns_count.size(); i++) {
		//     cout << "recieve_ns_count[" << i << "] is "
		//         << recieve_ns_count[i] << endl;
		//   }
		// }

		int totalNS = 0;

		//calculate displacement vector now
		vector<int> displacement_vector(total_no_of_processors);
		int total_ns_size = 0;
		int temp_displ = 0;
		if (processor_rank == 0) {
			for (int count = 1; count < recieve_ns_count.size(); ++count) {
				displacement_vector[count] = displacement_vector[count - 1]
						+ recieve_ns_count[count - 1];
			}
		}

		// if(processor_rank == 0){
		//   cout << "------------------ Printing displ vector buffer for rank 0 ---------------------------" << endl;
		//   for(int count =0; count < recieve_ns_count.size(); ++count){
		//     cout << "displ "   << displacement_vector[count] << endl;
		//   }
		// }

		MPI_Gatherv(NS.data(), NS.size(), MPI_INT, combinedNS.data(),
				recieve_ns_count.data(), displacement_vector.data(), MPI_INT, 0,
				MPI_COMM_WORLD);

		// / for verification printing combined NS
		// if(processor_rank == 0){
		// cout << "Printing combined NS in rank 0" << endl;
		// for(int i = 0; i < combinedNS.size(); i++){
		// cout << "combinedNS[" << i << "] = " << combinedNS[i] << endl;
		// }
		// }

		// Push data from combine NS to FS by rank 0
		FS.clear();

		// refilling FS from combined NS vector by rank 0
		if (processor_rank == 0) {

			for (int i = 0; i < combinedNS.size(); i++) {
				if (combinedNS[i] != -1) {
					FS.push_back(combinedNS[i]);
					// cout << "new FS[" << i << "]" << combinedNS[i] << endl;
				}
			}

			// cout << " FS new size is " << FS.size() << endl;
		}

		// now broadcasting new FS size

		newFSSize = FS.size();

		MPI_Bcast(&newFSSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

		if (processor_rank != 0) {
			FS.resize(newFSSize);
		}

		MPI_Bcast(&FS[0], FS.size(), MPI_INT, 0, MPI_COMM_WORLD);

		/*/ printing distance vector
		 cout << "Distance vector of rank " << processor_rank << endl;
		 for(int i = 0; i < distance_from_source.size(); i++){
		 cout << "distance_from_source" << processor_rank << "[" << i << "] = " <<  distance_from_source[i] << endl;
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

		if (FS.empty()) {
			if (processor_rank == 0) {
				endtime = MPI_Wtime();
				cout << endl;
				cout
						<< "------------------- resulting distance_from_source ---------------------"
						<< endl;
				for (int i = 0; i < distance_from_source.size(); i++) {
					cout << "distance_from_source[" << i << "] = "
							<< distance_from_source[i] << endl;
				}
				cout << endl;
				// printing elapsed time
				cout << "----------- total running time --------------"
						<< endl;
				cout << "total no of processors :"
						<< total_no_of_processors << " total no of vertices :"
						<< total_no_of_vertices << endl;
				// printing elapsed time
				cout << "total time";
				cout << endl;
				cout << fixed;
				cout << (endtime - starttime) << endl;
			}
		}

		// before loop continution incrementing level value
		level++;

		// clearing NS for all ranks
		NS.clear();

		MPI_Barrier(MPI_COMM_WORLD);

	}

}

