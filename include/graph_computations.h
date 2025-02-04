#ifndef _graph_computations_included_
#define _graph_computations_included_

#include "localization.h"
#include "solution.h"

//define multi layer graph class, capable of representing the needed bipartite graph with edge sets e1,e2,e3
class CandidateGraph {
	Graph g;
	std::vector<int> left_vs, right_vs;
	int num_polys1, num_polys2;
	//per vertex v, store the vertices, where the represented polygons area majorly included (> 0.5 of their area) in the polygons represented by v
	std::vector<std::vector<int>> inclusion_lookup;
	//per vertex v, store if it is known the vertex that it has to be in a match with 
	std::vector<int> fix_matches;

public:
	//initiates an empty multilayer graph
	CandidateGraph();
	//initiates the multilayer graph with a vertex for every input polygon
	CandidateGraph(std::vector<Polygon_wh>* polys1, std::vector<Polygon_wh>* polys2);
	//returns a deepcopy of the entire CandidateGraph
	CandidateGraph copy();
	//adds edge to multilayer graph in between the vertices 'source' and 'target' in layer 'layer' (in [0,2]) 
	void add_edge(int source, int target, double weight);
	//overloads add_edge with no given weight
	void add_edge(int source, int target);
	//deletes the edge from the layer
	void delete_edge(int source, int target);

	//adds a vertex to the multilayer graph and returns its index
	int add_vertex(bool ref_map, std::vector<int> ref_polys);
	//deletes the vertex 'v_id' and all its incident edges from the multilayergraph
	void delete_vertex(int v_id);
	//edits vertex information of vertex 'v_id' and applies to all layers
	void edit_vertex(int v_id, bool ref_map, std::vector<int> ref_polys);
	//returns a vertex of the graph
	VertexProps get_vertex(int v_id);

	std::vector<int> get_Vi(bool map);


	//get amount of vertices in the graph
	int num_vertices();

	//returns the amount of edges in layer l
	int num_edges();


	//returns the requested graph layer
	Graph* get_graph();

	//returns the total number of connected components and writes the component id per vertex in the given vector
	int get_components(std::vector<int>* components);

	//returns true, if the vertices 'source' and 'target' are connected in layer 'layer' 
	bool are_connected(int source, int target);

	//returns -1, if vertex representing the given components does not exist, else returns the index of the vertex
	int does_vertex_exist(bool map, const std::vector<int>& ref_polys);

	bool referenced_map(int v);

	void addIncludedVertex(int v_id, int included_vertex);
	void setIncludedVertices(int v_id, std::vector<int> included_vertices);
	std::vector<int> getIncludedVertices(int v_id);

	void setFixMatch(int v_1, int v_2);
	int getFixMatch(int v_id);

	//prints overview of all vertices and represented polygons
	void printVertexOverview();
};

//given a vector of graphs containing nodes only, computes the weighted edges
void computeEdges(std::vector<CandidateGraph>& g_vec, Localization osm_rtree, Localization atkis_rtree, const MapOverlay& mo, double lambda, const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2);

//a simple depth first search on the graph g, initialized on vertex v and excluding all vertices with index >= upper_limit
void DFS(Graph* g, int v, int upper_limit, std::vector<bool>* visited, std::vector<int>* component);

void completeCandidateGraph(CandidateGraph* cg, const MapOverlay& mo, double lambda, std::chrono::steady_clock::time_point start_time, int time_limit, int size_limit);

//uses a guessing technique to find cases fulfilling the 1:1-Matches-Lemma
//the matches are directly added to the solution and the corresponding vertices are removed from the CandidateGraph
void precomputeSimpleMatches(CandidateGraph* cg, const MapOverlay& mo, double lambda, Solution* sol);

// takes inclusion into consideration to form groups of polygons, which will always be matched together
// replaces single vertices representing those polygons by cumulative vertices
// does only need a graph with vertices and correct inclusion lookup (see setInclusionProperties), edges are not considered
// returns true if the graph was modified, else false (on modification, every edge will be deleted and E has to be recomputed)
bool pregroupIncludedPolygons(CandidateGraph* cg, const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2, const Localization& osm_rtree, const Localization& atkis_rtree, const MapOverlay& mo, double lambda);

//scans a graph for inclusion properties between polygons, as this info is needed for pregrouping
//adds the inclusion relations to the CandidateGraphs in g_vec
void setInclusionProperties(std::vector<CandidateGraph> &g_vec, const Localization& osm_rtree, const Localization& atkis_rtree, const MapOverlay& mo, double lambda, const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2);

#endif