#ifndef POLYGONMATCHING_TREE_COMPUTATIONS_H
#define POLYGONMATCHING_TREE_COMPUTATIONS_H

#include "cgal_includes.h"
#include "polygon_wh.h"
#include "localization.h"
#include "graph_computations.h"
#include "logger.h"

#include <boost/graph/copy.hpp>

#include "command_line_parser.h"

// Define the edge property to include a weight
struct BiGraphEdgeProperties {
    double weight;
};

// Define the vertex property to include inclusion_score and represented_polys
struct BiGraphVertexProperties {
    //remember an id for lookups as vertices are deleted on merging operations
    int id;
    double inclusion_score;
    bool represented_map;
    std::vector<int> represented_polys;
};

// Define bidirectional graph using an adjacency list with directed edges
// note that we use directed edges instead of bidirectional ones in order to be able to specify individual weights
typedef boost::adjacency_list<boost::listS, boost::vecS, boost::directedS,
        BiGraphVertexProperties, BiGraphEdgeProperties> BiGraph;

// Define the types for the vertices and edges
typedef boost::graph_traits<Graph>::vertex_descriptor BiVertex;
typedef boost::graph_traits<Graph>::edge_descriptor BiEdge;

// Define an enum class for the mode.
enum class PathtoRootMode {
    START_AT_VERTEX,
    START_AT_POLYGON_ID
};

//struct to store objective information
struct ObjectiveInfo {
    OBJECTIVE objective = OBJECTIVE::JACCARD;
    std::pair<double, double> weights = std::make_pair(0.0, 0.0);
    double lambda = 0.0;
};


//class to store the bipartite candidate graph including the tree constraints of both sides
class TreeConstrainedCandidateGraph {
    //stores the left and right tree modeling the constraints of pairs to be matched
    Graph g;

    //keep track of the current maximum edge id for new inserts
    int max_edge_id;

    //vector storing the corresponding tree vertex index for each vertex in g
    std::vector<int> corresponding_tree_vertex;
    std::vector<std::vector<int>> corresponding_graph_vertex;

    //store leafs of each tree and the represented polygon id
    //the vector is sorted w.r.t. the polygon ids represented by the vertex
    //note that a leaf can represent a set of polygons due to simplifying precomputations!
    std::vector<std::vector<Vertex>> leafs;

    std::vector<DiGraph> trees;
    //stores the root indices of both trees
    std::vector<int> roots;
public:
    //constructor based on poly sets
    TreeConstrainedCandidateGraph(const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2){

        //initialize trees with the vertices of the Candidate Graph, such that they have the same indices
        //(with an offset of polys1.size() for the vertices in the atkis tree)
        this->g = Graph();
        this->corresponding_tree_vertex = std::vector<int>();
        this->corresponding_graph_vertex = std::vector<std::vector<int>>(2);
        this->trees = std::vector<DiGraph>(2);
        this->roots = {-1,-1};
        this->leafs = std::vector<std::vector<Vertex>>(2);
        this->max_edge_id = 0;
        for(int o_p  = 0; o_p < polys1.size(); o_p++) {
            Vertex v_new = boost::add_vertex(this->g);
            this->g[v_new].referenced_map = false;
            this->g[v_new].referenced_polys = {o_p};
            //do not set attributes in the trees to avoid redundant data
            this->corresponding_graph_vertex[0].push_back(v_new);
            v_new = boost::add_vertex(this->trees[0]);
            this->corresponding_tree_vertex.push_back(v_new);
            this->leafs[0].push_back(v_new);
        }
        for(int a_p = 0; a_p < polys2.size(); a_p++) {
            Vertex v_new = boost::add_vertex(this->g);
            this->g[v_new].referenced_map = false;
            this->g[v_new].referenced_polys = {a_p};
            //do not set attributes in the trees to avoid redundant data
            this->corresponding_graph_vertex[1].push_back(v_new);
            v_new = boost::add_vertex(this->trees[1]);
            this->corresponding_tree_vertex.push_back(v_new);
            this->leafs[1].push_back(v_new);
        }

    }

    //constructor based on a CandidateGraph
    //note that the candidate graph is expected to for each polygon have at most one vertex, where the polygon
    //is amongst the represented polygons! (Otherwise building tree structures would already be implied or impossible)
    TreeConstrainedCandidateGraph(CandidateGraph& cg) {
        this->g = Graph();
        this->corresponding_tree_vertex = std::vector<int>();
        this->corresponding_graph_vertex = std::vector<std::vector<int>>(2);
        this->trees = std::vector<DiGraph>(2);
        this->roots = {-1,-1};
        this->leafs = std::vector<std::vector<Vertex>>(2);
        this->max_edge_id = 0;

        //candidate graph gives no information about hierarchy, simply copy all vertices
        for(int v_id = 0; v_id < cg.num_vertices(); v_id++) {
            Vertex v_new = boost::add_vertex(this->g);
            this->g[v_new].referenced_polys = cg.get_vertex(v_id).referenced_polys;
            this->g[v_new].referenced_map = cg.get_vertex(v_id).referenced_map;

            bool map = this->g[v_new].referenced_map;

            this->corresponding_graph_vertex[this->g[v_new].referenced_map].push_back(v_new);
            Vertex v_new_tree = boost::add_vertex(this->trees[this->g[v_new].referenced_map]);
            this->corresponding_tree_vertex.push_back(v_new_tree);

            //treat every vertex as a leaf
            this->leafs[map].push_back(v_new_tree);

        }

        //set roots if trivially possible
        for (int m=0; m<2; m++) {
            if (boost::num_vertices(this->trees[m]) == 1) this->roots[m] = 0;
        }

    }

    //specify copy operation
    TreeConstrainedCandidateGraph(TreeConstrainedCandidateGraph& cg_tree) {
        //copy simple data
        this->max_edge_id = cg_tree.max_edge_id;
        this->corresponding_tree_vertex = std::vector<int>(cg_tree.corresponding_tree_vertex.size());
        for(int i=0; i< cg_tree.corresponding_tree_vertex.size();i++) this->corresponding_tree_vertex[i]=cg_tree.corresponding_tree_vertex[i];
        this->corresponding_graph_vertex = std::vector<std::vector<int>>(2);
        for(int g=0; g<2; g++) {
            for(int i=0; i< cg_tree.corresponding_graph_vertex[g].size();i++) {
                this->corresponding_graph_vertex[g].push_back(cg_tree.corresponding_graph_vertex[g][i]);
            }
        }
        this->leafs = std::vector<std::vector<Vertex>>(2);
        for(int g=0; g<2; g++) {
            for(int i=0; i < cg_tree.leafs[g].size(); i++) {
                this->leafs[g].push_back(cg_tree.leafs[g][i]);
            }
        }
        this->roots = {cg_tree.roots[0],cg_tree.roots[1]};

        //copy graphs
        this->g = Graph();
        boost::copy_graph(cg_tree.g,this->g);

        this->trees = std::vector<DiGraph>(2);
        for(int i=0; i< 2; i++) boost::copy_graph(cg_tree.trees[i],this->trees[i]);

    }


    //DEBUG: constructor for debug instance
    TreeConstrainedCandidateGraph(bool dbg = false) {
        this->g = Graph();

        this->trees = std::vector<DiGraph>(2);
        this->corresponding_graph_vertex = std::vector<std::vector<int>>(2);
        this->corresponding_tree_vertex = std::vector<int>();

        this->roots = std::vector<int>(2);

        //create two stars with three leafs each
        for(int t=0; t<2;t++) {
            for (int i = 0; i < 4; i++) {
                auto v_new_g = boost::add_vertex(this->g);
                this->corresponding_tree_vertex.push_back(v_new_g);

                auto v_new_t = boost::add_vertex(this->trees[t]);
                this->corresponding_graph_vertex[t].push_back(v_new_t);
            }
            //edges to root
            for(int i=1; i<4;i++) {
                boost::add_edge(i,0,this->trees[t]);
            }

            this->roots[t] = 0;
        }

        //edges to other root
        for(int t=0; t<2;t++) {
            for(int i=1; i<4;i++) {
                auto e= boost::add_edge(this->corresponding_graph_vertex[t][i],this->corresponding_graph_vertex[(t+1)%2][0],this->g);
                this->g[e.first].weight = 1.0;
            }
        }



    }
    //END DEBUG


    //adds edge to bipartite graph in between the vertices 'source' and 'target'
    void add_edge(int source, int target, std::optional<double> weight);

    //adds edge to tree graph [map]
    //expects the vertex indices of the general graph, translates the indices into tree indices
    void add_tree_edge(bool map, int source, int target);

    //deletes the edge
    void delete_edge(int source, int target);

    //computes the weightes edges in the bipartite graph
    //note that it expects the tree constraints and roots to be defined as those are used
    //the polygons are optional arguments, as those are only needed for combined distances, not for the jaccard index
    void computeBipartiteWeightedEdges(const ObjectiveInfo& obj_info, const MapOverlay& mo, const std::vector<Polygon_wh>& polys1 = {}, const std::vector<Polygon_wh>& polys2 = {});

    //adds a vertex to the multilayer graph and returns its index
    int add_vertex(bool ref_map, std::vector<int> ref_polys);
    //deletes the vertex 'v_id' and all its incident edges from the multilayergraph
    void delete_vertex(int v_id);
    //returns a vertex of the graph
    VertexProps get_vertex(int v_id);

    Vertex get_root(bool map);
    //set an individual root based on the vertex id in the bipartite graph
    void set_root(bool map, int root);
    void set_roots(std::vector<int> roots);
    //check if a biparite vertex is a root in its corresponding tree
    bool is_root(int v_id);

    Graph& get_graph();

    //prints overview of all vertices and represented polygons
    void printVertexOverview();
    void printTrees();
    //return overviews as strings for logging
    std::string StringVertexOverview();
    std::string StringTrees();

    //checks for both trees if they have size = 1
    //if this is the case, the root is set to that vertex
    void setRootsforSingleVertexTrees();

    //recomputes the leafs of both trees and stores them in the local vectors
    void updateLeafs();

    //returns the tree ID that v is in
    int getTree(Vertex v);

    //returns the parent of the vertex in its tree
    Vertex getParent(Vertex v);

    //sets direction of directed edges in the trees
    //0 = downwards (root to leafs)
    //1 = upwards (leafs to root)
    //if the direction is not fulfilled, edges are inverted
    void setTreeDirection(bool dir);

    //returns a vector containing all vertices on the path from
    //the leaf representing source to the root in the corresponding tree
    std::vector<Vertex> getPathToRoot(bool tree, int source, PathtoRootMode mode);

    //returns a vector containing all vertices in the subtree rooted at source
    std::vector<Vertex> getSubtree(bool tree, int source);

    //returns maximum edge id + 1
    int getMaxEdgeID();

    //returns the vertex indices of the tree
    std::vector<Vertex> getLeafs(int tree);

    //returns the vertex indices of all pairs of leafs, meaning those pairs of vertices u,v that are both leafs and
    //have the same parent
    std::vector<std::pair<Vertex,Vertex>> getLeafPairs();
};


//takes an intersection graph (bipartite) and creates a bidirectional graph (bipartite), where the edge
// weights of a directed edge (u,v) correspond to intersection_area(u,v) / area(u)
// does use the edges of cg as information if those exist, else computes them from scratch using the overlay
BiGraph buildBiGraphFromIntersectionGraph(CandidateGraph& cg, const MapOverlay& mo);

//takes a bidirectional intersection graph with edge weights corresponding to the inclusion ratios
//note that the information about referenced polygons are provided by the candidate graph, which is why the vertex
//indices must be equal to those in g_bi, which is assured by buildBiGraphFromIntersectionGraph
//computes a tree constrained candidate graph from it via a merging process of vertices
TreeConstrainedCandidateGraph buildTreesFromBiGraph(BiGraph g_bi, CandidateGraph& cg, const MapOverlay& mo, double lambda);

//builds the tree constrained candidate graph via a simple kruskal strategy
TreeConstrainedCandidateGraph buildTreesViaKruskal(const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2, CandidateGraph& cg, MapOverlay& mo, double lambda);



#endif //POLYGONMATCHING_TREE_COMPUTATIONS_H
