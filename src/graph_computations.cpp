#include <utility>

#include "../include/graph_computations.h"

//MULTILAYERGRAPH CLASS
CandidateGraph::CandidateGraph() {
	this->g = Graph();
	this->left_vs = std::vector<int>();
	this->right_vs = std::vector<int>();
	this->inclusion_lookup = std::vector<std::vector<int>>();
	this->fix_matches = std::vector<int>();
}

CandidateGraph::CandidateGraph(std::vector<Polygon_wh>* polys1, std::vector<Polygon_wh>* polys2) {
	//create three empty graphs
	this->g = Graph();
	this->left_vs = std::vector<int>();
	this->right_vs = std::vector<int>();
	this->inclusion_lookup = std::vector<std::vector<int>>();
	this->fix_matches = std::vector<int>();

	//remember polygon count
	this->num_polys1 = polys1->size();
	this->num_polys2 = polys2->size();

	//insert one node per poly1 and poly2
	//the indexing is [0..,polys1_size-1] for the polys1 and [polys1_size, ... , polys2_size-1 + polys1_size] for the polys2
    //add polys1
    for (int j = 0; j < polys1->size(); j++) {
        int new_vertex = boost::add_vertex(this->g);
        this->g[new_vertex].referenced_map = false;
        this->g[new_vertex].referenced_polys = std::vector<int>{ j };
        //remember vertex location for quick retreival
        this->left_vs.push_back(new_vertex);
    }

    //add polys2
    for (int j = polys1->size(); j < polys1->size() + polys2->size(); j++) {
        int new_vertex = boost::add_vertex(this->g);
        this->g[new_vertex].referenced_map = true;
        this->g[new_vertex].referenced_polys = std::vector<int>{ j - (int)polys1->size()};
        //remember vertex location for quick retreival
        this->right_vs.push_back(new_vertex);
    }



}

CandidateGraph CandidateGraph::copy() {
	CandidateGraph g_copy;

	//iterate over all vertices and insert them into the new graph
	Graph g = *this->get_graph();
	for (int v = 0; v < this->num_vertices(); v++) {
		g_copy.add_vertex(g[v].referenced_map, g[v].referenced_polys);
	}

    //iterate over all edges and insert them into the new graph
    auto es = edges(g);
    for (auto eit = es.first; eit != es.second; ++eit) {
        g_copy.add_edge(source(*eit, g), boost::target(*eit, g),g[*eit].weight);
    }


	//add all fix matches to the copied graph
	for (int v = 0; v < this->num_vertices(); v++) {
		if (this->getFixMatch(v) != -1) g_copy.setFixMatch(v, this->getFixMatch(v));
	}

	//return copied graph
	return g_copy;

}

int CandidateGraph::add_vertex(bool ref_map, std::vector<int> ref_polys) {
	//add vertex in every graph
	int new_vertex = boost::add_vertex(this->g);
    this->g[new_vertex].referenced_map = ref_map;
    this->g[new_vertex].referenced_polys = std::move(ref_polys);

    if (!ref_map) this->left_vs.push_back(new_vertex);
    else this->right_vs.push_back(new_vertex);


	//add memory field to inclusion lookup
	this->inclusion_lookup.push_back(std::vector<int> {});

	//add memory field to fix matches
	this->fix_matches.push_back(-1);

	//this is possible, since every graph has the same size and this every new vertex index is identical for the three graphs
	return new_vertex;
}

void CandidateGraph::delete_vertex(int v_id) {

    //delete all edges incident to v
    clear_vertex(v_id, this->g);
    //remove_vertex
    remove_vertex(v_id, this->g);


	//remove vertex ID from left and right vs (we don't know which list it was in previously)
	this->left_vs.erase(std::remove(this->left_vs.begin(), this->left_vs.end(), v_id), this->left_vs.end());
	this->right_vs.erase(std::remove(this->right_vs.begin(), this->right_vs.end(), v_id), this->right_vs.end());
	//each vertex ID bigger than the v_id of the deleted vertex needs to be decreases now
	for (auto& l : this->left_vs) {
		if (l > v_id) l--;
	}
	for (auto& r : this->right_vs) {
		if (r > v_id) r--;
	}

	//delete field from lookup table
	for (auto& inc_vec : this->inclusion_lookup) {
		//erase all occurences of v_id
		inc_vec.erase(remove(inc_vec.begin(), inc_vec.end(), v_id), inc_vec.end());
		for (auto& inc : inc_vec) {
			//if a reference is smaller than the deleted vertex ID, adapt it 
			if (inc > v_id) inc--;
		}
	}
	this->inclusion_lookup.erase(std::next(this->inclusion_lookup.begin(), v_id));

	//delete field from fix matches memory
	for (auto& fm: this->fix_matches) {
		if (fm == v_id) fm = -1;
		else if (fm > v_id) fm--;
	}

	this->fix_matches.erase(std::next(this->fix_matches.begin(), v_id));

}

void CandidateGraph::edit_vertex(int v_id, bool ref_map, std::vector<int> ref_polys) {
    this->g[v_id].referenced_map = ref_map;
    this->g[v_id].referenced_polys = ref_polys;
}

VertexProps CandidateGraph::get_vertex(int v_id) {
	return this->g[v_id];
}

int CandidateGraph::num_vertices() {
	return this->g.vertex_set().size();
}

int CandidateGraph::num_edges() {
	return boost::num_edges(this->g);
}

std::vector<int> CandidateGraph::get_Vi(bool map) {
	return !map? this->left_vs : this->right_vs;
}

void CandidateGraph::add_edge(int source, int target, double weight) {
	auto e = boost::add_edge(source, target, this->g);
    this->g[e.first].weight = weight;
}

void CandidateGraph::add_edge(int source, int target) {
	boost::add_edge(source, target, this->g);
}

void CandidateGraph::delete_edge(int source, int target) {
	auto source_vertex = vertex(source, this->g);
	auto target_vertex = vertex(target, this->g);

	remove_edge(source_vertex, target_vertex, this->g);
	remove_edge(target_vertex, source_vertex, this->g);
}

Graph* CandidateGraph::get_graph() {
	return &(this->g);
}

int CandidateGraph::get_components(std::vector<int>* components) {
	std::vector<int> components_local = std::vector<int>(boost::num_vertices(this->g));
	int num_components = connected_components(this->g, &components_local[0]);
	components = &components_local;
	return num_components;
}

bool CandidateGraph::are_connected(int source, int target) {
	return boost::edge(source, target, this->g).second;
}

int CandidateGraph::does_vertex_exist(bool map, const std::vector<int>& ref_polys) {
	Graph::vertex_iterator v, vend;
	for (boost::tie(v, vend) = vertices(this->g); v != vend; ++v) {
		if (this->g[0].referenced_map == map && this->g[0].referenced_polys == ref_polys) return *v;
	}
	return -1;

}

bool CandidateGraph::referenced_map(int v) {
	return this->g[v].referenced_map;
}

void CandidateGraph::addIncludedVertex(int v_id, int included_vertex) {
	this->inclusion_lookup[v_id].push_back(included_vertex);
}

void CandidateGraph::setIncludedVertices(int v_id, std::vector<int> included_vertices) {
	this->inclusion_lookup[v_id] = std::move(included_vertices);
}

std::vector<int> CandidateGraph::getIncludedVertices(int v_id) {
	return this->inclusion_lookup[v_id];
}

void CandidateGraph::setFixMatch(int v_1, int v_2) {
	this->fix_matches[v_1] = v_2;
	this->fix_matches[v_2] = v_1;
}

int CandidateGraph::getFixMatch(int v_id) {
	return this->fix_matches[v_id];
}

void CandidateGraph::printVertexOverview() {
	cout << "\n---CANDIDATE GRAPH VERTEX OVERVIEW---\n\n";
	//iterate over all vertices
	for (int v = 0; v < this->g.vertex_set().size(); v++) {
		std::string map_name = !this->g[v].referenced_map ? "1" : "2";
		cout << "vertex " << v << " represents "<<map_name << " polygons: {";
		for (const auto& ref_poly : this->g[v].referenced_polys) cout << ref_poly << ",";
		cout << "}" << endl;
	}

	cout << "\n---END CANDIDATE GRAPH VERTEX OVERVIEW---\n\n";
}

//END MULTILAYER GRAPH CLASS

void computeEdges(std::vector<CandidateGraph>& g_vec, Localization rtree1, Localization rtree2, const MapOverlay& mo, double lambda, const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2) {
    //for all graphs
    for (auto& cg : g_vec) {
        //add edges to graph
        //iterate over all pairs of vertices to check, if an edge in between them should be inserted
        Graph::vertex_iterator vL, vLend, vR, vRend;
        //remember all substituting edges in order to be able to add them if the vertices did not get connected via a path beforehand
        std::vector<std::tuple<int, int, bool>> sub_edges;

        //for referenced polygons, an arbitrary layer of the multilayer graph can be accessed, take layer 0
        Graph* g = cg.get_graph();
        for (boost::tie(vL, vLend) = vertices(*g); vL != vLend; ++vL) {
            for (boost::tie(vR, vRend) = vertices(*g); vR != vRend; ++vR) {
                //do not compute edges twice
                if (*vL > *vR) continue;
                //cout << "checking " << *vL << " and " << *vR << endl;

                //EDGES IN E3
                //vertices are of different maps, make sure vL has map 0 aand vR map 1 to not compute edges twice
                if ((*g)[*vL].referenced_map != (*g)[*vR].referenced_map) {
                    //pair of vertices where vL and vR are in different parts of the bipartite graph
                    std::vector<Polygon_wh> Lpolys = !(*g)[*vL].referenced_map ? polys1 : polys2;
                    std::vector<Polygon_wh> Rpolys = !(*g)[*vR].referenced_map ? polys1 : polys2;
                    Localization rtree = !(*g)[*vL].referenced_map ? rtree2 : rtree1;

                    int vertex_index1 = !(*g)[*vL].referenced_map ? *vL : *vR;
                    int vertex_index2 = !(*g)[*vL].referenced_map ? *vR : *vL;

                    //check if edge should be inserted at all, which is only the case if the sets of polygons overlap
                    if (mo.doOverlap((*g)[vertex_index1].referenced_polys, (*g)[vertex_index2].referenced_polys, 0.0)) {

                        //retreive edge weight IntersectionOverUnion from MapOverlay
                        double IoU = mo.getIoU((*g)[vertex_index1].referenced_polys, (*g)[vertex_index2].referenced_polys);

                        //consider lambda
                        double weight = IoU - lambda;

                        //only add the edge if at least one polygon has an inside/outside ratio that can improve a match above the lambda threshold
                        if (IoU > 0.0) cg.add_edge(*vL, *vR, weight);

                    }
                }
            }
        }
        //check if the vertices are already connected via path, if not add the direct edge to the graph
        for (const auto& sub_edge : sub_edges) {

            int source = std::get<0>(sub_edge);
            int target = std::get<1>(sub_edge);
            bool map = std::get<2>(sub_edge);



            std::vector<bool> visited; visited.resize(cg.num_vertices(), false);
            std::vector<int> source_component;
            DFS(cg.get_graph(), source, cg.num_vertices(), &visited, &source_component);

            if (std::find(source_component.begin(), source_component.end(), target) == source_component.end()) {
                //vertices are not connected yet, add the edge
                cg.add_edge(source, target, map);
            }

        }
    }
}

void DFS(Graph* g, int v, int upper_limit, std::vector<bool>* visited, std::vector<int>* component) {
	//mark v as visited
	(*visited)[v] = true;

	//add v to solution
	component->push_back(v);

	//get all neighbors of v
	auto neighbors = boost::adjacent_vertices(v, *g);

	for (auto n : make_iterator_range(neighbors)) {
		//call DFS recursively on the neighbor, if it is not visited yet and is legal to explore (if it is above upper limit it is a cumulated node, which should not be explored)
		if (!(*visited)[n] && v < upper_limit) DFS(g, n, upper_limit, visited, component);

	}

}

//DFS to collect all adjacent nodes up to a limited depth
void DFS_limited(Graph* g, int v, int upper_limit, int depth, int depth_limit, std::vector<bool>* visited, std::vector<int>* component) {
	//mark v as visited
	(*visited)[v] = true;

	//add v to solution
	component->push_back(v);

	//get all neighbors of v
	auto neighbors = boost::adjacent_vertices(v, *g);

	for (auto n : make_iterator_range(neighbors)) {
		//call DFS recursively on the neighbor, if it is not visited yet and is legal to explore (if it is above upper limit it is a cumulated node, which should not be explored)
		if (!(*visited)[n] && v < upper_limit && depth < depth_limit) DFS_limited(g, n, upper_limit,depth+1,depth_limit, visited, component);

	}

}

//checks if a set of vertices induces a connected subgraph of G, not consider vertices with index >= upper_limit
bool is_connected(std::vector<int> set, Graph g, int upper_limit) {
	//rule out edge case
	if (set.size() == 0) return false;

	//start DFS at an arbitrary node of the set
	std::vector<bool> visited;
	visited.resize(g.vertex_set().size(), false);
	std::vector<int> dfs_found_nodes;
	DFS(&g, set[0], upper_limit, &visited, &dfs_found_nodes);
	std::sort(dfs_found_nodes.begin(),dfs_found_nodes.end());
	std::vector<int> intersection;
	std::set_intersection(set.begin(), set.end(), dfs_found_nodes.begin(), dfs_found_nodes.end(),std::back_inserter(intersection));

	//if all nodes were found, the intersection has the same size as the input set, else at least one node has not been found via the DFS, meaning the set is not connected
	return set.size() == intersection.size();

}

template<typename T>
std::vector<std::vector<T>> getPowerset(std::vector<T> vec)
{
	int set_size = vec.size();
	std::vector<std::vector<T>> powerset;
	// Set_size of power set of a set with set_size
	// n is (2^n-1)
	unsigned int pow_set_size = pow(2, set_size);
	int counter, j;

	// Run from counter 000..0 to 111..1
	for (counter = 1; counter < pow_set_size; counter++) {
		std::vector<T> cur_set;
		for (j = 0; j < set_size; j++) {
			// Check if jth bit in the counter is set
			// If set then print jth element from set
			if (counter & (1 << j))
				cur_set.push_back(vec[j]);
		}
		powerset.push_back(cur_set);
	}

	return powerset;
}

//computes and returns the all subsets of vertices in vec, which are connected in cg
std::vector<std::vector<int>> getConnectedPowerset(std::vector<int> vec, CandidateGraph* cg, int upper_limit)
{
	//make sure, vector is sorted
	std::sort(vec.begin(),vec.end());

	//decide, which side of the bipartite graph the set is on
	bool map = cg->referenced_map(vec[0]);

	int set_size = vec.size();
	std::vector<std::vector<int>> powerset;
	// Set_size of power set of a set with set_size
	// n is (2^n-1)
	unsigned int pow_set_size = pow(2, set_size);
	int counter, j;

	// Run from counter 000..0 to 111..1
	for (counter = 1; counter < pow_set_size; counter++) {
		std::vector<int> cur_set;
		for (j = 0; j < set_size; j++) {
			// Check if jth bit in the counter is set
			// If set then print jth element from set
			if (counter & (1 << j))
				cur_set.push_back(vec[j]);
		}

		//check if the current set is connected in cg, if yes add it to the powerset
		Graph g = *(cg->get_graph());

		std::sort(cur_set.begin(), cur_set.end());
		assert(is_connected(cur_set, g, upper_limit) && "FOUND NON CONNECTED SET!");
		if(is_connected(cur_set, g, upper_limit)) powerset.push_back(cur_set);

	}

	return powerset;
}

bool compareVectors(const std::vector<int>& v1, const std::vector<int>& v2) {
	return std::lexicographical_compare(v1.begin(), v1.end(), v2.begin(), v2.end());
}


//REGION GROWING
struct queue_element {
	std::vector<int> included_vertices;
	std::vector<boost::graph_traits<Graph>::edge_descriptor> included_edges;
	std::vector< boost::graph_traits<Graph>::edge_descriptor> guessed_onetoone_matching;
	double guessed_onetoone_matching_quality, sum_of_edge_weights;

	// Constructor for queue_element
	queue_element(const std::vector<int>& vertices, const std::vector<boost::graph_traits<Graph>::edge_descriptor>& edges, const std::vector<boost::graph_traits<Graph>::edge_descriptor>& matching, double& quality, double& edge_weights)
		: included_vertices(vertices), included_edges(edges),guessed_onetoone_matching(matching), guessed_onetoone_matching_quality(quality), sum_of_edge_weights(edge_weights) {}

	//comparator for queue_element
	inline bool operator==(queue_element q) {
		//two queue elements are equal, if they contain the same vertices
		std::sort(this->included_vertices.begin(), this->included_vertices.end());
		std::sort(q.included_vertices.begin(), q.included_vertices.end());
		return this->included_vertices == q.included_vertices;
	}

	inline bool operator<(queue_element q) {
		//first criterion: sort by size of included edge set:
		if (this->included_edges.size() < q.included_edges.size()) return true;
		else if (this->included_edges.size() > q.included_edges.size()) return false;
		else {
			if (this->sum_of_edge_weights < q.sum_of_edge_weights) return true;
			else return false;
		}
	}

	inline bool operator>(queue_element q) {
		return !(*this < q);
	}

	//for printing
	friend std::ostream& operator<<(std::ostream& os, const queue_element& obj) {
		os << "QUEUE_ELEMENT_VERTICES:{";
		for (const auto& v : obj.included_vertices) os << v << ",";
		os << "};";
		return os;
	}
};

//an object representing a component consisting of vertex indices of included vertices
struct component {
	//for comparison we assume the ids are always sorted!
	std::vector<int> represented_vertex_indices;

	component(const std::vector<int>& ids)
		: represented_vertex_indices(ids) {}

	inline bool operator==(component c) {
		if (this->represented_vertex_indices.size() == c.represented_vertex_indices.size()) {
			return this->represented_vertex_indices == c.represented_vertex_indices;
		}
		return false;
	}

	inline bool operator<(component c) {
		if (this->represented_vertex_indices.size() < c.represented_vertex_indices.size()) {
			return true;
		}
		else if (this->represented_vertex_indices.size() == c.represented_vertex_indices.size()) {
			return std::lexicographical_compare(this->represented_vertex_indices.begin(), this->represented_vertex_indices.end(), c.represented_vertex_indices.begin(), c.represented_vertex_indices.end());
		}
		return false;
	}

	inline bool operator>(component c) {
		return !(*this < c);
	}
};

//an object representing a cumulative node containing its vertex descriptor and an array of indices of vertices represented by it
struct cumulative_node {
	std::vector<int> represented_vertex_indices;
	Graph::vertex_descriptor vertex_descriptor;
	std::size_t hash;

	cumulative_node(const std::vector<int>& rep_ids, const Graph::vertex_descriptor& desc)
		: represented_vertex_indices(rep_ids), vertex_descriptor(desc), hash(0) {}

	inline bool operator== (cumulative_node c) {
		//two cumulative nodes are equal iff they represent the same vertex set, this way they can get compared to sets of vertices without a corresponding descriptor as well
		if (this->represented_vertex_indices.size() == c.represented_vertex_indices.size()) {
			return this->represented_vertex_indices == c.represented_vertex_indices;
		}
		return false;
	}

	inline bool operator< (cumulative_node c) {
		if (this->represented_vertex_indices.size() < c.represented_vertex_indices.size()) {
			return true;
		}
		else if (this->represented_vertex_indices.size() == c.represented_vertex_indices.size()) {
			return std::lexicographical_compare(this->represented_vertex_indices.begin(), this->represented_vertex_indices.end(), c.represented_vertex_indices.begin(), c.represented_vertex_indices.end());
		}
		return false;
		//return this->hash < c.hash;
	}

	inline bool operator> (cumulative_node c) {
		return !(*this < c);
	}

};

//define struct for memory

struct queue_element_vb {
	std::vector<int> vertices;
	std::vector<int> vertex_mem;
	std::vector<boost::graph_traits<Graph>::edge_descriptor> edge_mem;
	//store sum of weights of included edges as unique key
	double sum_of_edge_weights;
	std::size_t hash;

	// Constructor for queue_element
	queue_element_vb(const std::vector<int>& vertices, const std::vector<int>& vertices_mem, const std::vector<boost::graph_traits<Graph>::edge_descriptor>& edges, const double& weights)
		: vertices(vertices), vertex_mem(vertices_mem), edge_mem(edges), sum_of_edge_weights(weights), hash(0) {}

	//comparator for queue_element
	inline bool operator==(queue_element_vb q) {
		//two queue elements are equal, if they contain the same vertices
		//std::sort(this->vertices.begin(), this->vertices.end());
		//std::sort(q.vertices.begin(), q.vertices.end());
		return this->vertices == q.vertices && this->vertex_mem == q.vertex_mem;
	}

	inline bool operator<(queue_element_vb q) {
		//first criterion: sort by size of vertex memory
		if (this->vertex_mem.size() < q.vertex_mem.size()) return true;
		else {
			if (this->sum_of_edge_weights < q.sum_of_edge_weights) return true;
			else return false;
		}
		//compare the hashes instead
		//return this->hash < q.hash;
	}

	inline bool operator>(queue_element_vb q) {
		return !(*this < q);
	}

	//for printing
	friend std::ostream& operator<<(std::ostream& os, const queue_element_vb& obj) {
		os << "QUEUE_ELEMENT_VERTICES:{";
		for (const auto& v : obj.vertex_mem) os << v << ",";
		os << " || ";
		for (const auto& v : obj.vertices) os << v << ",";
		os << "};";
		return os;
	}
};

//BINARY SEARCH TREE DATA STRUCTURE
//Define a data structure for keeping track of already found sets
// Given Node
template <typename T>
struct TreeNode {
	T key;
	TreeNode* left;
	TreeNode* right;
	int height; // height of the subtree rooted at this node

	// Constructor
	TreeNode(T item) : key(item), left(nullptr), right(nullptr), height(1) {}
};

// Function to create a new BST node
template <typename T>
TreeNode<T>* newNode(T item) {
	return new TreeNode<T>(item);
}

// Function to perform a right rotation
template <typename T>
TreeNode<T>* rotateRight(TreeNode<T>* y) {
	TreeNode<T>* x = y->left;
	TreeNode<T>* T2 = x->right;

	// Perform rotation
	x->right = y;
	y->left = T2;

	// Update heights
	updateHeight(y);
	updateHeight(x);

	return x; // New root
}

// Function to perform a left rotation
template <typename T>
TreeNode<T>* rotateLeft(TreeNode<T>* x) {
	TreeNode<T>* y = x->right;
	TreeNode<T>* T2 = y->left;

	// Perform rotation
	y->left = x;
	x->right = T2;

	// Update heights
	updateHeight(x);
	updateHeight(y);

	return y; // New root
}

// Function to insert a new node with given key in AVL tree
template <typename T>
TreeNode<T>* insert(TreeNode<T>* node, T key) {
	// If the tree is empty, return a new node
	if (node == nullptr)
		return newNode(key);

	// Otherwise, recur down the tree
	if (key < node->key)
		node->left = insert(node->left, key);
	else if (key > node->key)
		node->right = insert(node->right, key);
	else // Duplicate keys not allowed
		return node;

	// Update height of this ancestor node
	updateHeight(node);

	// Get the balance factor of this ancestor node to check whether this node became unbalanced
	int balance = getBalanceFactor(node);

	// If this node becomes unbalanced, there are four cases

	// Left Left Case
	if (balance > 1 && key < node->left->key)
		return rotateRight(node);

	// Right Right Case
	if (balance < -1 && key > node->right->key)
		return rotateLeft(node);

	// Left Right Case
	if (balance > 1 && key > node->left->key) {
		node->left = rotateLeft(node->left);
		return rotateRight(node);
	}

	// Right Left Case
	if (balance < -1 && key < node->right->key) {
		node->right = rotateRight(node->right);
		return rotateLeft(node);
	}

	// Return the (unchanged) node pointer
	return node;
}

// Function to search for a key in the BST
template <typename T>
bool search(TreeNode<T>* node, T key) {
	// If the tree is empty or the key is not found, return false
	if (node == nullptr)
		return false;

	// If the key is found, return true
	if (node->key == key)
		return true;

	// Recur down the tree
	if (key < node->key)
		return search(node->left, key);
	else
		return search(node->right, key);
}

// Function to search for a key in the BST and return it
template <typename T>
T* retreive(TreeNode<T>* node, T key) {
	// If the tree is empty or the key is not found, return false
	if (node == nullptr)
		return nullptr;

	// If the key is found, return true
	if (node->key == key)
		return &(node->key);

	// Recur down the tree
	if (key < node->key)
		return retreive(node->left, key);
	else
		return retreive(node->right, key);
}

// Function to find the minimum value node in a given BST
template <typename T>
TreeNode<T>* minValueNode(TreeNode<T>* node) {
	TreeNode<T>* current = node;
	while (current && current->left != nullptr)
		current = current->left;
	return current;
}

// Function to remove a node with the given key from AVL tree
template <typename T>
TreeNode<T>* remove(TreeNode<T>* root, T key) {
	// Base case: If the tree is empty, return
	if (root == nullptr)
		return root;

	// Perform standard BST delete
	if (key < root->key)
		root->left = remove(root->left, key);
	else if (key > root->key)
		root->right = remove(root->right, key);
	else {
		// Node with only one child or no child
		if (root->left == nullptr || root->right == nullptr) {
			TreeNode<T>* temp = root->left ? root->left : root->right;

			// No child case
			if (temp == nullptr) {
				temp = root;
				root = nullptr;
			}
			else { // One child case
				*root = *temp; // Copy the contents of the non-empty child
			}
			delete temp;
		}
		else {
			// Node with two children: Get the inorder successor (smallest in the right subtree)
			TreeNode<T>* temp = minValueNode(root->right);

			// Copy the inorder successor's content to this node
			root->key = temp->key;

			// Delete the inorder successor
			root->right = remove(root->right, temp->key);
		}
	}

	// If the tree had only one node, then return
	if (root == nullptr)
		return root;

	// Update the height of the current node
	updateHeight(root);

	// Get the balance factor of this node to check whether this node became unbalanced
	int balance = getBalanceFactor(root);

	// If this node becomes unbalanced, there are four cases

	// Left Left Case
	if (balance > 1 && getBalanceFactor(root->left) >= 0)
		return rotateRight(root);

	// Left Right Case
	if (balance > 1 && getBalanceFactor(root->left) < 0) {
		root->left = rotateLeft(root->left);
		return rotateRight(root);
	}

	// Right Right Case
	if (balance < -1 && getBalanceFactor(root->right) <= 0)
		return rotateLeft(root);

	// Right Left Case
	if (balance < -1 && getBalanceFactor(root->right) > 0) {
		root->right = rotateRight(root->right);
		return rotateLeft(root);
	}

	// Return the (unchanged) node pointer
	return root;
}

// Function to remove the root node of the BST
template <typename T>
TreeNode<T>* removeRoot(TreeNode<T>* root) {
	return remove(root, root->key);
}

// Function to extract the key of the root node and delete it
template <typename T>
T extractRootKey(TreeNode<T>*& root) {
	auto root_key = root->key;
	removeRoot(root);
	return root_key;
}

// Function to get the total size of the BST (number of nodes)
template <typename T>
std::size_t getSize(TreeNode<T>* root) {
	if (root == nullptr) {
		return 0;
	}

	// Recursively count the nodes in the left and right subtrees
	std::size_t leftSize = getSize(root->left);
	std::size_t rightSize = getSize(root->right);

	// Return the total size of the tree (sum of nodes in left and right subtrees + 1 for the root)
	return leftSize + rightSize + 1;
}

// Function to traverse the BST and print its elements (inorder traversal)
template <typename T>
void printTree(TreeNode<T>* root) {
	if (root == nullptr) {
		return;
	}

	

	// Print current node's key
	std::cout << root->key << " ";

	// Traverse left subtree
	printTree(root->left);

	// Traverse right subtree
	printTree(root->right);
}

// Function to get the size of a subtree
template <typename T>
std::size_t getSubtreeSize(TreeNode<T>* node) {
	if (node == nullptr) {
		return 0;
	}
	// Recursively count the nodes in the left and right subtrees
	std::size_t leftSize = getSubtreeSize(node->left);
	std::size_t rightSize = getSubtreeSize(node->right);
	// Return the total size of the subtree (sum of nodes in left and right subtrees + 1 for the root)
	return leftSize + rightSize + 1;
}

// Function to print the size of the left and right subtrees separately
template <typename T>
void printSubtreeSizes(TreeNode<T>* root) {
	if (root == nullptr) {
		std::cout << "Left Subtree Size: 0\n";
		std::cout << "Right Subtree Size: 0\n";
		return;
	}

	std::size_t leftSize = getSubtreeSize(root->left);
	std::size_t rightSize = getSubtreeSize(root->right);

	std::cout << "Left Subtree Size: " << leftSize << std::endl;
	std::cout << "Right Subtree Size: " << rightSize << std::endl;
}

template <typename T>
int getBalanceFactor(TreeNode<T>* node) {
	if (node == nullptr) {
		return 0;
	}
	return getHeight(node->left) - getHeight(node->right);
}

// Function to get the height of a node
template <typename T>
int getHeight(TreeNode<T>* node) {
	if (node == nullptr) {
		return 0;
	}
	return node->height;
}

// Function to update the height of a node
template <typename T>
void updateHeight(TreeNode<T>* node) {
	if (node == nullptr) {
		return;
	}
	node->height = 1 + std::max(getHeight(node->left), getHeight(node->right));
}

void VertexBasedSetCollection1to1Guessing(CandidateGraph* cg, double lambda, const MapOverlay& mo, std::chrono::steady_clock::time_point start_time, int time_limit, int size_limit) {
	Graph g_2 = *cg->get_graph();

	//make sure not to explore vertices, that have been inserted by this function
	int upper_limit = cg->num_vertices();

	//create local memory of all components (=subsets) that have been found already, comps_all serves as a global memory, that stores if a certain set was already queued
	TreeNode<cumulative_node>* comps_all = nullptr, * comps0 = nullptr, * comps1 = nullptr;

	//start out at every vertex once
	Graph::vertex_iterator v, vend;
	for (tie(v, vend) = vertices(g_2); v != vend; v++) {
		//make sure not to explore backwards
		int lower_limit = *v;

		//TreeNode<queue_element_vb>* queue = nullptr;
		std::queue<queue_element_vb> queue;

		//queue = insert(queue, queue_element_vb(std::vector<int>{(int)*v}, std::vector<int>{(int)*v}, std::vector<Graph::edge_descriptor> {}, 0.0));
		queue.push(queue_element_vb(std::vector<int>{(int)*v}, std::vector<int>{(int)*v}, std::vector<Graph::edge_descriptor> {}, 0.0));


		//check if the current vertex is already a found component, else mark it as one such that it does not get visited again
		if (comps_all == nullptr || !search(comps_all, cumulative_node(std::vector<int>{(int)*v}, * v))) {
			//v is not yet in components, insert in a sorted way

			//the set is already included in the graph as a cumulated vertex, load its id
			//insert in sorted way into comps
			comps_all = insert(comps_all, cumulative_node(std::vector<int>{(int)*v}, * v));
		}

		while (!queue.empty()) {
			//kill exploration if necessary
			if (std::chrono::steady_clock::now() - start_time > std::chrono::seconds(time_limit)) return;

			//dequeue node
			auto queue_element = queue.front();
			queue.pop();

			auto current_node_set = queue_element.vertices;

			//get all neighbors of v
			std::vector<int> neighbors;

			for (const auto& node : current_node_set) {
				auto neighbors_boost = adjacent_vertices(node, g_2);

				for (auto n : make_iterator_range(neighbors_boost)) {
					//remember neighbor if it is not already a cumulated node or in the vertex memory
					if ((int)n > lower_limit && (int)n < upper_limit && !std::binary_search(queue_element.vertex_mem.begin(), queue_element.vertex_mem.end(), n)) neighbors.push_back(n);
				}
			}
			//make neighbor list unique
			std::sort(neighbors.begin(), neighbors.end());
			neighbors.erase(unique(neighbors.begin(), neighbors.end()), neighbors.end());

			//those neighbors, that allow for a 1:1 matching to be formed, that exceeds 1- lambda in quality, do not need to be considered, as this set
			//would never be in a m:n match

			//remember all edges to a neighbor n, if it should be explored further we need those for the edge memory
			std::vector<Graph::edge_descriptor> new_edges;
			for (int n_i = 0; n_i < neighbors.size(); n_i++) {
				int n = neighbors[n_i];
				//load all of the edges connecting the current node set with n
				Graph::vertex_descriptor n_desc = vertex(n, g_2);
				std::vector<Graph::edge_descriptor> edges_to_n;
				for (const auto& cur : current_node_set) {
					Graph::vertex_descriptor cur_desc = vertex(cur, g_2);
					//if an edge between the vertices exists, remember it
					if (edge(cur_desc, n_desc, g_2).second) {
						edges_to_n.push_back(edge(cur_desc, n_desc, g_2).first);
					}
				}
				//remember, if there is any edge that allows ruling out n for further exploration
				bool n_should_be_explored = true;

				//we guess the 1:1 matching as follows: first, we collect the set of edges and sort it w.r.t the weights
				std::vector<Graph::edge_descriptor>  all_edges = edges_to_n;
				for (const auto& e_mem_pair : queue_element.edge_mem) all_edges.push_back(edge(source(e_mem_pair, g_2), target(e_mem_pair, g_2), *(cg->get_graph())).first);

				std::vector<double> edge_weights;
				//while collecting weights, directly remember max weight, if this is less than 0, we do not need to guess as then the sum of 1:1 match qualities can never be positive
				double max_weight = std::numeric_limits<double>::min();
				for (const auto& ae : all_edges) {
					double w = g_2[ae].weight;
					edge_weights.push_back(w);

					if (w > max_weight) max_weight = w;
				}

				if (max_weight > 0.0) {
					std::vector<int> edge_indices(edge_weights.size());
					std::iota(edge_indices.begin(), edge_indices.end(), 0);
					std::sort(edge_indices.begin(), edge_indices.end(), [&](int A, int B) -> bool {
						return edge_weights[A] < edge_weights[B];
						});

					//the edge_indices vector now contains the indices in the sorted order w.r.t. the edge weights
					//we will traverse the edges in descending order and try to add as many good edges to the solution as possible
					//remember, which vertices are already contained in the solution, no vertex should be added twice as it is a 1:1 matching
					std::vector<int> contained_vertices;
					double matching_quality = 0.0;

					for (int i = edge_weights.size() - 1; i >= 0; i--) {
						if (edge_weights[edge_indices[i]] < 0.0) break;
						int s = source(all_edges[edge_indices[i]], g_2), t = target(all_edges[edge_indices[i]], g_2);
						if (std::find(contained_vertices.begin(), contained_vertices.end(), s) == contained_vertices.end() && std::find(contained_vertices.begin(), contained_vertices.end(), t) == contained_vertices.end()) {
							//edge is legal and should be added to the solution
							contained_vertices.push_back(s);
							contained_vertices.push_back(t);
							matching_quality += edge_weights[edge_indices[i]];
						}
					}
					if (matching_quality > 1 - lambda) {
						n_should_be_explored = false;
					}

				}


				if (!n_should_be_explored) {
					//remove n from neighbors
					neighbors.erase(remove(neighbors.begin(), neighbors.end(), n), neighbors.end());
					n_i--;
				}
				else {
					//remember new edges
					new_edges.insert(new_edges.end(), edges_to_n.begin(), edges_to_n.end());
				}
			}


			//only proceed, if new neighbors have been found
			if (neighbors.size() != 0) {
                //compute powerset of neighbors
				std::vector<std::vector<int>> powerset = getPowerset<int>(neighbors);

				int queue_element_size = queue_element.vertex_mem.size();
				std::vector<std::vector<int>> powerset_filtered;
				for (const auto& p : powerset) {
                    // if size limitation is set, do not explore further when it is exceeded
					if (size_limit == -1 || queue_element_size + p.size() <= size_limit) powerset_filtered.push_back(p);
				}
				powerset = powerset_filtered;

				//some elements of the powerset can be ruled out, if the polygon has a fixed match, of which we know it has to be contained in the same match in any optimal solution
				std::vector<int> fix_matches;
				for (auto& node : current_node_set) {
					int fm = cg->getFixMatch(node);
                    //filter the fix matches to those not already satisfied by the vertex memory
					if (fm > -1 && std::find(queue_element.vertex_mem.begin(),queue_element.vertex_mem.end(),fm) == queue_element.vertex_mem.end()) fix_matches.push_back(fm);
				}
				std::vector<std::vector<int>> reduced_powerset;
				if (fix_matches.size() > 0) {
					std::sort(fix_matches.begin(), fix_matches.end());

					for (auto& set : powerset) {
                        //sort set for later processing
                        std::sort(set.begin(), set.end());

						std::vector<int> set_inter;
						std::set_intersection(fix_matches.begin(), fix_matches.end(), set.begin(), set.end(), std::back_inserter(set_inter));
						if (set_inter.size() > 0) reduced_powerset.push_back(set);
					}
				}
				else reduced_powerset = powerset;


				//add every element within the powerset combined with the path memory to the solution

				for (auto& set : reduced_powerset) {

					std::sort(set.begin(), set.end());

					std::vector<int> new_set;
					new_set.insert(new_set.end(), queue_element.vertex_mem.begin(), queue_element.vertex_mem.end());
					new_set.insert(new_set.end(), set.begin(), set.end());

					std::sort(new_set.begin(), new_set.end());

					auto it = std::unique(new_set.begin(), new_set.end());
					new_set.resize(std::distance(new_set.begin(), it));

					//if this new set has already been visited, we do not need to check for existing vertices or edges, as we can be certain they exist
					if (comps_all == nullptr || !search(comps_all, cumulative_node(new_set, 0))) {


						//split up new set into sets for each side
						std::vector<std::vector<int>> new_sets_split(2);
						for (const auto& n : new_set) {
							new_sets_split[g_2[n].referenced_map].push_back(n);
						}

						std::sort(new_sets_split[0].begin(), new_sets_split[0].end());
						std::sort(new_sets_split[1].begin(), new_sets_split[1].end());


						//remember vertex indices of new sets
						std::vector<int> ids = { new_sets_split[0][0],new_sets_split[1][0] };
						std::vector<std::vector<int>> rep_polys(2);

						//note, if new set was inserted (both sets could have size 1)
						bool new_set_inserted = false;
						for (int l = 0; l < 2; l++) {
							TreeNode<cumulative_node>** comps = l == 0 ? &comps0 : &comps1;
							//only push back if not yet included
							if (new_sets_split[l].size() > 1 && (comps == nullptr || !search(*comps, cumulative_node(new_sets_split[l], 0)))) {
								//found a new set, this should be added to the graph
								std::vector<int> rep_polys_of_new_set;
								for (const auto& ns : new_sets_split[l]) rep_polys_of_new_set.insert(rep_polys_of_new_set.end(), g_2[ns].referenced_polys.begin(), g_2[ns].referenced_polys.end());
								std::sort(rep_polys_of_new_set.begin(), rep_polys_of_new_set.end());
								ids[l] = cg->add_vertex(l, rep_polys_of_new_set);
								//reset local graph
								g_2 = *cg->get_graph();
								rep_polys[l] = rep_polys_of_new_set;

								//insert in sorted way into comps
								*comps = insert(*comps, cumulative_node(new_sets_split[l], ids[l]));

								new_set_inserted = true;

							}
							else if (new_sets_split[l].size() == 1) {
								//the set includes only one vertex, collect rep polys
								rep_polys[l] = g_2[new_sets_split[l][0]].referenced_polys;

							}
							else if (comps != nullptr && search(*comps, cumulative_node(new_sets_split[l], 0))) {
								//the set is already included in the graph as a cumulated vertex, load its id
								cumulative_node* c_node = retreive(*comps, cumulative_node(new_sets_split[l], 0));

								// value should have been found
								if (c_node != nullptr) {
									//set id to vertex index
									ids[l] = c_node->vertex_descriptor;
									//set polys
									rep_polys[l] = g_2[ids[l]].referenced_polys;

								}
								else {
									cout << "weird error at finding an element that should have been in components." << endl;
								}


							}
						}

						//insert to comps
						comps_all = insert(comps_all, cumulative_node(new_set, 0));

						//prepare all edges to insert
						std::vector<Graph::edge_descriptor> all_new_edges;
						all_new_edges.insert(all_new_edges.end(), queue_element.edge_mem.begin(), queue_element.edge_mem.end());
						all_new_edges.insert(all_new_edges.end(), new_edges.begin(), new_edges.end());
						sort(all_new_edges.begin(), all_new_edges.end());


						//compute the sum of edge weights as key
						double sum_of_edge_weights = 0.0;
						for (const auto& e : all_new_edges) {
							sum_of_edge_weights += g_2[e].weight;
						}

						//push to queue
						queue.push(queue_element_vb(set, new_set, all_new_edges, sum_of_edge_weights));

						//insert edge if it does not exists yet
						if (!edge(ids[0], ids[1], g_2).second) {
							//insert the edge between the two sets
							double IoU = !g_2[ids[0]].referenced_map ? mo.getIoU(rep_polys[0], rep_polys[1]) : mo.getIoU(rep_polys[1], rep_polys[0]);
							cg->add_edge(ids[0], ids[1], IoU - lambda);
						}
					}
				}

			}

		}

	}
}

void completeCandidateGraph(CandidateGraph* cg, const MapOverlay& mo, double lambda, std::chrono::steady_clock::time_point start_time, int time_limit, int size_limit) {
	VertexBasedSetCollection1to1Guessing(cg, lambda, mo, start_time, time_limit, size_limit);
}

void precomputeSimpleMatches(CandidateGraph* cg, const MapOverlay& mo, double lambda, Solution* sol) {
	Graph g = *(cg->get_graph());
	//store number of vertices in side of bipartite graph, this might change
	int num_vertices = cg->get_Vi(0).size();

	//initialize tracker, if the computation did change the graph, it might have to be executed multiple times since
	//deletions may open up possibilities for further simple matches
	int graph_size_before_loop = cg->num_vertices();
	
	//iterate over both sets of vertices in the bipartite graph as starting points
	do {
		graph_size_before_loop = cg->num_vertices();
		for (bool map : {0, 1}) {
			num_vertices = cg->get_Vi(map).size();
			for (int i = 0; i < num_vertices; i++) {
				g = *(cg->get_graph());
				//guess 1:1 matching in the following way: we fix the current polygon to be matched to its best compatible partner q
				int p = cg->get_Vi(map)[i];

				//find best partner
				double max_weight = -1.0;
				auto outg_edges = g.out_edge_list(p);

				if (g.out_edge_list(p).size() == 0) continue;

				//while checking for the most expensive incident edge, also check if the referenced polygon(s) might have such little oversection with every neighbor,
				//that the vertex can be deleted
				std::vector<int> neighbor_ref_polys;

				Graph::edge_descriptor ei_max;

                for (auto ei = out_edges(p,g); ei.first != ei.second; ++ei.first) {
                    auto e = *ei.first;
					//remember referenced polygons of neighbor
					for (const auto& np : g[e.m_target].referenced_polys) neighbor_ref_polys.push_back(np);

					double e_weight = g[e].weight;
					if (e_weight > max_weight) {
						max_weight = e_weight;
						ei_max = e;
					}

				}


				//compute intersection area with every neighbor
				double intersection_area = 0.0;
				if (!map) {
					intersection_area += mo.getIntersectionArea(g[p].referenced_polys, neighbor_ref_polys);
				}
				else {
					intersection_area += mo.getIntersectionArea(neighbor_ref_polys, g[p].referenced_polys);
				}


				//if intersection area in relation to not intersected area is small, the vertex can be deleted as it will never be contained in any match
				if (intersection_area / (mo.getArea(map, g[p].referenced_polys) - intersection_area) < lambda) {
						 
					cg->delete_vertex(p);

					//update amount of vertices
					num_vertices = cg->get_Vi(map).size();

					//set p taking deletions into account
					i--;

					//skip to next iteration
					continue;
				}

				int q = ei_max.m_target;
				std::pair<int, int> m = { p, q };
				double IoU_m = !map ? mo.getIoU(g[p].referenced_polys, g[q].referenced_polys) : mo.getIoU(g[q].referenced_polys, g[p].referenced_polys);

				//if the match is worse than 0.0, we can skip p,q as a potential pair right away
				if (IoU_m < lambda) continue;

				//if the three properties are fulfilled, we know that p,q are in the same match in any optimal solution and can flag them accordingly:
				// 1. IoU(m) > lambda | 2. I(p,q) / area(p\q) | 3. I(p, other neighbors) / area (p \other neighbors) < lambda
				double I_area_pq = !map ? mo.getIntersectionArea(g[p].referenced_polys, g[q].referenced_polys) : mo.getIntersectionArea(g[q].referenced_polys, g[p].referenced_polys);
				double area_p = mo.getArea(map, g[p].referenced_polys);
				//intersection area with every neighbor has been computed before, the Intersection of p and q only has to be substracted from that
				double p_inter_with_other_neighbors = intersection_area - I_area_pq;
				if (IoU_m > lambda && I_area_pq > area_p - I_area_pq && p_inter_with_other_neighbors / (p - p_inter_with_other_neighbors) < lambda) {
					cg->setFixMatch(p, q);
				}


				//from there on, for each other polygon, which is either connected to p or q in the graph, we find its max available neighbor and match it to it
				//collect all polygons, that need to be considered
				std::vector<int> neighbors_p, neighbors_q;

				Graph::adjacency_iterator ai, ai_end;
				for (tie(ai, ai_end) = adjacent_vertices(p, g); ai != ai_end; ++ai) {
					if (*ai != q) neighbors_p.push_back(*ai);
				}
				for (tie(ai, ai_end) = adjacent_vertices(q, g); ai != ai_end; ++ai) {
					if (*ai != p) neighbors_q.push_back(*ai);
				}

				//create lists to remember which polygons are already covered in a match
				std::vector<int> considered;
				considered.resize(g.vertex_set().size(), false);
				considered[p] = true; considered[q] = true;

				//indicator: is property fulfilled and match p,q can be deleted from the graph?
				bool pq_is_optimal_match = true;

				//first iterate through all p polygons and assign the best still unconsidered q polygon
				for (int p_prime : neighbors_q) {
					//find best partner
					double max_weight = -1.0;
					Graph::out_edge_iterator ei, ei_end, ei_max;
					//collect all other partners as well
					std::vector<Graph::out_edge_iterator> ei_other;
					for (tie(ei, ei_end) = out_edges(p_prime, g); ei != ei_end; ++ei) {
						double e_weight = g[*ei].weight;
						if (e_weight > max_weight && target(*ei, g) != q) {
							max_weight = e_weight;
							ei_max = ei;
						}

						ei_other.push_back(ei);
					}

					//remove max edge from ei_other
					ei_other.erase(remove(ei_other.begin(), ei_other.end(), ei_max), ei_other.end());

					//check if partner has been found
					std::pair<int, int> m_prime;
					if (max_weight > 0.0) {
						int q_prime = target(*ei_max, g);
						m_prime = { p_prime, q_prime };


						double IoU_m_prime = !map ? mo.getIoU(g[p_prime].referenced_polys, g[q_prime].referenced_polys) : mo.getIoU(g[q_prime].referenced_polys, g[p_prime].referenced_polys);
						double IoU_sum = !map ? mo.getIoU(g[p].referenced_polys, g[q].referenced_polys) + IoU_m_prime : mo.getIoU(g[q].referenced_polys, g[p].referenced_polys) + IoU_m_prime;

					

						//also check, that the second needed Lemma actually holds ( I(p_prime,q_prime) > A(p_prime \ q_prime)  |  IoU(p_prime,prime) > lambda | I(p_prime,q_stars) / A(p_prime \ q_stars) < lambda)
						double I = !map ? mo.getIntersectionArea(g[p_prime].referenced_polys, g[q_prime].referenced_polys) : mo.getIntersectionArea(g[q_prime].referenced_polys, g[p_prime].referenced_polys);
						double A_without_q_prime = mo.getArea(map, g[p_prime].referenced_polys) - I;

						//collect intersection with every other neighbor
						std::vector<int> ei_other_ref_polys;
						for (const auto& eio : ei_other) ei_other_ref_polys.insert(ei_other_ref_polys.end(), g[target(*eio, g)].referenced_polys.begin(), g[target(*eio, g)].referenced_polys.end());
						double I_with_other_neighbors = !map ? mo.getIntersectionArea(g[p_prime].referenced_polys, ei_other_ref_polys) : mo.getIntersectionArea(ei_other_ref_polys, g[p_prime].referenced_polys);
						double A_without_other_neighbors = mo.getArea(map, g[p_prime].referenced_polys) - I_with_other_neighbors;

						//COMPUTE IOU OF BIG MATCH INCLUDING ALL 4 POLYGONS
						std::vector<int> ref_polys_p, ref_polys_q;
						ref_polys_p.insert(ref_polys_p.end(), g[p].referenced_polys.begin(), g[p].referenced_polys.end());
						ref_polys_p.insert(ref_polys_p.end(), g[p_prime].referenced_polys.begin(), g[p_prime].referenced_polys.end());
						ref_polys_q.insert(ref_polys_q.end(), g[q].referenced_polys.begin(), g[q].referenced_polys.end());
						ref_polys_q.insert(ref_polys_q.end(), g[q_prime].referenced_polys.begin(), g[q_prime].referenced_polys.end());

						//double IoU_all = !map ? mo.getIoU(ref_polys_p, ref_polys_q) : mo.getIoU(ref_polys_q, ref_polys_p);
						//END COMPUTING BIG MATCH IOU

						if (1 - lambda > IoU_sum - 2 * lambda || I <= A_without_q_prime || IoU_m_prime <= lambda || I_with_other_neighbors / A_without_other_neighbors >= lambda) {
						//if (1 - lambda > IoU_sum - 2 * lambda) {
						//if this happens, the bigger match would actually improve upon the singular matches and our property is hurt
							pq_is_optimal_match = false;
							continue;
						}
						else {
							//in this case, forming the two matches {p,q} and {p',q'} is actually better than one match could ever be, therefore we know that {p,q'} 
							//and {q',p} can never be in the same match and we can delete the edges from the graph
							cg->delete_edge(p, q_prime);
							cg->delete_edge(q, p_prime);

						}
					}
					else {
						pq_is_optimal_match = false;
						continue;
					}
				}
				if (!pq_is_optimal_match) continue;

				//the same needs to be done for all q polygons
				for (int q_prime : neighbors_p) {
					//find best partner
					double max_weight = -1.0;
					Graph::out_edge_iterator ei, ei_end, ei_max;
					std::vector<Graph::out_edge_iterator> ei_other;
					for (tie(ei, ei_end) = out_edges(q_prime, g); ei != ei_end; ++ei) {
						double e_weight = g[*ei].weight;
						if (e_weight > max_weight && target(*ei, g) != p) {
							max_weight = e_weight;
							ei_max = ei;
						}
						ei_other.push_back(ei);
					}


					//remove max edge from ei_other
					ei_other.erase(remove(ei_other.begin(), ei_other.end(), ei_max), ei_other.end());

					//check if partner has been found
					std::pair<int, int> m_prime;
					if (max_weight > 0.0) {
						int p_prime = target(*ei_max, g);
						m_prime = { p_prime, q_prime };

						double IoU_m_prime = !map ? mo.getIoU(g[p_prime].referenced_polys, g[q_prime].referenced_polys) : mo.getIoU(g[q_prime].referenced_polys, g[p_prime].referenced_polys);
						double IoU_sum = !map ? mo.getIoU(g[p].referenced_polys, g[q].referenced_polys) + IoU_m_prime : mo.getIoU(g[q].referenced_polys, g[p].referenced_polys) + IoU_m_prime;

						//also check, that the second needed Lemma actually holds ( I(p_prime,q_prime) > A(p_prime \ q_prime)  |  IoU(p_prime,prime) > lambda | I(p_prime,q_stars) / A(p_prime \ q_stars) < lambda)
						double I = !map ? mo.getIntersectionArea(g[p_prime].referenced_polys, g[q_prime].referenced_polys) : mo.getIntersectionArea(g[q_prime].referenced_polys, g[p_prime].referenced_polys);
						double A_without_p_prime = mo.getArea(map, g[p_prime].referenced_polys) - I;

						//collect intersection with every other neighbor
						std::vector<int> ei_other_ref_polys;
						for (const auto& eio : ei_other) ei_other_ref_polys.insert(ei_other_ref_polys.end(), g[target(*eio, g)].referenced_polys.begin(), g[target(*eio, g)].referenced_polys.end());
						double I_with_other_neighbors = map ? mo.getIntersectionArea(g[q_prime].referenced_polys, ei_other_ref_polys) : mo.getIntersectionArea(ei_other_ref_polys, g[q_prime].referenced_polys);
						double A_without_other_neighbors = mo.getArea(!map, g[q_prime].referenced_polys) - I_with_other_neighbors;

						//COMPUTE IOU OF BIG MATCH INCLUDING ALL 4 POLYGONS
						std::vector<int> ref_polys_p, ref_polys_q;
						ref_polys_p.insert(ref_polys_p.end(), g[p].referenced_polys.begin(), g[p].referenced_polys.end());
						ref_polys_p.insert(ref_polys_p.end(), g[p_prime].referenced_polys.begin(), g[p_prime].referenced_polys.end());
						ref_polys_q.insert(ref_polys_q.end(), g[q].referenced_polys.begin(), g[q].referenced_polys.end());
						ref_polys_q.insert(ref_polys_q.end(), g[q_prime].referenced_polys.begin(), g[q_prime].referenced_polys.end());

						double IoU_all = !map ? mo.getIoU(ref_polys_p, ref_polys_q) : mo.getIoU(ref_polys_q, ref_polys_p);
						//END COMPUTING BIG MATCH IOU
						if (1 - lambda > IoU_sum - 2 * lambda || I <= A_without_p_prime || IoU_m_prime <= lambda || I_with_other_neighbors / A_without_other_neighbors >= lambda) {
							//if this happens, the bigger match would actually improve upon the singular matches and our property is hurt
							pq_is_optimal_match = false;
							continue;
						}
						else {
							//in this case, forming the two matches {p,q} and {p',q'} is actually better than one match could ever be, therefore we know that {p,q'} 
							//and {q',p} can never be in the same match and we can delete the edges from the graph
							cg->delete_edge(p, q_prime);
							cg->delete_edge(q, p_prime);

						}
					}
					else {
						pq_is_optimal_match = false;
						continue;
					}
				}
				if (!pq_is_optimal_match) continue;

				//if this point is reached, it is confirmed that for each adjacent polygon, there is a partner other than p/q, where forming two matches improves the target function
				//therefore p,q form a fix match and can be assigned in the solution right away
				if (!map) sol->addMatch(g[p].referenced_polys, g[q].referenced_polys, IoU_m - lambda);
				else  sol->addMatch(g[q].referenced_polys, g[p].referenced_polys, IoU_m - lambda);

				//delete the vertex with the highest index first so that the lower index is still the same after the deletion
				cg->delete_vertex(std::max(p, q));
				cg->delete_vertex(std::min(p, q));

				//update amount of vertices
				num_vertices = cg->get_Vi(map).size();

				//set p taking deletions into account
				i--;


			}
		}
	} while (graph_size_before_loop != cg->num_vertices());
	
}

bool pregroupIncludedPolygons(CandidateGraph* cg, const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2, const Localization& rtree1, const Localization& rtree2, const MapOverlay& mo, double lambda) {
	//indicator to remember if graph was modified
	bool function_modified_g = false;
	
	

	//iterate over every vertex
	for (int v = 0; v < cg->num_vertices(); v++) {
		//reset pointer
		Graph g = *cg->get_graph();

		//get included poylgons
		bool map = g[v].referenced_map;

		std::vector<int> inc = cg->getIncludedVertices(v);


		//pregrouping can only be performed if there are at least two included vertices
		std::vector<int> group;
		if (inc.size() > 1) {
			const std::vector<Polygon_wh>& polys = map ? polys1 : polys2;
			const Localization& rtree_opp = !map ? rtree1 : rtree2;
			//every vertex has to be checked for it to be inserted in to the new group
			for (const auto& v_inc : inc) {
				std::vector<int> neighbors;
				for (const auto& p_id : g[v_inc].referenced_polys) {
					rtree_opp.get_neighbors(polys[p_id], &neighbors);
				}
				//make neighbors vector unique and delete v from it 
				std::sort(neighbors.begin(), neighbors.end());
				neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
				for(const auto& rp : g[v].referenced_polys) neighbors.erase(std::remove(neighbors.begin(), neighbors.end(), rp), neighbors.end());

				//compute inside and outside area
				double inside_area = 0.0, outside_area = 0.0;
				if (!map) inside_area = mo.getIntersectionArea(neighbors, g[v_inc].referenced_polys);
				else inside_area = mo.getIntersectionArea(g[v_inc].referenced_polys, neighbors);
				outside_area = mo.getArea(!map, g[v_inc].referenced_polys) - inside_area;

				//if relation is less than lambda, neighbor can already be grouped as by Lemma it will always be matched to q
				if (outside_area == 0.0 || inside_area / outside_area < lambda) {
					//v_inc qualifies for group
					group.push_back(v_inc);
				}
			}

		}
		//a group should be formed if more than one vertex qualifies for group building
		if (group.size() > 1) {
			std::sort(group.begin(), group.end());

			std::vector<int> ref_polys;
			for (const auto& v_g : group) ref_polys.insert(ref_polys.end(), g[v_g].referenced_polys.begin(), g[v_g].referenced_polys.end());

			//delete all vertices from the group in decreasing order to have no index shifts on deletion
			for (int v_g = group.size() - 1; v_g >= 0; v_g--) {
				//reset pointer
				g = *cg->get_graph();
				cg->delete_vertex(group[v_g]);

				//reset current vertex index
				if (v_g < v) v--;
			}
			//insert new vertex representing the group
			cg->add_vertex(!map, ref_polys);

			function_modified_g = true;

		}

	}

	return function_modified_g;

}

void setInclusionProperties(std::vector<CandidateGraph> &g_vec, const Localization& rtree1, const Localization& rtree2, const MapOverlay& mo, double lambda, const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2) {
    //for all graphs
    for (auto& cg : g_vec) {
        //iterate over all pairs of vertices to check, if inclusion property occurs
        Graph::vertex_iterator vL, vLend, vR, vRend;
        //remember all substituting edges in order to be able to add them if the vertices did not get connected via a path beforehand
        std::vector<std::tuple<int, int, bool>> sub_edges;

        //for referenced polygons, an arbitrary layer of the multilayer graph can be accessed, take layer 0
        Graph* g = cg.get_graph();
        for (boost::tie(vL, vLend) = vertices(*g); vL != vLend; ++vL) {
            for (boost::tie(vR, vRend) = vertices(*g); vR != vRend; ++vR) {
                //do not compute edges twice
                if (*vL > *vR) continue;

                //vertices are of different maps, check for inclusion
                if ((*g)[*vL].referenced_map != (*g)[*vR].referenced_map) {
                    //pair of vertices where vL and vR are in different parts of the bipartite graph
                    std::vector<Polygon_wh> Lpolys = !(*g)[*vL].referenced_map ? polys1 : polys2;
                    std::vector<Polygon_wh> Rpolys = !(*g)[*vR].referenced_map ? polys1 : polys2;
                    Localization rtree = !(*g)[*vL].referenced_map ? rtree2 : rtree1;

                    int vertex_index1 = !(*g)[*vL].referenced_map ? *vL : *vR;
                    int vertex_index2 = !(*g)[*vL].referenced_map ? *vR : *vL;

                    //check if edge should be inserted at all, which is only the case if the sets of polygons overlap
                    if (mo.doOverlap((*g)[vertex_index1].referenced_polys, (*g)[vertex_index2].referenced_polys, 0.0)) {
                        double inter_area = mo.getIntersectionArea((*g)[vertex_index1].referenced_polys, (*g)[vertex_index2].referenced_polys);
                        if (inter_area > mo.getArea(0, (*g)[vertex_index1].referenced_polys) - inter_area) {
                            cg.addIncludedVertex(vertex_index2, vertex_index1);
                        }
                        if (inter_area > mo.getArea(1, (*g)[vertex_index2].referenced_polys) - inter_area) {
                            cg.addIncludedVertex(vertex_index1, vertex_index2);
                        }
                    }
                }
            }
        }
    }


}