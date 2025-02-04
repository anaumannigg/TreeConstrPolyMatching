#ifndef _cgal_includes_included_
#define _cgal_includes_included_

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/draw_polygon_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/draw_polygon_set_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>

//for graph computations
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/connected_components.hpp>

//for convenience
#include <sstream>
#include <string>

//for convex hull
#include <CGAL/convex_hull_2.h>

//for centroid
#include <CGAL/centroid.h>

//for stability
#include <CGAL/Polygon_repair/repair.h>

//for computation
#include <math.h>

//for arrangements
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_non_caching_segment_traits_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_batched_point_location.h>

//for collision detection
#include <CGAL/Quadtree.h>

//for intersection points
#include <CGAL/Surface_sweep_2_algorithms.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/aff_transformation_tags.h>

//to write into files/output
#include <fstream>
#include <iostream>
#include <filesystem>

//for timing measurements
#include <chrono>

//for performance
#include <thread>

#define PI  3.14159265358979323846 

typedef CGAL::Exact_predicates_exact_constructions_kernel K;


typedef K::Point_2                              Point;
typedef K::Segment_2                            Segment;
typedef K::Line_2                               Line;
typedef K::Vector_2                             Vector;
typedef CGAL::Polygon_with_holes_2<K>           Polygon_wh;
typedef CGAL::Polygon_2<K>                      Polygon;
//typedef CGAL::Polygon_with_holes_2<K>			Polygon_wh;
//typedef CGAL::General_polygon_with_holes_2<K>	General_Polygon_wh_2;
typedef CGAL::Polygon_set_2<K>                  Polygon_set_2;

//for arrangement
typedef CGAL::Arr_segment_traits_2<K> ArrTraits;
typedef CGAL::Arr_face_extended_dcel<ArrTraits,int> Face_labeled_DCEL;
typedef ArrTraits::Segment_2 ArrSegment;
typedef CGAL::Arrangement_2<ArrTraits,Face_labeled_DCEL > Arrangement;

//for graph computations
using namespace boost;
struct VertexProps {
	//referenced map, =0 for OSM, =1 for ATKIS
	bool referenced_map;
	std::vector<int> referenced_polys;
    int id = -1;
};
struct EdgeProps {
    int id = -1;
    double weight;
};
//define undirected graph as default
typedef adjacency_list<vecS, vecS, undirectedS, VertexProps, EdgeProps> Graph;
typedef std::pair<int, int> Edge;
typedef typename graph_traits<Graph>::vertex_descriptor Vertex;

//also define directed Graph
typedef adjacency_list<vecS, vecS, directedS, VertexProps, EdgeProps> DiGraph;
typedef typename graph_traits<DiGraph>::vertex_descriptor DiVertex;


using std::cout; using std::endl;

//standard vector operations
template <typename T>
auto normalize(T const& V)
{
    auto const slen = V.squared_length();
    auto const d = CGAL::approximate_sqrt(slen);
    return V / d;
}

template <typename T>
auto rotate(T const& V, double angle) {
    return T(cos(angle) * V.x()- sin(angle)*V.y(), sin(angle)*V.x() + cos(angle)*(V.y()));
}

//computes angle between two vectors, input in COUNTERCLOCKWISE order
template <typename T>
auto get_angle(T const& w, T const& v) {
    //double angle = atan2(to_double(v.y() * w.x() - v.x() * w.y()), to_double(v.x() * w.x() + v.y() * w.y()));
    double angle = atan2(v.y() * w.x() - v.x() * w.y(), v.x() * w.x() + v.y() * w.y());
    if (angle < 0) angle += 2 * PI;
    return angle;
}




#endif