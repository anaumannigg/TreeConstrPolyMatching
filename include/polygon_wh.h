#ifndef POLYGON_WH_H
#define POLYGON_WH_H

#include <boost/geometry/index/detail/rtree/visitors/iterator.hpp>

#include "cgal_includes.h"

//create custom class for handling polygons with holes
class Polygon_wh {
public:
    Polygon outer_boundary_polygon;
    std::vector<Polygon> hole_polygons = std::vector<Polygon>();

    //every polygon knows its ID in the global set
    int global_id = -1;

    //base constructor
    Polygon_wh(): outer_boundary_polygon(Polygon()) {};

    //constructor from an outer boundary
    Polygon_wh(const Polygon& _outer_boundary_polygon): outer_boundary_polygon(_outer_boundary_polygon) {};

    //constructor with holes
    Polygon_wh(const Polygon& _outer_boundary_polygon, std::vector<Polygon>::iterator _holes_begin, std::vector<Polygon>::iterator _holes_end):
        outer_boundary_polygon(_outer_boundary_polygon), hole_polygons(_holes_begin,_holes_end) {};

    // Copy constructor
    Polygon_wh(const Polygon_wh& other)
        : outer_boundary_polygon(other.outer_boundary_polygon),
          hole_polygons(other.hole_polygons),
          global_id(other.global_id) {};


    //return outer boundary polygon
    Polygon& outer_boundary();
    const Polygon& outer_boundary() const;

    //return holes
    const std::vector<Polygon>& holes() const;
    std::vector<Polygon>& holes();

    //return number of holes
    const int number_of_holes() const;

    //return hole iterators
    std::vector<Polygon>::const_iterator holes_begin() const;
    std::vector<Polygon>::const_iterator holes_end() const;

    //add a hole
    void add_hole(const Polygon& hole);

    //check if a point is on the bounded side of the polygon
    bool has_on_bounded_side(Point p);

    //return if this polygon does intersect the given polygon_wh p
    bool does_intersect(const Polygon_wh& p) const;

    //return if this polygon does intersect the given polygon p
    //invokes does_intersect(Polygon& p)
    bool does_intersect(Polygon& p);

    //return if the given polygon p lies entirely within a hole of the polygon
    bool lies_entirely_within_hole(const Polygon& p) const;

    //compute intersection area of polygon and given polygon p
    double intersection_area(const Polygon_wh& p);

    //compute hausdorff distance between two (sets of) polygons with holwws
    static double hausdorff_distance(const std::vector<const Polygon_wh*>& polys1, const std::vector<const Polygon_wh*>& polys2);

    //provide bounding box
    const CGAL::Bbox_2 bbox() const;

    //check if the polygon (and all of its holes) are valid
    bool is_valid();

    //uses CGAL's repair functionality to repair outer boundary and holes
    void repair();



    //reset memory
    void clear();


};


#endif //POLYGON_WH_H
