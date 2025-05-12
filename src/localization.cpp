#include "../include/localization.h"

//LOCALIZATION CLASS

Localization::Localization(std::vector<Polygon_wh>& polys) : polys(polys) {
    //init r-tree
    this->rtree = bgi::rtree< Value_rtree, bgi::quadratic<16>>();

    for (int i = 0; i < polys.size(); i++) {

        //input the bounding boxes of the polygons, connected with their index, into the r-tree
        auto bbox = polys[i].bbox();
        //get lowerleft and upperright point of the bounding box
        Point_rtree ll = Point_rtree(bbox.xmin(), bbox.ymin());
        Point_rtree ur = Point_rtree(bbox.xmax(), bbox.ymax());
        //create box for rtree
        Box_rtree b(ll, ur);
        //insert polygon and original index into r-tree
        rtree.insert(std::make_pair(b, i));


    }
}

//locates the input polygon p within the r-tree of the instance and writes neighbors into given vector
//returns true if p has intersecting neighbors within the r tree, else false
bool Localization::get_neighbors(const Polygon_wh& p, std::vector<Polygon_wh>* neighbors) const {
    //compute bounding box of circular arc
    auto poly_bbox = p.bbox();
    Point_rtree ll = Point_rtree(poly_bbox.xmin(), poly_bbox.ymin());
    Point_rtree ur = Point_rtree(poly_bbox.xmax(), poly_bbox.ymax());
    Box_rtree query_box = Box_rtree(ll, ur);

    //query r-tree
    std::vector<Value_rtree> query_result = this->query(query_box);

    for (const auto& b : query_result) {
        //remember every polygon 
        Polygon_wh query_poly = this->polys[b.second];
        if (query_poly.does_intersect(p)) {
            neighbors->push_back(query_poly);
        }
    }

    if (neighbors->size() > 0) return true;
    else return false;
}

//locates the input polygon p within the r-tree of the instance and writes neighbor indices into given vector
//returns true if p has intersecting neighbors within the r tree, else false
bool Localization::get_neighbors(const Polygon_wh& p, std::vector<int>* neighbors) const {
    //compute bounding box of input polygon
    auto poly_bbox = p.bbox();
    Point_rtree ll = Point_rtree(poly_bbox.xmin(), poly_bbox.ymin());
    Point_rtree ur = Point_rtree(poly_bbox.xmax(), poly_bbox.ymax());
    Box_rtree query_box = Box_rtree(ll, ur);

    //query r-tree
    std::vector<Value_rtree> query_result = this->query(query_box);


    auto p_boundary = p.outer_boundary();

    for (const auto& b : query_result) {
        //remember every polygon
        Polygon_wh query_poly = this->polys[b.second];

        //fix query poly if necessary
        if (p.does_intersect(query_poly)) {
            neighbors->push_back(b.second);
        }
    }
    //sort output vector
    std::sort(neighbors->begin(), neighbors->end());
    //make output vector unique
    neighbors->erase(std::unique(neighbors->begin(), neighbors->end()), neighbors->end());

    if (neighbors->size() > 0) return true;
    else return false;
}

//overloads the get neighbors function for simple polygons, handles them as outer boundaries
bool Localization::get_neighbors(const Polygon& p, std::vector<int>* neighbors) const {
    return get_neighbors(Polygon_wh(p), neighbors);
}

//locates the input point and writes intersected polygons into neighbors
bool Localization::get_neighbors(const Point& p, std::vector<int>* neighbors, bool consider_holes) const {
    auto poly_bbox = p.bbox();
    //buffer bbox with some epsilon for stability
    double buffer = 1e-2;
    Point_rtree ll = Point_rtree(poly_bbox.xmin()-buffer, poly_bbox.ymin()-buffer);
    Point_rtree ur = Point_rtree(poly_bbox.xmax()+buffer, poly_bbox.ymax()+buffer);
    Box_rtree query_box = Box_rtree(ll, ur);

    //query r-tree
    std::vector<Value_rtree> query_result = this->query(query_box);

    for (const auto& b : query_result) {
        //remember every polygon 
        Polygon_wh query_poly = this->polys[b.second];
        if (!query_poly.outer_boundary().has_on_unbounded_side(p)) {
            //point lies inside the outer boundary, make sure the point does not lie inside a hole
            bool p_is_in_hole = false;

            //only consider holes if consider_holes = true (else only the outer boundary is considered)
            if(consider_holes) {
                for (const auto &hole: query_poly.holes()) {
                    if (!hole.has_on_unbounded_side(p)) {
                        p_is_in_hole = true;
                        break;
                    }
                }
            }


            if(!p_is_in_hole) neighbors->push_back(b.second);
        }
    }
    //sort output vector
    std::sort(neighbors->begin(), neighbors->end());
    //make output vector unique
    neighbors->erase(std::unique(neighbors->begin(), neighbors->end()), neighbors->end());
    if (neighbors->size() > 0) return true;
    else return false;
}

//locates the input segment and writes intersected polygons into neighbors
//returns true if p has intersecting neighbors, else false
bool Localization::get_neighbors(const Segment& p, std::vector<int>* neighbors) const {
    //if the segment represents a single point, return
    if (p.source() == p.target()) return this->get_neighbors(p.source(),neighbors);

    //build poly from segment
    Polygon p_poly;
    p_poly.push_back(p.source()); p_poly.push_back(p.target());

    //compute bounding box of segment
    auto poly_bbox = p.bbox();
    Point_rtree ll = Point_rtree(poly_bbox.xmin(), poly_bbox.ymin());
    Point_rtree ur = Point_rtree(poly_bbox.xmax(), poly_bbox.ymax());
    Box_rtree query_box = Box_rtree(ll, ur);

    //query r-tree
    std::vector<Value_rtree> query_result = this->query(query_box);

    for (const auto& b : query_result) {
        //retreive query polygon
        Polygon_wh query_poly = this->polys[b.second];

        //make sure the input segment is not one segment on the boundary of the polygon, as then the intersection check will crash
        bool is_boundary_edge = false;
        for(int part = 0; part < 1 + query_poly.number_of_holes(); part++) {
            Polygon qp = query_poly.outer_boundary();
            if(part > 0) {
                auto hole_it = query_poly.holes_begin();
                hole_it = std::next(hole_it, part-1);
                qp = *hole_it;
            }

            for (const auto &e: qp.edges()) {
                if (e.source() == p.source() && e.target() == p.target() ||
                    e.target() == p.source() && e.source() == p.target()) {
                    is_boundary_edge = true;
                    break;
                }
            }
        }

        if (is_boundary_edge) {
            //collision found
            neighbors->push_back(b.second);
        }
        //input segment is not an edge of the polygon, check for intersection
        else {
            if (!query_poly.does_intersect(p_poly)) {
                neighbors->push_back(b.second);
            }
        }
    }
    //sort output vector
    std::sort(neighbors->begin(), neighbors->end());
    //make output vector unique
    neighbors->erase(std::unique(neighbors->begin(), neighbors->end()), neighbors->end());
    if (neighbors->size() > 0) return true;
    else return false;
}

//neighbors function with constraint: an intersecting polygon is only returned as a neighbor, iff union / smaller_poly_area is >= the argument threshold
bool Localization::get_neighbors(const Polygon_wh& p, std::vector<int>* neighbors, double neighboring_threshold) const {
    //compute bounding box of circular arc
    auto poly_bbox = p.bbox();
    Point_rtree ll = Point_rtree(poly_bbox.xmin(), poly_bbox.ymin());
    Point_rtree ur = Point_rtree(poly_bbox.xmax(), poly_bbox.ymax());
    Box_rtree query_box = Box_rtree(ll, ur);

    //query r-tree
    std::vector<Value_rtree> query_result = this->query(query_box);

    for (const auto& b : query_result) {
        //remember every polygon 
        Polygon_wh query_poly = this->polys[b.second];

        
        if (query_poly.does_intersect(p)) {
            //polys do intersect,but further requirements have to be met
            //compute intersection area
            double inter_area = query_poly.intersection_area(p);

            //get smaller polygon area
            double smaller_poly_area = std::min(to_double(query_poly.outer_boundary().area()),
                                                to_double(p.outer_boundary().area()));

            if(inter_area/smaller_poly_area > neighboring_threshold) neighbors->push_back(b.second);
        }
    }
    //sort output vector
    std::sort(neighbors->begin(),neighbors->end());
    if (neighbors->size() > 0) return true;
    else return false;
}

//neighbors function with constraint: an intersecting polygon is only returned as a neighbor, iff union / smaller_poly_area is >= the argument threshold, speedup via arrangement
bool Localization::get_neighbors(const Polygon_wh& p, int p_index, bool p_map, std::vector<int>* neighbors, double neighboring_threshold, MapOverlay mo) const {
    //compute bounding box of circular arc
    auto poly_bbox = p.bbox();
    Point_rtree ll = Point_rtree(poly_bbox.xmin(), poly_bbox.ymin());
    Point_rtree ur = Point_rtree(poly_bbox.xmax(), poly_bbox.ymax());
    Box_rtree query_box = Box_rtree(ll, ur);

    //query r-tree
    std::vector<Value_rtree> query_result = this->query(query_box);


    

    for (const auto& b : query_result) {
        //remember every polygon 
        Polygon_wh query_poly = this->polys[b.second];

        if (p.does_intersect(query_poly)) {

            std::vector<int> indices1 = !p_map ? std::vector<int> {p_index} : std::vector<int> {(int)b.second};
            std::vector<int> indices2 = !p_map ? std::vector<int> {(int)b.second} : std::vector<int>{ p_index };

            double inter_area = mo.getIoU(indices1, indices2);

            if (inter_area > 0) {

                //get smaller polygon area
                double smaller_poly_area = std::min(to_double(query_poly.outer_boundary().area()),
                                                    to_double(p.outer_boundary().area()));

                if (inter_area / smaller_poly_area > neighboring_threshold) neighbors->push_back(b.second);
            }
        }
    }
    //sort output vector
    std::sort(neighbors->begin(), neighbors->end());
    if (neighbors->size() > 0) return true;
    else return false;
}

bool Localization::get_neighbors_fully_included(const Polygon_wh& p, std::vector<int>* neighbors) const {
    //compute bounding box of polygon
    auto poly_bbox = p.bbox();
    Point_rtree ll = Point_rtree(poly_bbox.xmin(), poly_bbox.ymin());
    Point_rtree ur = Point_rtree(poly_bbox.xmax(), poly_bbox.ymax());
    Box_rtree query_box = Box_rtree(ll, ur);

    //query r-tree
    std::vector<Value_rtree> query_result = this->query(query_box);
    for (const auto& b : query_result) {
        //remember every polygon, but only if it has not been considered yet and the query poly lies entirely within it
        Polygon_wh query_poly = this->polys[b.second];
        if (std::find(neighbors->begin(),neighbors->end(),b.second) == neighbors->end()) {
            //polygon has not been considered yet as neighbor
            //check if every point of p lies inside the polygon (more stable and efficient than intersection area check)
            bool all_points_lie_inside = true;
            for (const auto& point : p.outer_boundary().vertices()) {
                //also consider holes
                for(int part=0; part < 1 + query_poly.number_of_holes(); part++) {
                    Polygon query_poly_part = query_poly.outer_boundary();
                    if(part > 0) {
                        auto hole_it = query_poly.holes_begin();
                        hole_it = std::next(hole_it,part-1);
                        query_poly_part = *hole_it;
                    }
                    //outer boundary case
                    if (part == 0 && query_poly_part.has_on_unbounded_side(point)) {
                        //point lies outside of the polygon
                        all_points_lie_inside = false;
                        break;
                    }
                    //hole case
                    else if(part > 0 && !query_poly_part.has_on_unbounded_side(point)) {
                        //point lies in hole
                        all_points_lie_inside = false;
                        break;
                    }
                }
            }

            if (all_points_lie_inside) neighbors->push_back(b.second);

        }
    }
    //sort output vector
    std::sort(neighbors->begin(), neighbors->end());
    //make output vector unique
    neighbors->erase(std::unique(neighbors->begin(), neighbors->end()), neighbors->end());
    if (neighbors->size() > 0) return true;
    else return false;
}

bool Localization::get_neighbors_majorly_included(const Polygon_wh& p, std::vector<int>* neighbors, double inclusion_threshold) const {
    //compute bounding box of circular arc
    auto poly_bbox = p.bbox();
    Point_rtree ll = Point_rtree(poly_bbox.xmin(), poly_bbox.ymin());
    Point_rtree ur = Point_rtree(poly_bbox.xmax(), poly_bbox.ymax());
    Box_rtree query_box = Box_rtree(ll, ur);

    //query r-tree
    std::vector<Value_rtree> query_result = this->query(query_box);

    for (const auto& b : query_result) {
        //remember every polygon, but only if it has not ben considered yet and the query poly lies entirely within it
        Polygon_wh query_poly = this->polys[b.second];

        if (std::find(neighbors->begin(), neighbors->end(), b.second) == neighbors->end()) {
            //polygon has not been considered yet as neighbor
            
            //compute o_tilde for p and b
            double inter_area = query_poly.intersection_area(p);

            //get smaller polygon area
            double smaller_poly_area = std::min(to_double(query_poly.outer_boundary().area()),
                                                to_double(p.outer_boundary().area()));

            if (inter_area / smaller_poly_area > inclusion_threshold) neighbors->push_back(b.second);
     
        }
    }
    //sort output vector
    std::sort(neighbors->begin(), neighbors->end());
    //make output vector unique
    neighbors->erase(std::unique(neighbors->begin(), neighbors->end()), neighbors->end());
    if (neighbors->size() > 0) return true;
    else return false;
}

bool Localization::get_neighbors_with_minimum_intersection_proportion(const Polygon_wh& p, std::vector<int>* neighbors, double inclusion_threshold) const {
    //compute bounding box of circular arc
    auto poly_bbox = p.bbox();
    Point_rtree ll = Point_rtree(poly_bbox.xmin(), poly_bbox.ymin());
    Point_rtree ur = Point_rtree(poly_bbox.xmax(), poly_bbox.ymax());
    Box_rtree query_box = Box_rtree(ll, ur);

    //query r-tree
    std::vector<Value_rtree> query_result = this->query(query_box);

    for (const auto& b : query_result) {
        //remember every polygon, but only if it has not ben considered yet and the query poly lies entirely within it
        Polygon_wh query_poly = this->polys[b.second];

        if (std::find(neighbors->begin(), neighbors->end(), b.second) == neighbors->end()) {
            //polygon has not been considered yet as neighbor

            //compute o_tilde for p and b
            double inter_area = query_poly.intersection_area(p);
            //check if at last the necessary threshold of the query poly area is in the intersection
            if (inter_area / to_double(query_poly.outer_boundary().area()) > inclusion_threshold) neighbors->push_back(b.second);

        }
    }
    //sort output vector
    std::sort(neighbors->begin(), neighbors->end());
    //make output vector unique
    neighbors->erase(std::unique(neighbors->begin(), neighbors->end()), neighbors->end());
    if (neighbors->size() > 0) return true;
    else return false;
}

std::vector<Value_rtree> Localization::query(const Box_rtree& query_box) const {
    std::vector<Value_rtree> query_result;
    this->rtree.query(bgi::intersects(query_box),std::back_inserter(query_result));
    return query_result;
}

bool Localization::are_adjacent(const Polygon_wh& polyA, int polyB_index) const {
    //get neighbors of poly A
    std::vector<int> neighbors;
    //note: here, different get_neighbors modes could be applied in order to adapt overlap checks
    this->get_neighbors(polyA, &neighbors);

    //check, if polyB is in the neighbors list of A
    if (std::find(neighbors.begin(), neighbors.end(), polyB_index) != neighbors.end()) return true;
    else return false;
}

//END LOCALIZATION CLASS

// MAP OVERLAY CLASS

MapOverlay::MapOverlay(std::vector<Polygon_wh> polys1, std::vector<Polygon_wh> polys2) {
    this->arr = Arrangement();

    //insert all segments of polygons into the arrangement
    std::vector<ArrSegment> segments;
    for (const auto& polys : { polys1,polys2 }) {
        for (const auto& poly : polys) {
            //add outer boundary edges
            for (const auto& e : poly.outer_boundary().edges()) {
                segments.push_back(e);
            }

            //add edges in holes
            for(const auto hole : poly.holes()) {
                for(const auto& e : hole.edges()) {
                    segments.push_back(e);
                }
            }
        }
    }
    CGAL::insert(arr, segments.begin(), segments.end());

    this->areas = std::vector<double>();
    this->areas.resize(arr.number_of_faces(), 0.0);
    this->poly1_covered_faces = std::vector<std::vector<int>>(polys1.size());
    this->poly2_covered_faces = std::vector<std::vector<int>>(polys2.size());


}

void MapOverlay::assignFaces(std::vector<Polygon_wh> polys1, std::vector<Polygon_wh> polys2, Localization rtree1, Localization rtree2) {
    int face_id = 0;

    for (auto fit = this->arr.faces_begin(); fit != this->arr.faces_end(); fit++) {


        //if the face is the outer face or has no outer ccb, do not compute the face area
        if (fit->is_unbounded() || !fit->has_outer_ccb()) {
            face_id++;
            continue;
        }

        //compute face area
        double face_area = 0.0;

        //first compute and add outer ccb
        Polygon outer_ccb;
        Arrangement::Ccb_halfedge_const_circulator e = fit->outer_ccb();
        auto loop = e;
        do {
            auto x = e->source();
            outer_ccb.push_back(e->source()->point());
        } while (++e != loop);
        //make sure that the area can actually be computed, there might be artifacts leading to a self-intersecting outer ccb
        if (outer_ccb.is_simple()) {
            face_area += to_double(outer_ccb.area());

            //write area into vector
            this->areas[face_id] = face_area;


            //check, which polygons overlap the current face
            std::vector<int> neighbors1, neighbors2;
            rtree1.get_neighbors(outer_ccb, &neighbors1);
            rtree2.get_neighbors(outer_ccb, &neighbors2);

            //remember for those polygons, that the facet is included within them 
            for (const auto& n1 : neighbors1) this->poly1_covered_faces[n1].push_back(face_id);
            for (const auto& n2 : neighbors2) this->poly2_covered_faces[n2].push_back(face_id);


            face_id++;
        }
    }

}

std::vector<int> MapOverlay::getCoveredFaces(bool map, int poly_index) const {
    std::vector<std::vector<int>> polys = !map ? this->poly1_covered_faces : this->poly2_covered_faces;
    return polys[poly_index];
}

std::vector<int> MapOverlay::getCoveredFaces(bool map, std::vector<int> poly_indices) const {
    std::vector<int> covered_faces;
    for (const auto& poly : poly_indices) {
        std::vector<int> cov_faces_of_poly = getCoveredFaces(map, poly);
        if (cov_faces_of_poly.size() > 0) covered_faces.insert(covered_faces.end(), cov_faces_of_poly.begin(), cov_faces_of_poly.end());
    }
    return covered_faces;
}

double MapOverlay::getIoU(std::vector<int> indices1, std::vector<int> indices2) const {
    //get common intersected faces of polygons of both sets
    std::vector<int> faces1, faces2;
    for (const auto& i1 : indices1) {
        std::vector<int> covered_faces = this->getCoveredFaces(0, i1);
        faces1.insert(faces1.end(), covered_faces.begin(), covered_faces.end());
    }
    for (const auto& i2 : indices2) {
        std::vector<int> covered_faces = this->getCoveredFaces(1, i2);
        faces2.insert(faces2.end(), covered_faces.begin(), covered_faces.end());
    }

    return this->getIoUviaCoveredFaces(faces1, faces2);
}

double MapOverlay::getIntersectionArea(std::vector<int> indices1, std::vector<int> indices2) const {
    //get common intersected faces of polygons of both sets
    std::vector<int> faces1, faces2;
    for (const auto& i1 : indices1) {
        std::vector<int> covered_faces = this->getCoveredFaces(0, i1);
        faces1.insert(faces1.end(), covered_faces.begin(), covered_faces.end());
    }
    for (const auto& i2 : indices2) {
        std::vector<int> covered_faces = this->getCoveredFaces(1, i2);
        faces2.insert(faces2.end(), covered_faces.begin(), covered_faces.end());
    }

    return this->getIntersectionviaCoveredFaces(faces1, faces2);
}

double MapOverlay::getIntersectionviaCoveredFaces(std::vector<int> faces1, std::vector<int> faces2) const {
    //sort both lists and make unique
    std::sort(faces1.begin(), faces1.end());
    faces1.erase(unique(faces1.begin(), faces1.end()), faces1.end());
    std::sort(faces2.begin(), faces2.end());
    faces2.erase(unique(faces2.begin(), faces2.end()), faces2.end());

    //get intersection and union of face indices of both sets
    std::vector<int> inter_faces;
    std::set_intersection(faces1.begin(), faces1.end(), faces2.begin(), faces2.end(), std::back_inserter(inter_faces));


    double inter_area = 0.0;
    for (const auto& inter_face : inter_faces) inter_area += this->areas[inter_face];
    return inter_area;
}

double MapOverlay::getUnionviaCoveredFaces(std::vector<int> faces1, std::vector<int> faces2) const {
    //sort both lists and make unique
    std::sort(faces1.begin(), faces1.end());
    faces1.erase(unique(faces1.begin(), faces1.end()), faces1.end());
    std::sort(faces2.begin(), faces2.end());
    faces2.erase(unique(faces2.begin(), faces2.end()), faces2.end());

    //get intersection and union of face indices of both sets
    std::vector<int> union_faces;
    union_faces.insert(union_faces.end(), faces1.begin(), faces1.end());
    union_faces.insert(union_faces.end(), faces2.begin(), faces2.end());
    sort(union_faces.begin(), union_faces.end());
    union_faces.erase(unique(union_faces.begin(), union_faces.end()), union_faces.end());


    double union_area = 0.0;
    for (const auto& union_face : union_faces) union_area += this->areas[union_face];
    return union_area;
}

double MapOverlay::getIoUviaCoveredFaces(std::vector<int> faces1, std::vector<int> faces2) const {
    //sort both lists and make unique
    std::sort(faces1.begin(), faces1.end());
    faces1.erase(unique(faces1.begin(), faces1.end()), faces1.end());
    std::sort(faces2.begin(), faces2.end());
    faces2.erase(unique(faces2.begin(), faces2.end()), faces2.end());

    //get intersection and union of face indices of both sets
    std::vector<int> inter_faces, union_faces;
    std::set_intersection(faces1.begin(), faces1.end(), faces2.begin(), faces2.end(), std::back_inserter(inter_faces));
    union_faces.insert(union_faces.end(), faces1.begin(), faces1.end());
    union_faces.insert(union_faces.end(), faces2.begin(), faces2.end());
    sort(union_faces.begin(), union_faces.end());
    union_faces.erase(unique(union_faces.begin(), union_faces.end()), union_faces.end());


    double inter_area = 0.0, union_area = 0.0;
    for (const auto& inter_face : inter_faces) inter_area += this->areas[inter_face];

    for (const auto& union_face : union_faces) union_area += this->areas[union_face];
    return inter_area / union_area;
}

double MapOverlay::getArea(bool map, std::vector<int> indices) const {
    std::vector<int> faces;
    for (const auto& i : indices) {
        std::vector<int> covered_faces = this->getCoveredFaces(map, i);
        faces.insert(faces.end(), covered_faces.begin(), covered_faces.end());
    }

    return this->getAreaofCoveredFaces(faces);

}

double MapOverlay::getAreaofCoveredFaces(std::vector<int> covered_faces) const {
    double area = 0.0;
    for (const auto& f : covered_faces) area += this->areas[f];
    return area;
}

bool MapOverlay::doOverlap(std::vector<int> indices1, std::vector<int> indices2, double epsilon) const {
    //get common intersected faces of polygons of both sets
    std::vector<int> faces1, faces2;

    for (const auto& i1 : indices1) {
        std::vector<int> covered_faces = this->getCoveredFaces(0, i1);
        faces1.insert(faces1.end(), covered_faces.begin(), covered_faces.end());
    }
    for (const auto& i2 : indices2) {
        std::vector<int> covered_faces = this->getCoveredFaces(1, i2);
        faces2.insert(faces2.end(), covered_faces.begin(), covered_faces.end());
    }


    return this->doOverlapviaCoveredFaces(faces1, faces2, epsilon);

}

bool MapOverlay::doOverlapviaCoveredFaces(std::vector<int> covered_faces1, std::vector<int> atkis_covered_faces, double epsilon) const {
    double inter_area = this->getIntersectionviaCoveredFaces(covered_faces1, atkis_covered_faces);

    double min_area = std::min(this->getAreaofCoveredFaces(covered_faces1), this->getAreaofCoveredFaces(atkis_covered_faces));

    return inter_area / min_area > epsilon;
}

// END MAP OVERLAY CLASS

//helper struct for the comparison of segments to keep a unique set
struct SegmentComparator {
    bool operator()(const Segment& s1, const Segment& s2) const {
        // Compare the x values of the source points first
        if (s1.source().x() != s2.source().x()) {
            return s1.source().x() < s2.source().x();
        }
        // If x values are equal, compare the y values of the source points
        if (s1.source().y() != s2.source().y()) {
            return s1.source().y() < s2.source().y();
        }
        // If the source points are the same, compare the target points
        if (s1.target().x() != s2.target().x()) {
            return s1.target().x() < s2.target().x();
        }
        return s1.target().y() < s2.target().y();
    }
};

// Helper function to get the minimum x-value of a polygon
double min_x_value(const Polygon_wh& polygon) {
    //note that only outer boundary vertices need to be considered for min x value computation
    // as holes per definition are inside the polygon
    double min_x = to_double(polygon.outer_boundary().vertices_begin()->x());
    for (auto it = polygon.outer_boundary().vertices_begin(); it != polygon.outer_boundary().vertices_end(); ++it) {
        if (it->x() < min_x) {
            min_x = to_double(it->x());
        }
    }



    return min_x;
}

// Function to compute the union of two vectors of polygons
std::vector<Polygon_wh> merge(const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2) {

    //if the sets of the polygons are too big, the arrangement might exceed the RAM-size
    //we choose a limit of a total of 10^5 polygons. If this is exceeded, we split the
    //sets into subsets, which we iteratively union to reduce the individual arrangement size
    int split_size = 1e5;
    //vectors to store interval limits of polys
    std::vector<int> intervals1 = {0},intervals2 = {0};
    //vectors to store ids of polys sorted w.r.t. their minimal x values
    std::vector<int> ids1(polys1.size()),ids2(polys2.size());
    //initialize id vectors
    std::iota(ids1.begin(),ids1.end(),0); std::iota(ids2.begin(),ids2.end(),0);

    //check if split needs to be performed
    if(polys1.size() + polys2.size() > 2*split_size ) {
        //total number of polygons exceeds threshold, split is needed

        //sort polygons w.r.t. x-axis to form meaningful subsets

        // Precompute min x-values and store indices
        std::vector<double> min_x_values1(polys1.size());
        for (std::size_t i = 0; i < polys1.size(); ++i) {
            min_x_values1[i] = min_x_value(polys1[i]);
        }


        // Sort indices based on the precomputed min x-values
        std::sort(ids1.begin(), ids1.end(),
                  [&min_x_values1](std::size_t i1, std::size_t i2) {
                      return min_x_values1[i1] < min_x_values1[i2];
                  });

        //same for polys2
        std::vector<double> min_x_values2(polys2.size());
        for (std::size_t i = 0; i < polys2.size(); ++i) {
            min_x_values2[i] = min_x_value(polys2[i]);
        }


        // Sort indices based on the precomputed min x-values
        std::sort(ids2.begin(), ids2.end(),
                  [&min_x_values2](std::size_t i1, std::size_t i2) {
                      return min_x_values2[i1] < min_x_values2[i2];
                  });


        //now compute the subsets
        //we want them to be reasonably small so we choose x splits depending on the bigger set
        auto& min_x_values = polys1.size()> polys2.size()? min_x_values1 : min_x_values2;
        auto& intervals = polys1.size()> polys2.size()? intervals1 : intervals2;
        auto& ids = polys1.size() > polys2.size()? ids1 : ids2;
        auto size = polys1.size() > polys2.size()? polys1.size() : polys2.size();


        int i=split_size;
        //limits to split the second set by
        //they are determined by every split_size polygons of the bigger set
        std::vector<double> x_limits;

        while(i < min_x_values.size()) {
            x_limits.push_back(min_x_values[ids[i]]);
            intervals.push_back(i);
            i+=split_size;
        }
        //if size % split_size != 0, define last smaller set
        if( i != min_x_values.size()-1) {
            intervals.push_back(size-1);
        }

        //now the smaller set needs to be split into intervals depending on the x limits
        auto& min_x_values_2 = polys1.size()> polys2.size()? min_x_values2 : min_x_values1;
        auto& intervals_2 = polys1.size() > polys2.size() ? intervals2:intervals1;
        auto& ids_2 = polys1.size() > polys2.size() ? ids2:ids1;
        auto size_2 = polys1.size() > polys2.size() ? polys2.size() : polys1.size();


        i=0;
        for(const auto& x_lim : x_limits) {
            while(i < size_2 && min_x_values_2[ids_2[i]] < x_lim) i++;
            if(i < size_2) {
                intervals_2.push_back(i);
            }
            else break;
        }
        if(i!= min_x_values_2.size()-1) intervals_2.push_back(size_2-1);

        //fill up with empty intervals to assure equal interval count to the bigger set
        auto& bigger_intervals = polys1.size() > polys2.size() ? intervals1 : intervals2;
        while(intervals.size() != bigger_intervals.size()) intervals.push_back(--i);




    }
    else {
        //split does not need to be performed, simply create one interval
        intervals1.push_back(polys1.size()-1);
        intervals2.push_back(polys2.size()-1);
    }

    assert(intervals1.size() == intervals2.size() && "intervals of split polygons should have same size!");

    //get x-values of interval borders
    std::vector<double> interval_borders;
    for(int i=1; i<intervals1.size()-1;i++) interval_borders.push_back(min_x_value(polys1[ids1[intervals1[i]]]));
    //set the last limit to max
    interval_borders.push_back(std::numeric_limits<double>::max());

    //remember already computed unions
    std::vector<Polygon_wh> merged_polys;

    //remember incomplete unions (those polygons, that intersect the x_lim, which might be incomplete and need to be reconsidered
    //in the next step)
    std::vector<Polygon> tmp_polys;

    //iterate over all intervals
    for(int interval=1; interval < intervals1.size(); interval++) {
        //create a set of edges of the corresponding intervals of both sets of polygons as well as already merged ones
        std::set<Segment, SegmentComparator> segments;
        //start with the already merged polys
        for(const auto& mp : tmp_polys){
            for(const auto& e : mp.edges()) {
                segments.insert(e);
            }
        }
        tmp_polys.clear();

        //segments of polys1
        for(int i1 = intervals1[interval-1]; i1 <= intervals1[interval]; i1++) {
            for(const auto& e : polys1[ids1[i1]].outer_boundary().edges()) {
                segments.insert(e);
            }
        }
        //segments of polys2
        for(int i2 = intervals2[interval-1]; i2 <= intervals2[interval]; i2++) {
            for(const auto& e : polys2[ids2[i2]].outer_boundary().edges()) {
                segments.insert(e);
            }
        }

        //create an arrangement using all segments
        Arrangement arr;
        CGAL::insert(arr, segments.begin(), segments.end());

        //label the faces of the arrangement with their ids
        int f_id = 0;
        for (auto f = arr.faces_begin(); f != arr.faces_end(); f++) {
            f->set_data(f_id++);
        }

        //create an r-tree with all polygons to check, if a face lies on a polygon
        std::vector<Polygon_wh> all_polys;
        for (auto &polys: {polys1, polys2}) all_polys.insert(all_polys.end(), polys.begin(), polys.end());
        Localization idx(all_polys);


        //collect all faces, that lie within at least one polygon of the two sets
        std::vector<Arrangement::Face_const_handle> poly_faces;
        std::vector<int> poly_face_ids;
        for (auto f = arr.faces_begin(); f != arr.faces_end(); f++) {
            //get centroid of the face
            if (f->has_outer_ccb()) {
                Polygon f_p;

                auto e = f->outer_ccb();
                auto e_loop = e;
                do {
                    f_p.push_back(e->source()->point());
                } while (++e != e_loop);

                //due to artifacts (specifically holes that share a single point with the outer boundary, that CGAL
                //cannot handle), we might need to repair the polygon
                f_p = PolygonOperations::fixIfNotSimple(f_p);

                //create a point sample inside the polygon for r-tree query
                //(from a vertex, turn inside the polygon and sample a point)
                Point sample;
                auto v_c = f_p.vertices_circulator();

                //remember if valid sample could be found
                //due to artifacts, there may be invalid faces, which should be ignored
                bool sample_is_valid = true;
                auto v_c_loop = v_c;



                do {
                    Point source = *v_c;
                    Point target = *(v_c + 1);
                    Vector dir(source, target);

                    auto dx = target.x() - source.x();
                    auto dy = target.y() - source.y();
                    //turn vector inside perpendicularly
                    Vector inside(-dy,dx);
                    double eps = 1e-5;
                    sample = source + 0.5 * dir + eps * normalize(inside);

                    v_c++;

                    //check if complete loop has been performed, if so, break
                    if(v_c == v_c_loop) {
                        sample_is_valid = false;
                        break;
                    }

                } while (!f_p.has_on_bounded_side(sample));

                if(sample_is_valid) {

                    //if the sampled point lies at least within one polygon, remember the face
                    std::vector<int> neighbors;
                    if (idx.get_neighbors(sample, &neighbors)) {
                        poly_faces.push_back(f);
                        poly_face_ids.push_back(f->data());
                    }
                }
                else std::cerr << "warning: invalid sample in arrangement traversal." << endl;


            }
        }

        //sort poly face id vector for quick check if a face is a poly face
        std::sort(poly_face_ids.begin(), poly_face_ids.end());

        //DEBUG:export arrangement
        //writeToShapeFile(arr,"../export/debug/arr",poly_face_ids);

        //initialize a vector to keep track of visited faces
        std::vector<bool> visited_faces;
        visited_faces.resize(arr.number_of_faces(), false);

        std::sort(poly_face_ids.begin(), poly_face_ids.end());

        //start out with a poly face and trace the boundary as long as finding other selected faces and not arriving at the starting edge again
        for (const auto &p_face: poly_faces) {
            if (!visited_faces[p_face->data()] && p_face->number_of_outer_ccbs() == 1) {

                //mark face as selected
                visited_faces[p_face->data()] = true;

                //check if a selected face, which lies on the boundary of a connected subset is found
                bool is_on_boundary = false;
                Arrangement::Halfedge_const_handle curr;
                auto e_it = p_face->outer_ccb();
                auto e_loop = e_it;
                do {
                    Arrangement::Halfedge_const_handle e = e_it;
                    if (!std::binary_search(poly_face_ids.begin(), poly_face_ids.end(), e->twin()->face()->data())) {
                        curr = e;
                        is_on_boundary = true;
                        break;
                    }
                } while (++e_it != e_loop);

                //if edge on boundary is found, begin trace
                if (is_on_boundary) {

                    //prepare polygon to push back into solution
                    Polygon poly;

                    //remember max x value of polygon, to check if it intersects the limit
                    double max_x_value = std::numeric_limits<double>::lowest();

                    auto loop = curr;
                    do {
                        if(to_double(curr->source()->point().x()) > max_x_value) max_x_value = to_double(
                                    curr->source()->point().x());
                        //push back next point
                        poly.push_back(curr->source()->point());
                        //continue trace of the boundary
                        curr = curr->next();
                        while (std::binary_search(poly_face_ids.begin(), poly_face_ids.end(),
                                                  curr->twin()->face()->data())) {
                            curr = curr->twin()->next();
                            //mark face as visited
                            visited_faces[curr->face()->data()] = true;
                        }

                        //cout << "LINESTRING((" << std::fixed << to_double(curr->source()->point().x()) << " " << to_double(curr->source()->point().y()) << " , " << to_double(curr->target()->point().x()) << " " << to_double(curr->target()->point().y()) << "))" << endl;

                    } while (curr != loop);

                    //only holes will get found in clockwise order, we want to ignore them
                    if (poly.is_counterclockwise_oriented()) {
                        if(max_x_value>=interval_borders[interval-1]) {
                            //is intersecting divison line, needs to be considered again in the next loop
                            tmp_polys.push_back(poly);
                        }
                        else {
                            merged_polys.push_back(Polygon_wh(poly));
                        }
                    }

                }
            }
        }
    }
    return merged_polys;
}