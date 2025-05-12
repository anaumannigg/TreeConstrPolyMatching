#include "../include/polygon_wh.h"

//CLASS POLYGON_WH

Polygon& Polygon_wh::outer_boundary() {
    return this->outer_boundary_polygon;
}

const Polygon& Polygon_wh::outer_boundary() const {
    return this->outer_boundary_polygon;
}

const std::vector<Polygon>& Polygon_wh::holes() const {
    return this->hole_polygons;
}

std::vector<Polygon>& Polygon_wh::holes() {
    return this->hole_polygons;
}

const int Polygon_wh::number_of_holes() const {
    return this->hole_polygons.size();
}

std::vector<Polygon>::const_iterator Polygon_wh::holes_begin() const {
    return this->hole_polygons.begin();
}

std::vector<Polygon>::const_iterator Polygon_wh::holes_end() const {
    return this->hole_polygons.end();
}

void Polygon_wh::add_hole(const Polygon& hole) {
    this->hole_polygons.push_back(hole);
}

bool Polygon_wh::has_on_bounded_side(Point p) {
    if (this->outer_boundary().has_on_bounded_side(p)) {
        for (const auto& h : this->hole_polygons) {
            if (h.has_on_bounded_side(p)) {
                return false;
            }
        }
        return true;
    }
    return false;
}

bool Polygon_wh::does_intersect(const Polygon_wh& p) const {
    //first check, if the outer boundaries do intersect
    if (!CGAL::do_intersect(this->outer_boundary(),p.outer_boundary())) {
        return false;
    }

    bool this_is_larger = this->outer_boundary().area() > p.outer_boundary().area();

    if (this_is_larger) {
        if (this->lies_entirely_within_hole(p.outer_boundary()))
            return false;
    }
    else {
        if (p.lies_entirely_within_hole(this->outer_boundary()))
            return false;
    }


    return true;
}

bool Polygon_wh::does_intersect(Polygon& p) {
    Polygon_wh p_wh(p);
    return this->does_intersect(p_wh);
}


bool Polygon_wh::lies_entirely_within_hole(const Polygon& p) const {
    //check all vertices for containment, only holes that contain
    //the vertices are candidates
    std::vector<const Polygon*> hole_candidates;
    for (const auto& h : this->holes()) {
        bool all_vertices_are_contained = true;
        for (const auto& v : p.vertices()) {
            if (h.has_on_unbounded_side(v)) {
                all_vertices_are_contained = false;
                break;
            }
        }
        if (all_vertices_are_contained) hole_candidates.push_back(&h);
    }
    if (hole_candidates.empty()) return false;

    //now it suffices for at least one candidate to have no intersecting edges

    //all vertices are contained, make sure no edge intersects the polygon
    for (const auto& e : p.edges()) {
        double e_minx = e.bbox().xmin(), e_miny = e.bbox().ymin();
        double e_maxx = e.bbox().xmax(), e_maxy = e.bbox().ymax();

        for (const auto& h : hole_candidates) {
            bool is_hole_candidate_intersected = false;

            //for efficiency, first check with bbox
            double h_minx = h->bbox().xmin(), h_miny = h->bbox().ymin();
            double h_maxx = h->bbox().xmax(), h_maxy = h->bbox().ymax();

            if (e_maxx < h_minx ||
                e_maxy < h_miny ||
                e_minx > h_maxx ||
                e_miny > h_maxy)
                return true;


            for (const auto& e_h : h->edges()) {
                if (CGAL::do_intersect(e_h, e)) {
                    is_hole_candidate_intersected = true;
                }
            }
            if (!is_hole_candidate_intersected) return true;
        }
    }

    return false;
}

double Polygon_wh::intersection_area(const Polygon_wh& p) {
    //init intersection area
    double area = 0.0;

    //first get intersection area of outer boundaries
    std::vector<CGAL_Polygon_wh> inter;
    CGAL::intersection(this->outer_boundary(),p.outer_boundary(),std::back_inserter(inter));
    for (const auto& inter_p : inter) {
        area += to_double(inter_p.outer_boundary().area());
        for (const auto& h : inter_p.holes()) {
            area += to_double(h.area());
        }
    }

    //subtract all hole areas (considering hole areas occuring in both polygons
    for (const auto& h : this->holes()) {
        area -= to_double(h.area());

        //check for intersections with holes of the other polygon, add intersection area as it will
        //be subtracted twice
        for (const auto& h_p : p.holes()) {
            if (CGAL::do_intersect(h_p,h)) {
                std::vector<CGAL_Polygon_wh> inter_h;
                CGAL::intersection(h_p,h_p,std::back_inserter(inter_h));
                for (const auto& inter_p_h : inter_h) {
                    area += to_double(inter_p_h.outer_boundary().area());
                    for (const auto& h_p_h : inter_p_h.holes()) {
                        area -= to_double(h_p_h.area());
                    }
                }
            }
        }
    }

    for (const auto & h : p.holes()) {
        area += to_double(h.area());
    }

    return area;
}

double Polygon_wh::hausdorff_distance(const std::vector<const Polygon_wh*>& polys1, const std::vector<const Polygon_wh*>& polys2) {
    double hausdorff = 0.0;

    //create CGAL objects to compute union
    std::vector<CGAL_Polygon_wh> union_input1,union_input2;
    for (auto& p : polys1) {union_input1.emplace_back(p->outer_boundary());}
    for (auto& p : polys2) {union_input2.emplace_back(p->outer_boundary());}

    std::vector<CGAL_Polygon_wh> union1,union2;
    CGAL::join(union_input1.begin(), union_input1.end(), std::back_inserter(union1));
    CGAL::join(union_input2.begin(), union_input2.end(), std::back_inserter(union2));

    //collect all edges per polygon set
    std::vector<Segment> polys1_segments,polys2_segments;
    for (const auto& p : union1) {
        for (const auto& e : p.outer_boundary().edges()) polys1_segments.push_back(e);
        for (const auto& h : p.holes()) for (const auto& e : h.edges()) polys1_segments.push_back(e);
    }
    for (const auto& p : union2) {
        for (const auto& e : p.outer_boundary().edges()) polys2_segments.push_back(e);
        for (const auto& h : p.holes()) for (const auto& e : h.edges()) polys2_segments.push_back(e);
    }

    //compare all edges pairwise
    for (const auto& e_this : polys1_segments) {
        for (const auto& e_p : polys2_segments) {
            double min_dist = to_double(CGAL::squared_distance(e_p,e_this));
            if (min_dist > hausdorff) hausdorff = min_dist;
        }
    }

    //comparison needs to be performed in both directions to make the hausdorff distance symmetric
    for (const auto& e_this : polys2_segments) {
        for (const auto& e_p : polys1_segments) {
            double min_dist = to_double(CGAL::squared_distance(e_p,e_this));
            if (min_dist > hausdorff) hausdorff = min_dist;
        }
    }

    return hausdorff;
}

const CGAL::Bbox_2 Polygon_wh::bbox() const {
    return this->outer_boundary_polygon.bbox();
}

bool Polygon_wh::is_valid() {
    if (!this->outer_boundary_polygon.is_simple()) return false;
    for (const auto& h : this->hole_polygons) {
        if (!h.is_simple()) return false;
    }
    return true;
}

void Polygon_wh::repair() {
    //TODO
}

void Polygon_wh::clear() {
    this->outer_boundary_polygon.clear();
    this->hole_polygons.clear();

    auto b = this->outer_boundary_polygon.bbox();
}


//END CLASS POLYGON_WH