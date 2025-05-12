#include "../include/polygon_operations.h"

namespace PolygonOperations {

    Polygon fixIfNotSimple(Polygon& p) {
        Polygon p_return = p;
        if(!p.is_simple()) {
            auto p_repaired = CGAL::Polygon_repair::repair(p);
            //every face should only be represented by a single polygon, if that cannot be done something is wrong
            assert(p_repaired.number_of_polygons_with_holes() >= 1 && "cannot repair invalid geometry, contains more than one polygon!");
            p_return = p_repaired.polygons_with_holes()[0].outer_boundary();
        }
        return p_return;
    }

    void print(Polygon& p) {
        cout << std::fixed << std::setprecision(3) << "POLYGON((";
        for(const auto& v : p.vertices()) cout << v << ", ";
        cout << "))\n";
    }
}