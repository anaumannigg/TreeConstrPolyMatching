#ifndef TCPOLYGONMATCHING_POLYGON_OPERATIONS_H
#define TCPOLYGONMATCHING_POLYGON_OPERATIONS_H
#include "cgal_includes.h"

namespace PolygonOperations {
    // Declare the function to buffer a hole inward
    void buffer_holes(Polygon_wh& p, double epsilon = 1e-9);

    void buffer(Polygon& p, double epsilon=1e-4);

    Polygon fixIfNotSimple(Polygon& p);

    void print(Polygon& p);

}



#endif //TCPOLYGONMATCHING_POLYGON_OPERATIONS_H
