#ifndef TCPOLYGONMATCHING_POLYGON_OPERATIONS_H
#define TCPOLYGONMATCHING_POLYGON_OPERATIONS_H
#include "cgal_includes.h"

namespace PolygonOperations {
    Polygon fixIfNotSimple(Polygon& p);

    void print(Polygon& p);

}



#endif //TCPOLYGONMATCHING_POLYGON_OPERATIONS_H
