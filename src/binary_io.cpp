#include "../include/binary_io.h"


void writePolysToBinaryFile(const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2 , const std::string& file_path, size_t set_id) {
    std::ofstream ofs(file_path, std::ios::binary | std::ios::app);
    if (!ofs) {
        throw std::runtime_error("Failed to open file for writing");
    }
    // Write set ID
    ofs.write(reinterpret_cast<const char*>(&set_id), sizeof(set_id));

    //for both sets of polygons
    for(const auto& polys : {polys1, polys2}) {
        // Write number of polygons
        size_t polys_count = polys.size();
        ofs.write(reinterpret_cast<const char *>(&polys_count), sizeof(polys_count));
        // Write each polygon of polys
        for (const auto &polygon: polys) {
            //first write the count of parts (1 means the polygon has no holes)
            size_t part_count = 1 + polygon.number_of_holes();
            ofs.write(reinterpret_cast<const char* >(&part_count), sizeof(part_count));

            for(size_t part = 0; part < part_count; part ++ ) {
                Polygon part_poly = polygon.outer_boundary();
                if(part > 0) {
                    auto holes_it = polygon.holes_begin();
                    holes_it = std::next(holes_it,part-1);
                    part_poly = *holes_it;
                }

                size_t vertex_count = part_poly.size();
                ofs.write(reinterpret_cast<const char *>(&vertex_count), sizeof(vertex_count));
                for (const auto &vertex: part_poly.vertices()) {
                    double x = to_double(vertex.x());
                    double y = to_double(vertex.y());
                    ofs.write(reinterpret_cast<const char *>(&x), sizeof(x));
                    ofs.write(reinterpret_cast<const char *>(&y), sizeof(y));
                }
            }
        }
    }
    ofs.close();
}

Polygon_wh readPolygonFromBinaryFile(std::ifstream& ifs) {
    //init poly memory
    Polygon_wh polygon;

    size_t part_count;
    ifs.read(reinterpret_cast<char*>(&part_count),sizeof(part_count));
    for(size_t part = 0; part < part_count; part++) {
        size_t vertex_count;
        ifs.read(reinterpret_cast<char *>(&vertex_count), sizeof(vertex_count));

        Polygon poly_part;
        for (size_t i = 0; i < vertex_count; ++i) {
            double x, y;
            ifs.read(reinterpret_cast<char *>(&x), sizeof(x));
            ifs.read(reinterpret_cast<char *>(&y), sizeof(y));
            poly_part.push_back(Point(x, y));
        }

        if(part == 0) polygon = Polygon_wh(poly_part);
        else polygon.add_hole(poly_part);
    }

    return polygon;
}

std::pair<std::vector<Polygon_wh>,std::vector<Polygon_wh>> readPolysFromBinaryFile(const std::string& file_path, size_t set_id) {
    std::ifstream ifs(file_path, std::ios::binary);
    if (!ifs) {
        throw std::runtime_error("Failed to open file for reading");
    }


    std::vector<Polygon_wh> polys1,polys2;
    while (ifs) {
        size_t read_set_id;
        ifs.read(reinterpret_cast<char*>(&read_set_id), sizeof(read_set_id));
        if (read_set_id != set_id) {
            // Skip this set's section

            // Read the first set count
            size_t first_set_count;
            ifs.read(reinterpret_cast<char*>(&first_set_count), sizeof(first_set_count));
            // Calculate how much data to skip for the first set
            for (size_t i = 0; i < first_set_count; ++i) {
                size_t vertex_count;
                ifs.read(reinterpret_cast<char*>(&vertex_count), sizeof(vertex_count));
                size_t first_set_size = sizeof(double) * 2 * vertex_count; // 2 doubles (x, y) per vertex

                //skip the vertices
                std::streampos current_pos = ifs.tellg();
                ifs.seekg(current_pos + std::streampos(first_set_size));
            }

            //move to second set


            // Read the second set count
            size_t second_set_count;
            ifs.read(reinterpret_cast<char*>(&second_set_count), sizeof(second_set_count));
            // Calculate how much data to skip for the second set
            size_t second_set_size = 0;
            for (size_t i = 0; i < second_set_count; ++i) {
                size_t vertex_count;
                ifs.read(reinterpret_cast<char*>(&vertex_count), sizeof(vertex_count));
                size_t second_set_size = sizeof(double) * 2 * vertex_count; // 2 doubles (x, y) per vertex

                //skip the vertices
                std::streampos current_pos = ifs.tellg();
                ifs.seekg(current_pos + std::streampos(second_set_size));
            }

            continue;
        }

        // Read the first set of polygons
        size_t first_set_count;
        ifs.read(reinterpret_cast<char*>(&first_set_count), sizeof(first_set_count));
        for (size_t i = 0; i < first_set_count; ++i) {
            polys1.push_back(readPolygonFromBinaryFile(ifs));
        }

        // Read the second set of polygons
        size_t second_set_count;
        ifs.read(reinterpret_cast<char*>(&second_set_count), sizeof(second_set_count));
        for (size_t i = 0; i < second_set_count; ++i) {
            polys2.push_back(readPolygonFromBinaryFile(ifs));
        }

        // Once the data for this set_id is read, exit the loop
        break;
    }
    ifs.close();


    return {polys1,polys2};
}


//CLASS BINARYPOLYGONFILEREADER

bool BinaryPolygonFileReader::readNextSet(size_t &set_id, std::vector<Polygon_wh> &first_polygon_set,
                                           std::vector<Polygon_wh> &second_polygon_set) {
    if (!ifs || ifs.eof()) {
        return false; // End of file or error
    }


    // Read thread ID
    ifs.read(reinterpret_cast<char*>(&set_id), sizeof(set_id));
    if (ifs.eof()) return false;

    // Read first set count
    size_t first_set_count;
    ifs.read(reinterpret_cast<char*>(&first_set_count), sizeof(first_set_count));
    if (ifs.eof()) return false;

    // Read first set of polygons
    first_polygon_set.clear();
    for (size_t i = 0; i < first_set_count; ++i) {
        first_polygon_set.push_back(readPolygonFromBinaryFile(ifs));
    }

    // Read second set count
    size_t second_set_count;
    ifs.read(reinterpret_cast<char*>(&second_set_count), sizeof(second_set_count));
    if (ifs.eof()) return false;

    // Read second set of polygons
    second_polygon_set.clear();
    for (size_t i = 0; i < second_set_count; ++i) {
        second_polygon_set.push_back(readPolygonFromBinaryFile(ifs));
    }

    return true; // Successfully read one line
}

//END CLASS BINARYPOLYGONFILEREADER