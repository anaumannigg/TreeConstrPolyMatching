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
            //first write the original ID of the polygon
            int global_id = polygon.global_id;
            ofs.write(reinterpret_cast<const char *>(&global_id), sizeof(global_id));

            //then write the count of parts (1 means the polygon has no holes)
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

    //read and set global ID of polygon
    int global_id;
    ifs.read(reinterpret_cast<char*>(&global_id), sizeof(global_id));

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

    //set global ID and return
    polygon.global_id = global_id;
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

//write integer vector to binary file
void writeIntVecToBinaryFile(std::ofstream& out, const std::vector<int>& vec) {
    int size = vec.size();
    out.write(reinterpret_cast<const char*>(&size), sizeof(int));
    out.write(reinterpret_cast<const char*>(vec.data()), size * sizeof(int));
}

//write double vector to binary file
void writeDoubleVecToBinaryFile(std::ofstream& out, const std::vector<double>& vec) {
    int size = vec.size();
    out.write(reinterpret_cast<const char*>(&size), sizeof(int));
    out.write(reinterpret_cast<const char*>(vec.data()), size * sizeof(double));
}

//write solution object to binary file
void writeSolutionToBinaryFile(std::ofstream& out, const Solution& sol) {
    int match_size = sol.matching.size();

    out.write(reinterpret_cast<const char*>(&match_size), sizeof(int));
    for (const auto& pair : sol.matching) {
        writeIntVecToBinaryFile(out, pair.first);
        writeIntVecToBinaryFile(out, pair.second);
    }

    //write if the solution contains references
    out.write(reinterpret_cast<const char*>(&sol.has_references), sizeof(bool));
    if (sol.has_references) {
        for (const auto& pair : sol.matching_references) {
            writeIntVecToBinaryFile(out, pair.first);
            writeIntVecToBinaryFile(out, pair.second);
        }
    }

    writeIntVecToBinaryFile(out, sol.set1_match_index);
    writeDoubleVecToBinaryFile(out, sol.set1_match_weight);
    writeIntVecToBinaryFile(out, sol.set2_match_index);
    writeDoubleVecToBinaryFile(out, sol.set2_match_weight);

    out.write(reinterpret_cast<const char*>(&sol.match_count), sizeof(int));
    out.write(reinterpret_cast<const char*>(&sol.match_ids.first), sizeof(int));
    out.write(reinterpret_cast<const char*>(&sol.match_ids.second), sizeof(int));
    out.write(reinterpret_cast<const char*>(&sol.target_value), sizeof(double));

    out.write(reinterpret_cast<const char*>(&sol.numUnmatched), sizeof(int));
    out.write(reinterpret_cast<const char*>(&sol.num1to1), sizeof(int));
    out.write(reinterpret_cast<const char*>(&sol.num1toN), sizeof(int));
    out.write(reinterpret_cast<const char*>(&sol.numMtoN), sizeof(int));
    out.write(reinterpret_cast<const char*>(&sol.obj1to1), sizeof(double));
    out.write(reinterpret_cast<const char*>(&sol.obj1toN), sizeof(double));
    out.write(reinterpret_cast<const char*>(&sol.objMtoN), sizeof(double));
}

//read integer vector from binary file
std::vector<int> readIntVecFromBinaryFile(std::ifstream& in) {
    int size;
    in.read(reinterpret_cast<char*>(&size), sizeof(int));
    std::vector<int> vec(size);
    in.read(reinterpret_cast<char*>(vec.data()), size * sizeof(int));
    return vec;
}

//read double vector from binary file
std::vector<double> readDoubleVecFromBinaryFile(std::ifstream& in) {
    int size;
    in.read(reinterpret_cast<char*>(&size), sizeof(int));
    std::vector<double> vec(size);
    in.read(reinterpret_cast<char*>(vec.data()), size * sizeof(double));
    return vec;
}

//read solution from binary file
Solution readSolutionFromBinaryFile(std::ifstream& in) {
    Solution sol;
    int match_size;
    in.read(reinterpret_cast<char*>(&match_size), sizeof(int));

    sol.matching.resize(match_size);
    for (size_t i = 0; i < match_size; ++i) {
        sol.matching[i].first = readIntVecFromBinaryFile(in);
        sol.matching[i].second = readIntVecFromBinaryFile(in);
    }

    //write if the solution contains references
    bool has_references;
    in.read(reinterpret_cast<char*>(&has_references), sizeof(bool));
    sol.has_references = has_references;
    if (sol.has_references) {
        sol.matching_references.resize(match_size);
        for (size_t i = 0; i < match_size; ++i) {
            sol.matching_references[i].first = readIntVecFromBinaryFile(in);
            sol.matching_references[i].second = readIntVecFromBinaryFile(in);
        }
    }

    sol.set1_match_index = readIntVecFromBinaryFile(in);
    sol.set1_match_weight = readDoubleVecFromBinaryFile(in);
    sol.set2_match_index = readIntVecFromBinaryFile(in);
    sol.set2_match_weight = readDoubleVecFromBinaryFile(in);

    in.read(reinterpret_cast<char*>(&sol.match_count), sizeof(int));
    in.read(reinterpret_cast<char*>(&sol.match_ids.first), sizeof(int));
    in.read(reinterpret_cast<char*>(&sol.match_ids.second), sizeof(int));
    in.read(reinterpret_cast<char*>(&sol.target_value), sizeof(double));

    in.read(reinterpret_cast<char*>(&sol.numUnmatched), sizeof(int));
    in.read(reinterpret_cast<char*>(&sol.num1to1), sizeof(int));
    in.read(reinterpret_cast<char*>(&sol.num1toN), sizeof(int));
    in.read(reinterpret_cast<char*>(&sol.numMtoN), sizeof(int));
    in.read(reinterpret_cast<char*>(&sol.obj1to1), sizeof(double));
    in.read(reinterpret_cast<char*>(&sol.obj1toN), sizeof(double));
    in.read(reinterpret_cast<char*>(&sol.objMtoN), sizeof(double));

    return sol;
}

