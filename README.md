# TreeConstrPolyMatching
A C++ implementation for the matching of polygons based on a tree-constrained bipartite matching solution.

## Project Setup

Note that this projects requires the shapelib (https://github.com/OSGeo/shapelib). It should be cloned and placed into a subfolder `lib/shapelib` and have a subdirectory `build`, where it is built.

Further, CGAL (6.0.1), Boost (1.83.0) and Gurobi (12.0) are required libraries that should be installed on your system.

After that the project can be built by navigating to the project folder and executing:
```bash 
mkdir build && cd ./build
cmake ..
make
``` 

## Usage

After making the project, it can be executed by:
```bash 
./TCPolygonMatching 
``` 

It supports multiple Command Line Options.

### Necessary Command Line Options
1. `-d`: specify the `dataset_name`. The code expects it to be located within the directory `input/dataset_name` as two shapefiles `dataset_name_1.shp`and `dataset_name_2.shp`.
2. Specify the value of $\lambda \in [0,1]$ used in the objective $f(\mu) = \mathrm{IoU}(\mu) - \lambda$ 
    - `-l`: allows to specify one value for lambda
    - `-lr`: allows to specify a range of multiple lambda values `start stepsize end`

### Optional Command Line Options
1. `-s`: activates exploit of properties of optimal solutions in preprocessing, **recommended**.
2. `-t`: specify the number of threads used by the code, default is 1.
3. `-r`: specify the mode that is used to build the hierachical groupings (informed,kruskal), default is informed.
4. `-m`: specify the solution mode (opt,3approx), default is 3approx.
5. `-e`: specify a name for the log-file generated on output, default is log.

## Input Data
For our experiments, we used input data from the following sources:
- https://download.geofabrik.de/ for OSM data
- https://www.opengeodata.nrw.de/produkte/geobasis/lk/akt/gru_xml/ for cadastral data
- https://github.com/microsoft/GlobalMLBuildingFootprints for AI generated building footprints

## Further remarks
Note that the code requires shapefiles containing **single** polygons as inputs. Multipolygons are not supported. Downloaded data in other formats needs to be converted.
We recommend to use GIS-Tools to repair invalid geometries before using this code, as those might lead to instabilities.
The code can be provided a third dataset `input/dataset_name/dataset_name_merged.shp`. This is used to parallelize the decomposition of the input polygons into connected components. 
It should contain the dissolved union of both input datasets. If not provided, the code will initially compute this union itself. Note that due to invalid geometries, this
might not always be stable. We recommend to use a GIS-Tool to fix problematic geometries after the computation by the code. The `merged` file needs to be computed only once per input pair of datasets.
