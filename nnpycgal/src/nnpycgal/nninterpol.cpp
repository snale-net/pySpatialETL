// Natural Neighbor Interpolation following the CGAL User Manual:
// https://doc.cgal.org/latest/Interpolation/index.html

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <stdint.h>
#include <ctime>

#include <omp.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K>    Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>      Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>        Delaunay_triangulation;
typedef Delaunay_triangulation::Vertex_handle         Vertex_handle;
typedef K::FT                                         Coord_type;
typedef K::Point_2                                    Point;
typedef std::vector< std::pair<Point, Coord_type> >   Coordinate_vector;
typedef std::map<Point, Coord_type, K::Less_xy_2>     Point_value_map;
typedef std::vector<int>                              Vint;
typedef std::vector<double>                           Vdouble;
typedef std::vector<std::vector<double>>              Vdouble2d;

// https://stackoverflow.com/questions/32102340/how-can-a-create-a-define-to-create-a-2d-vector
template<typename T>
std::vector<std::vector<T>> make_2d_vector(std::size_t rows, std::size_t cols)
{
    return std::vector<std::vector<T>>(rows, std::vector<T>(cols));
}

namespace py=pybind11;

Delaunay_triangulation triangulate(Vdouble x, Vdouble y) {

    int ntrainpoints = size(x);
    std::vector< std::pair<Point,unsigned> > points;
    points.reserve(ntrainpoints);

    Delaunay_triangulation T;
    for(int i=0; i<ntrainpoints ; i++){
        Point p(x[i], y[i]);
        points.push_back(std::make_pair(p,i));
    }
    T.insert( points.begin(),points.end() );
    return T;
}

Vdouble2d nninterpol(Delaunay_triangulation T, Vdouble z, Vdouble2d xx, Vdouble2d yy, double fill_value) {
    int dst_x_size = size(xx);
    int dst_y_size = size(xx[0]);
    Vdouble2d raster = make_2d_vector<double>(dst_x_size, dst_y_size);
    Point_value_map values;

    for (Delaunay_triangulation::Vertex_handle v : T.finite_vertex_handles()) {
          values.insert(std::make_pair(v->point(), z[v->info()]));
    }

    // Optimization (Face_handle look-up + OpenMP) following:
    // https://stackoverflow.com/questions/30354284/cgal-natural-neighbor-interpolation
    Delaunay_triangulation::Face_handle fh;
    double start = omp_get_wtime();
    #pragma omp parallel for num_threads(1) private(fh) collapse(2)
    for(int i=0; i<dst_x_size; i++) {
        for(int j=0; j<dst_y_size; j++) {
            Point p(xx[i][j], yy[i][j]);
            fh = T.locate(p, fh);
            std::vector< std::pair< Point, Coord_type > > coords;
            CGAL::Triple<std::back_insert_iterator<Coordinate_vector>, K::FT, bool> natneighbor = CGAL::natural_neighbor_coordinates_2(T, p, std::back_inserter(coords),fh);

            // error checking following:
            // https://github.com/remotesensinginfo/spdlib/blob/cf88633bd068638b13fb7701d93f01a28a8cd488/src/spd/SPDPointInterpolation.cpp
            if(!natneighbor.third)
            {
                raster[i][j] = fill_value;
            }
            else
            {
                raster[i][j] = CGAL::linear_interpolation(coords.begin(), coords.end(), natneighbor.second, CGAL::Data_access<Point_value_map>(values));
            }
        }
    }

    double end = omp_get_wtime();
    double elapsed = double(end - start);
    std::cout << "Time : " << elapsed << " seconds." << std::endl;

    return raster;
}

PYBIND11_MODULE(nninterpol, m) {
    m.doc() = "Natural Neighbor Interpolation using CGAL";
    py::class_<Delaunay_triangulation>(m, "Delaunay_triangulation");
    m.def("triangulate", &triangulate, "Triangulate the source grid using CGAL");
    m.def("nninterpol", &nninterpol,py::call_guard<py::gil_scoped_release>(), "Linear Interpolation using CGAL");
}
