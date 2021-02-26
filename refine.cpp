//
// Created by Aonymity404 on 10/30/19.
//
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/squared_distance_3.h> //for 3D functions
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <tuple>
#include <CGAL/Tetrahedron_3.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K> Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<K> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
//Use the Fast_location tag. Default or Compact_location works too.
typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Triangulation;
typedef Triangulation::Point_3 Point;
typedef Triangulation::Finite_cells_iterator Finite_cells_iterator;
typedef Triangulation::Finite_facets_iterator Finite_facets_iterator;
typedef Triangulation::Cell_handle Cell_handle;

std::map<Triangulation::Cell_handle, int> id;
std::map<unsigned, double> f;


int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage ./Refine <dataset> " << std::endl;
        return -1;
    }
    std::string data = argv[1];

    std::string vert_file = data + "/" + data + "_vert_Q.txt";
    std::string pers_file = data + "/" + data + ".txt";
    std::string out_vert_file = data + "/" + data + "_vert_refined.txt";

    std::vector<std::tuple<int, int, double>> persistence;
    double nx = 0.0, ny = 0.0, nz = 0.0;

    std::ifstream inVert(vert_file);
    std::vector<Point> pts, newpts;
    unsigned i = 0;

    // Reading vertices data
    double a, b, c;

    while (inVert >> a >> b >> c)
    {
        Point p(a, b, c);

        pts.push_back(p);
        if (a > nx)
            nx = a;
        if (b > ny)
            ny = b;
        if (c > nz)
            nz = c;

        i++;
    }
    inVert.close();
    std::cout << "# Points " << pts.size() << std::endl;
    Triangulation T(pts.begin(), pts.end());
    std::cout << "Vertices: " << T.number_of_vertices() << " Cells: "
              << T.number_of_finite_cells() << std::endl;

    std::vector<std::pair<Triangulation::Facet, Point> >face_pts;
    for(auto fit = T.finite_facets_begin();fit!=T.finite_facets_end();++fit){
        Triangulation::Cell_handle cit = fit->first;
        int ind = fit->second;
        Point u = cit->vertex((ind + 1) % 4)->point();
        Point v = cit->vertex((ind + 2) % 4)->point();
        Point w = cit->vertex((ind + 3) % 4)->point();
        const Triangulation::Triangle s(u, v, w);
        K::Vector_3 k(Point(0.0001, 0.0001, 0.0001), CGAL::ORIGIN);
        Point p = CGAL::circumcenter(s);

        face_pts.emplace_back(std::make_pair(*fit, p));

    }
    for(auto fp: face_pts){
        T.insert_in_facet(std::get<1>(fp), std::get<0>(fp));
    }
    std::cout << "Points added in facet " << face_pts.size() << " Number of vertices: " << T.number_of_vertices()
        << std::endl;
    std::cout << "Vertices: " << T.number_of_vertices() << " Cells: "
              << T.number_of_finite_cells() << std::endl;
    std::ofstream outvert(out_vert_file);
    for (auto vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); vit++)
    {
        outvert << vit->point() << std::endl;
    }
    outvert.close();
}

