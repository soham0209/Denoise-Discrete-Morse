#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <tuple>
#include <boost/timer/progress_display.hpp>
#include <chrono>
#include <algorithm>
#include "UF.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K>    Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<K>                 Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                Tds;
//Use the Fast_location tag. Default or Compact_location works too.
typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Triangulation;
typedef Triangulation::Point_3 Point;
typedef Triangulation::Finite_cells_iterator Finite_cells_iterator;
typedef Triangulation::Finite_facets_iterator Finite_facets_iterator;

std::map<Triangulation::Cell_handle, int> id;
std::map<int, Triangulation::Cell_handle> rev_id;
std::map<unsigned, double> f;
std::map<int, double> f_val;
double getFiltration(Triangulation::Facet fi);
double getFiltrationC(Triangulation::Cell_handle c);
unsigned getindices(Triangulation::Facet fi);
void printInfo(std::ofstream &file_out, Triangulation::Facet fi, const int &root, const int &merged, const double persistence);
void writeFaces(const std::string &f_out, UF *reg);
Triangulation* T;

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        std::cerr << "Usage ./Region3D <dataset> <delta>" << std::endl;
        return -1;
    }
    const std::string data = argv[1];
    const double delta = std::stod(argv[2]);
    const std::string vert_file = data + "/" + data + "_vert.txt";
    const std::string pers_file = data + "/" + data + ".txt";
    const std::string bd_file = data + "/" + data + "_regions" + argv[2] + ".txt";
    const std::string tri_file = data + "/" + data + "_triangulated.txt";
    const std::string ind_file = data + "/" + data + "_indices";
    std::ifstream t_T(tri_file);
    std::ifstream ind(ind_file);
    
    std::vector<std::tuple <int, int, double>> persistence;  
    double nx = 0.0, ny = 0.0, nz = 0.0;
    Triangulation T_temp;
    if (t_T.good() && ind.good())
    {   std::cout << "Reading saved triangulation "<< std::endl;
        t_T >> T_temp;
        T = &T_temp;
        unsigned i;
        for (auto vit = T->finite_vertices_begin(); vit != T->finite_vertices_end(); vit++)
        {
            ind.read(reinterpret_cast<char*>(&i), sizeof(i));
			vit->info() = i;
            double a = CGAL::to_double(vit->point().x());
            double b = CGAL::to_double(vit->point().y());
            double c = CGAL::to_double(vit->point().z());
             if (a > nx)
                nx = a;
            if (b > ny)
                ny = b;
            if (c > nz)
                nz = c;
        }
        std::cout << "Triangulation read " << std::endl;
        std::cout << "Triangulation contains " << T->number_of_vertices() << " points" << std::endl;
        t_T.close();
        ind.close();
    }
    else
    {
        std::ifstream inVert(vert_file);
        std::vector<std::pair<Point, unsigned>> pts;
        unsigned i = 0;

        // Reading vertices data
        double a, b, c, w;
        
        while (inVert >> a >> b >> c)
        {
            Point p(a, b, c);
            
            pts.emplace_back(p, i);
             if (a > nx)
                nx = a;
            if (b > ny)
                ny = b;
            if (c > nz)
                nz = c;
        
            i++;
        }
        inVert.close();
        std::cout << "Read " << i << " points from file" << std::endl;
        std::cout << "Triangulating points " << std::endl;
        T = new Triangulation(pts.begin(), pts.end());
        std::cout << "Saving Triangulation for future use " << std::endl;
        std::ofstream outT(tri_file);
        outT << *T;
        std::ofstream outind(ind_file);
        for(auto vit = T->finite_vertices_begin();vit!=T->finite_vertices_end();vit++){
            unsigned x = vit->info();
            outind.write(reinterpret_cast<char*>(&x),sizeof(x));
        }
        std::cout << "Done\n";
    }
    std::cout << "Reading filtration data ..." << std::endl;
    // Reading filtration data
    std::ifstream pers(pers_file);

    std::string temp;
    pers >> temp;
    int dim = std::stoi(temp);
    int total_vert = 1;
    for (int d = 0; d < dim; d++)
    {
        pers >> temp;
        int each_dim = std::stoi(temp);
        total_vert = total_vert * each_dim;
    }
    std::cout << "Total ver on T " << T->number_of_vertices() << " total val read "
    << total_vert << std::endl; 
    assert(total_vert >= T->number_of_vertices());
    for (unsigned i = 0; i < total_vert; i++)
    {
        pers >> temp;
        //double offset = (double)i/(double)((nx+1.0)*(ny+1.0));
        f[i] = -std::stod(temp); //- offset;
    }
    pers.close();
    std::vector<std::pair<Triangulation::Facet, double>> faces;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (auto fit = T->finite_facets_begin(); fit != T->finite_facets_end(); fit++)
    {
        faces.emplace_back(std::make_pair(*fit, getFiltration(*fit)));
    }
    int id_cell = 0;
    UF *pers_tree = new UF(T->number_of_finite_cells()+1);
    UF *reg_tree = new UF(T->number_of_finite_cells()+1);
    pers_tree->make(-1);
    reg_tree->make(-1);
    f_val[-1] = std::numeric_limits<double>::max();
    
    for (Finite_cells_iterator cit = T->finite_cells_begin(); cit != T->finite_cells_end(); cit++)
    {   
        id[cit] = id_cell;
        rev_id[id_cell] = cit;
        pers_tree->make(id_cell);
        reg_tree->make(id_cell);
        f_val[id_cell] = getFiltrationC(cit);
        id_cell++;
    }

    pers_tree->f_val = &f_val;
    reg_tree->f_val = &f_val;

      std::sort(faces.begin(), faces.end(), [](std::pair<Triangulation::Facet, double> a, std::pair<Triangulation::Facet, double> b) {
        if (a.second > b.second)
            return true;
        else if (a.second == b.second){
            unsigned ind_a = getindices(a.first);
            unsigned ind_b = getindices(b.first);
            return ind_a > ind_b;
        }
        return false;
    });
    boost::timer::progress_display *show_progress = nullptr;
    show_progress = new boost::timer::progress_display(faces.size());
    unsigned num_simplification = 0;

   for (const auto fi : faces)
    {
        ++(*show_progress);
        Triangulation::Facet face = fi.first;
        double val = fi.second;
        int id_c = T->is_infinite(face.first) ? -1 : id[face.first];
        int id_cnbr = T->is_infinite(face.first->neighbor(face.second)) ? -1 : id[face.first->neighbor(face.second)];
        int root = pers_tree->root(id_c);
        int merged = pers_tree->root(id_cnbr);
        if (root == merged)
                 continue;


      if(f_val[merged] > f_val[root]){
         std::swap(merged, root);
      }
      else if (std::abs(f_val[merged]-f_val[root]) < 1e-16){
          if (merged > root)
              std::swap(merged, root);
      }
    
        double pers_val = std::abs(f_val[merged] - val);


        if (pers_val < delta)
        {
            int root_reg = T->is_infinite(face.first) ? -1 : reg_tree->root(id[face.first]);
            int merged_reg = T->is_infinite(face.first->neighbor(face.second)) ? -1 : reg_tree->root(id[face.first->neighbor(face.second)]);
            if(f_val[merged_reg] > f_val[root_reg]){
                std::swap(merged_reg, root_reg);
            }
            else if (std::abs(f_val[merged_reg]-f_val[root_reg]) < 1e-16){
                if (merged_reg > root_reg)
                    std::swap(merged_reg, root_reg);
            }
            reg_tree->merge_y_to_x(root_reg, merged_reg);
            ++num_simplification;
        
        }
        pers_tree->merge_y_to_x(root, merged);
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "Number of regions " << reg_tree->cnt << "\nNumber of components " << pers_tree->cnt << std::endl;
    std::cout << "Writing ... " << std::endl;
    writeFaces(bd_file, reg_tree);
    std::cout << "Time difference (sec) = " <<  (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) /1000000.0  <<std::endl;

}
double getFiltration(Triangulation::Facet fi)
{
    int ind = fi.second;
    unsigned v0 = fi.first->vertex((ind + 1) % 4)->info();
    unsigned v1 = fi.first->vertex((ind + 2) % 4)->info();
    unsigned v2 = fi.first->vertex((ind + 3) % 4)->info();
    return std::max(f[v0], std::max(f[v1], f[v2]));
}
double getFiltrationC(Triangulation::Cell_handle c)
{
    std::vector<double> vals;
    vals.reserve(4);
for (int i = 0; i < 4; i++)
    {
        vals.push_back(f[c->vertex(i)->info()]);
       
    }
    return *std::max_element(vals.begin(), vals.end());
}
unsigned getindices(Triangulation::Facet fi){
    Triangulation::Cell_handle c = fi.first;
    int ind = fi.second;
    unsigned v0 = c->vertex((ind+1)%4)->info();
    unsigned v1 = c->vertex((ind+2)%4)->info();
    unsigned v2 = c->vertex((ind+3)%4)->info();
    return std::max(std::max(v0, v1), v2);

}
void printInfo(std::ofstream &file_out, Triangulation::Facet fi, const int &root, const int &merged, const double persistence){
    Triangulation::Cell_handle  t = fi.first;
    int ind = fi.second;
    Triangulation::Cell_handle tnbr = t->neighbor(ind);
    int id_c = T->is_infinite(t) ? -1 : id[t];
    int id_cnbr = T->is_infinite(tnbr) ? -1 : id[tnbr];
    unsigned v0 = t->vertex((ind+1)%4)->info();
    unsigned v1 = t->vertex((ind+2)%4)->info();
    unsigned v2 = t->vertex((ind+3)%4)->info();
    file_out << " Facet vertices " << v0 << " " << v1 << " " << v2 << " f: "<< getFiltration(fi) << std::endl;
    file_out << "N0: Tetrahedron id: " << id_c << std::endl;
    if (id_c != -1){
        for(int i =0; i< 4;i++){
            file_out << t->vertex(i)->info() << " : " << t->vertex(i)->point() << std::endl;
        }
    }
    file_out << "N1: Tetrahedron id: " << id_cnbr << std::endl;
    if (id_cnbr != -1){
        for(int i =0; i< 4;i++){
            file_out << tnbr->vertex(i)->info() << " : " << tnbr->vertex(i)->point() << std::endl;
        }
    }

    file_out << "Root: Tetrahedron id: " << root << std::endl;
    if (root != -1){
        t = rev_id[root];
        for(int i =0; i< 4;i++){
            file_out << t->vertex(i)->info() << " : " << t->vertex(i)->point() << std::endl;
        }
    }
    file_out << "Merged: Tetrahedron id: " << merged << std::endl;
    if (merged != -1){
        tnbr = rev_id[merged];
        for(int i =0; i< 4;i++){
            file_out << tnbr->vertex(i)->info() << " : " << tnbr->vertex(i)->point() << std::endl;
        }
    }
    file_out << "Persistence " << persistence << std::endl;

}
void writeFaces(const std::string &f_out, UF *reg){
    std::ofstream ofile(f_out);
    boost::timer::progress_display *show_progress = nullptr;
    show_progress = new boost::timer::progress_display(std::distance(T->finite_facets_begin(), T->finite_facets_end()));
    int count = 0;
    for (Finite_facets_iterator fit = T->finite_facets_begin(); fit != T->finite_facets_end(); fit++)
    {
        ++(*show_progress);
        Triangulation::Cell_handle c = fit->first;
        Triangulation::Cell_handle cnbr = fit->first->neighbor(fit->second);

        int idr1 = T->is_infinite(c) ? -1 : id[c];
        int idr2 = T->is_infinite(cnbr) ? -1 : id[cnbr];
        auto r1 = reg->root(idr1);
        auto r2 = reg->root(idr2);
        int ii = fit->second;
        //printDict(dict);
        //std::cout << "\n\n" << idr1 << " " << idr2 << std::endl;
        //std::cout << r1 << " " << r2 << std::endl;
        if (r1 != r2)
        {
            count++;
            unsigned v0 = c->vertex((ii + 1) % 4)->info();
            unsigned v1 = c->vertex((ii + 2) % 4)->info();
            unsigned v2 = c->vertex((ii + 3) % 4)->info();
            ofile << v0 << " " << v1 << " " << v2 << std::endl;
        }
    }
    ofile.close();
    std::cout << "Wrote " << count << " faces to " << f_out << std::endl;


}
