#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <sstream>

#define REAL double
#define VOID void
#define ANSI_DECLARATORS
#include "triangle/triangle.h"
#undef REAL
#undef VOID
#undef ANSI_DECLARATORS


#include "constants.hpp"
#include "parameters.hpp"
#include "utils.hpp"
#include "mesh.hpp"

#ifdef WIN32
#ifndef _MSC_VER
namespace std {
  static std::string to_string(long double t) 
  {
    char temp[32];
    sprintf(temp,"%f",double(t));
    return std::string(temp);
  }
}
#endif //_MSC_VER
#endif // WIN32


namespace { // anonymous namespace

void set_verbosity_str(std::string &verbosity, int meshing_verbosity)
{
    switch (meshing_verbosity) {
    case -1:
        verbosity = "Q";
        break;
    case 0:
        verbosity = "";
        break;
    case 1:
        verbosity = "V";
        break;
    case 2:
        verbosity = "VV";
        break;
    case 3:
        verbosity = "VVV";
        break;
    default:
        verbosity = "";
        break;
    }
}


void set_volume_str(std::string &vol, double max_volume)
{
    vol.clear();
    if (max_volume > 0) {
        vol += 'a';
        vol += std::to_string((long double)max_volume);
    }
    else if (max_volume == 0) {
        // max_volume is set inside regional attributes (usu. read from poly files)
        vol += 'a';
    }
}


void set_2d_quality_str(std::string &quality, double min_angle)
{
    quality.clear();
    if (min_angle > 0) {
        quality += 'q';
        quality += std::to_string((long double)min_angle);
    }
}


void triangulate_polygon
(double min_angle, double max_area,
 int meshing_verbosity,
 int npoints, int nsegments,
 const double *points, const int *segments, const int *segflags,
 const int nregions, const double *regionattributes,
 int *noutpoints, int *ntriangles, int *noutsegments,
 double **outpoints, int **triangles,
 int **outsegments, int **outsegflags, double **outregattr)
{
    char options[255];
    triangulateio in, out;

    std::string verbosity, vol, quality;
    set_verbosity_str(verbosity, meshing_verbosity);
    set_volume_str(vol, max_area);
    set_2d_quality_str(quality, min_angle);

    if( nregions > 0 )
        std::sprintf(options, "%s%spjz%sA", verbosity.c_str(), quality.c_str(), vol.c_str());
    else
        std::sprintf(options, "%s%spjz%s", verbosity.c_str(), quality.c_str(), vol.c_str());

    if( meshing_verbosity >= 0 )
        std::cout << "The meshing option is: " << options << '\n';

    in.pointlist = const_cast<double*>(points);
    in.pointattributelist = NULL;
    in.pointmarkerlist = NULL;
    in.numberofpoints = npoints;
    in.numberofpointattributes = 0;

    in.trianglelist = NULL;
    in.triangleattributelist = NULL;
    in.trianglearealist = NULL;
    in.numberoftriangles = 0;
    in.numberofcorners = 3;
    in.numberoftriangleattributes = 0;

    in.segmentlist = const_cast<int*>(segments);
    in.segmentmarkerlist = const_cast<int*>(segflags);
    in.numberofsegments = nsegments;

    in.holelist = NULL;
    in.numberofholes = 0;

    in.numberofregions = nregions;
    if( nregions > 0 )
        in.regionlist = const_cast<double*>(regionattributes);
    else
        in.regionlist = NULL;

    out.pointlist = NULL;
    out.pointattributelist = NULL;
    out.pointmarkerlist = NULL;
    out.trianglelist = NULL;
    out.triangleattributelist = NULL;
    out.neighborlist = NULL;
    out.segmentlist = NULL;
    out.segmentmarkerlist = NULL;
    out.edgelist = NULL;
    out.edgemarkerlist = NULL;

    /*******************************/
    triangulate(options, &in, &out, NULL);
    /*******************************/

    *noutpoints = out.numberofpoints;
    *outpoints = out.pointlist;

    *ntriangles = out.numberoftriangles;
    *triangles = out.trianglelist;

    *noutsegments = out.numberofsegments;
    *outsegments = out.segmentlist;
    *outsegflags = out.segmentmarkerlist;
    *outregattr = out.triangleattributelist;

    trifree(out.pointmarkerlist);
}


void points_to_mesh(const Param &param, Variables &var,
                    int npoints, const double *points,
                    int n_init_segments, const int *init_segments, const int *init_segflags,
                    int nregions, const double *regattr,
                    double max_elem_size, int vertex_per_polygon)
{
    double *pcoord, *pregattr;
    int *pconnectivity, *psegment, *psegflag;

    points_to_new_mesh(param.mesh, npoints, points,
                       n_init_segments, init_segments, init_segflags,
                       nregions, regattr,
                       max_elem_size, vertex_per_polygon,
                       var.nnode, var.nelem, var.nseg,
                       pcoord, pconnectivity, psegment, psegflag, pregattr);

    var.coord = new array_t(pcoord, var.nnode);
    var.connectivity = new conn_t(pconnectivity, var.nelem);
    var.segment = new segment_t(psegment, var.nseg);
    var.segflag = new segflag_t(psegflag, var.nseg);
    var.regattr = new regattr_t(pregattr, var.nelem);
}


void my_fgets(char *buffer, int size, std::FILE *fp,
              int &lineno, const std::string &filename)
{
    char *s;
    while (1) {
        ++ lineno;
        s = std::fgets(buffer, size, fp);
        if (! s) {
            std::cerr << "Error: reading line " << lineno
                      << " of '" << filename << "'\n";
            std::exit(2);
        }

        // check for blank lines and comments
        if (buffer[0] != '\n' && buffer[0] != '#') break;
    }
}


void new_mesh_from_polyfile(const Param& param, Variables& var)
{
    /* The format specifiction for the poly file can be found in:
     *   2D:  http://www.cs.cmu.edu/~quake/triangle.poly.html
     *   3D:  http://wias-berlin.de/software/tetgen/fformats.poly.html
     *
     * Note that the poly file in 3D has a more complicated format.
     */

    const double std_elem_size = 1.5 * param.mesh.resolution * param.mesh.resolution;

    std::FILE *fp = std::fopen(param.mesh.poly_filename.c_str(), "r");
    if (! fp) {
        std::cerr << "Error: Cannot open poly_filename '" << param.mesh.poly_filename << "'\n";
        std::exit(2);
    }

    int lineno = 0;
    int n;
    char buffer[255];

    // get header of node list
    int npoints;
    {
        my_fgets(buffer, 255, fp, lineno, param.mesh.poly_filename);

        int dim, nattr, nbdrym;
        n = std::sscanf(buffer, "%d %d %d %d", &npoints, &dim, &nattr, &nbdrym);
        if (n != 4) {
            std::cerr << "Error: parsing line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            std::exit(1);
        }

        if (dim != NDIMS ||
            nattr != 0 ||
            nbdrym != 0) {
            std::cerr << "Error: unsupported value in line " << lineno
                      << " of '" << param.mesh.poly_filename << "'\n";
            std::exit(1);
        }
    }

    // get node list
    double *points = new double[npoints * NDIMS];
    for (int i=0; i<npoints; i++) {
        my_fgets(buffer, 255, fp, lineno, param.mesh.poly_filename);

        int k;
        double *x = &points[i*NDIMS];
        n = std::sscanf(buffer, "%d %lf %lf", &k, x, x+1);

        if (n != NDIMS+1) {
            std::cerr << "Error: parsing line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            std::exit(1);
        }
        if (k != i) {
            std::cerr << "Error: node number is continuous from 0 at line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            std::exit(1);
        }
    }

    // get header of segment (facet) list
    int n_init_segments;
    {
        my_fgets(buffer, 255, fp, lineno, param.mesh.poly_filename);

        int has_bdryflag;
        n = std::sscanf(buffer, "%d %d", &n_init_segments, &has_bdryflag);
        if (n != 2) {
            std::cerr << "Error: parsing line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            std::exit(1);
        }

        if (has_bdryflag != 1) {
            std::cerr << "Error: unsupported value in line " << lineno
                      << " of '" << param.mesh.poly_filename << "'\n";
            std::exit(1);
        }
    }

    
    int *init_segments = new int[n_init_segments * NODES_PER_FACET];
    int *init_segflags = new int[n_init_segments];
    for (int i=0; i<n_init_segments; i++) {
        my_fgets(buffer, 255, fp, lineno, param.mesh.poly_filename);

        int *x = &init_segments[i*NODES_PER_FACET];
        int junk, bdryflag;
        n = std::sscanf(buffer, "%d %d %d %d", &junk, x, x+1, &bdryflag);
        if (n != NODES_PER_FACET+2) {
            std::cerr << "Error: parsing line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            std::exit(1);
        }
        if (bdryflag == 0) goto flag_ok;
        for (int j=0; j<nbdrytypes; j++) {
            if (bdryflag == 1 << j) goto flag_ok;
        }
        std::cerr << "Error: bdry_flag has multiple bits set in line " << lineno
                  << " of '" << param.mesh.poly_filename << "'\n";
        std::exit(1);
    flag_ok:
        init_segflags[i] = bdryflag;
    }
    for (int i=0; i<n_init_segments; i++) {
        int *x = &init_segments[i*NODES_PER_FACET];
        for (int j=0; j<NODES_PER_FACET; j++) {
            if (x[j] < 0 || x[j] >= npoints) {
                std::cerr << "Error: segment contains out-of-range node # [0-" << npoints
                          <<"] in line " << lineno << " of '"
                          << param.mesh.poly_filename << "'\n";
                std::exit(1);
            }
        }
    }

    // get header of hole list
    {
        my_fgets(buffer, 255, fp, lineno, param.mesh.poly_filename);

        int nholes;
        n = std::sscanf(buffer, "%d", &nholes);
        if (n != 1) {
            std::cerr << "Error: parsing line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            std::exit(1);
        }

        if (nholes != 0) {
            std::cerr << "Error: unsupported value in line " << lineno
                      << " of '" << param.mesh.poly_filename << "'\n";
            std::exit(1);
        }
    }

    // get header of region list
    int nregions;
    {
        my_fgets(buffer, 255, fp, lineno, param.mesh.poly_filename);

        n = std::sscanf(buffer, "%d", &nregions);
        if (n != 1) {
            std::cerr << "Error: parsing line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            std::exit(1);
        }
        if (nregions <= 0) {
            std::cerr << "Error: nregions <= 0, at line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            std::exit(1);
        }
    }

    // get region list
    double *regattr = new double[nregions * (NDIMS+2)]; // each region has these data fields: x, (y,) z, region marker (mattype), and volume.
    bool has_max_size = false;
    for (int i=0; i<nregions; i++) {
        my_fgets(buffer, 255, fp, lineno, param.mesh.poly_filename);

        int junk;
        double *x = &regattr[i*(NDIMS+2)];
        n = std::sscanf(buffer, "%d %lf %lf %lf %lf", &junk, x, x+1, x+2, x+3);

        if (n != NDIMS+3) {
            std::cerr << "Error: parsing line " << lineno << " of '"
                      << param.mesh.poly_filename << "'. "<<NDIMS+3<<" values should be given but only "<<n<<" found.\n";
            std::exit(1);
        }

        if ( x[NDIMS] < 0 || x[NDIMS] >= param.mat.nmat ) {
            std::cerr << "Error: "<<NDIMS+2<<"-th value in line "<<lineno<<" should be >=0 and < "<<param.mat.nmat<<" (=mat.num_materials) but is "<<x[NDIMS]<<"\n";
            std::cerr << "Note that this parameter is directly used as the index of mat. prop. arrays.\n";
            std::exit(1);
        }

        if ( x[NDIMS+1] > 0 ) {
            has_max_size = true; // max area is set for this region
        }
    }

    double max_elem_size = std_elem_size;
    if ( has_max_size ) max_elem_size = 0; // special value, see set_volume_str() above.

    std::cerr << " About to start triangulation\n";
    points_to_mesh(param, var, npoints, points,
                   n_init_segments, init_segments, init_segflags, nregions, regattr,
                   max_elem_size, NODES_PER_FACET);

    std::cerr << " Finished triangulation\n";

    delete [] points;
    delete [] init_segments;
    delete [] init_segflags;
    delete [] regattr;
}

}


void points_to_new_mesh(const Mesh &mesh, int npoints, const double *points,
                        int n_init_segments, const int *init_segments, const int *init_segflags,
                        int n_regions, const double *regattr,
                        double max_elem_size, int vertex_per_polygon,
                        int &nnode, int &nelem, int &nseg, double *&pcoord,
                        int *&pconnectivity, int *&psegment, int *&psegflag, double *&pregattr)
{
    triangulate_polygon(mesh.min_angle, max_elem_size,
                        mesh.meshing_verbosity,
                        npoints, n_init_segments, points,
                        init_segments, init_segflags,
                        n_regions, regattr,
                        &nnode, &nelem, &nseg,
                        &pcoord, &pconnectivity,
                        &psegment, &psegflag, &pregattr);

    if (nelem <= 0) {
        std::cerr << "Error: triangulation failed\n";
        std::exit(10);
    }

}


void create_boundary_flags2(uint_vec &bcflag, int nseg,
                            const int *psegment, const int *psegflag)
{
    for (int i=0; i<nseg; ++i) {
        uint flag = static_cast<uint>(psegflag[i]);
        const int *n = psegment + i * NODES_PER_FACET;
        for (int j=0; j<NODES_PER_FACET; ++j) {
            bcflag[n[j]] |= flag;
        }
    }
}


void create_boundary_flags(Variables& var)
{
    // allocate and init to 0
    if (var.bcflag) delete var.bcflag;
    var.bcflag = new uint_vec(var.nnode);

    create_boundary_flags2(*var.bcflag, var.segment->size(),
                           var.segment->data(), var.segflag->data());
}


void create_boundary_nodes(Variables& var)
{
    /* var.bnodes[i] contains a list of nodes on the i-th boundary.
     * (See constants.hpp for the order of boundaries.)
     */
    for (std::size_t i=0; i<var.bcflag->size(); ++i) {
        uint f = (*var.bcflag)[i];
        for (int j=0; j<nbdrytypes; ++j) {
            if (f & (1<<j)) {
                // this node belongs to a boundary
                (var.bnodes[j]).push_back(i);
            }
        }
    }

    // for (int j=0; j<nbdrytypes; ++j) {
    //     std::cout << "boundary " << j << '\n';
    //     print(std::cout, var.bnodes[j]);
    //     std::cout << '\n';
    // }
}


namespace {

    struct OrderedInt
    {
#ifdef THREED
        int a, b, c;
        OrderedInt(int x, int y, int z)
        {
            // ensure that a < b < c
            if (x < y) {
                if (y < z)
                    a = x, b = y, c = z;
                else {
                    if (x < z)
                        a = x, b = z, c = y;
                    else
                        a = z, b = x, c = y;
                }
            }
            else {
                if (x < z)
                    a = y, b = x, c = z;
                else {
                    if (y < z)
                        a = y, b = z, c = x;
                    else
                        a = z, b = y, c = x;
                }
            }
        }

        bool operator==(OrderedInt &rhs)
        {
            return a==rhs.a && b==rhs.b && c==rhs.c;
        }

#else

        int a, b;
        OrderedInt(int x, int y)
        {
            // ensure that a < b
            if (x < y)
                a = x, b = y;
            else
                a = y, b = x;
        }

        bool operator==(OrderedInt &rhs)
        {
            return a==rhs.a && b==rhs.b;
        }
#endif
    };
}


void create_boundary_facets(Variables& var)
{
    /* var.bfacets[i] contains a list of facets (or segments in 2D)
     * on the i-th boundary. (See constants.hpp for the order of boundaries.)
     */

    // Looping through var.segment
    for (int i=0; i<var.nseg; ++i) {
        uint flag = static_cast<uint>((*var.segflag)[i][0]);
        if ((flag & BOUND_ANY) == 0) continue; // not a boundary facet
        // Nodes of this facet
        OrderedInt af((*var.segment)[i][0], (*var.segment)[i][1]);

        // Finding the corresponding element and facet #
        for (int e=0; e<var.nelem; ++e) {
            const int *conn = (*var.connectivity)[e];
            for (int f=0; f<FACETS_PER_ELEM; ++f) {
                if ((flag & (*var.bcflag)[conn[NODE_OF_FACET[f][0]]]
                    & (*var.bcflag)[conn[NODE_OF_FACET[f][1]]]
                    ) == 0U) continue; // skip

                OrderedInt bf(conn[NODE_OF_FACET[f][0]], conn[NODE_OF_FACET[f][1]]);
                if (af == bf) {
                    for (int k=0; k<nbdrytypes; ++k) {
                        if (flag == (1U << k)) {
                            var.bfacets[k].push_back(std::make_pair(e,f));
                            goto found_facet; // break out of nested loops
                        }
                    }
                }
            }
        }
        // not found
        std::cerr << "Error: " << i << "-th segment is not on any element\n";
        std::exit(12);

    found_facet:
        continue;
    }
}



void create_support(Variables& var)
{
    var.support = new std::vector<int_vec>(var.nnode);

    // create the inverse mapping of connectivity
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            (*var.support)[conn[i]].push_back(e);
        }
    }
    // std::cout << "support:\n";
    // print(std::cout, *var.support);
    // std::cout << "\n";
}


void create_new_mesh(const Param& param, Variables& var)
{
    
    new_mesh_from_polyfile(param, var);

}