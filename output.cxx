#include <algorithm>  // For std::max_element
#include <cmath>
#include <cstdio>
#include <iterator>  // For std::distance
#include <iostream>

#ifdef USE_OMP
#include <omp.h>
#else
#include <ctime>
#endif

#include "constants.hpp"
#include "parameters.hpp"
#include "binaryio.hpp"
#include "geometry.hpp"
#include "matprops.hpp"
#include "output.hpp"

#ifdef WIN32
#ifdef _MSC_VER
#define snprintf _snprintf
#endif // _MSC_VER
namespace std { using ::snprintf; }
#endif // WIN32

Output::Output(const Param& param, double start_time, int start_frame) :
    modelname(param.sim.modelname),
    start_time(start_time),
    average_interval(param.sim.output_averaged_fields),
    has_marker_output(param.sim.has_marker_output),
    frame(start_frame),
    time0(0)
{}


Output::~Output()
{}


void Output::write_info(const Variables& var, double dt)
{
#ifdef USE_OMP
    double run_time = omp_get_wtime() - start_time;
#else
    double run_time = double(std::clock()) / CLOCKS_PER_SEC;
#endif
    char buffer[256];
    std::snprintf(buffer, 255, "%6d\t%10d\t%12.6e\t%12.4e\t%12.6e\t%8d\t%8d\t%8d\n",
                  frame, var.steps, var.time, dt, run_time,
                  var.nnode, var.nelem, var.nseg);

    std::string filename(modelname + ".info");
    std::FILE* f;
    if (frame == 0)
        f = std::fopen(filename.c_str(), "w");
    else
        f = std::fopen(filename.c_str(), "a");

    if (f == NULL) {
        std::cerr << "Error: cannot open file '" << filename << "' for writing\n";
        std::exit(2);
    }

    if (std::fputs(buffer, f) == EOF) {
        std::cerr << "Error: failed writing to file '" << filename << "'\n";
        std::cerr << "\tbuffer written:\n";
        std::cerr << buffer << '\n';
        std::exit(2);
    }

    std::fclose(f);
}


void Output::write(const Variables& var, bool is_averaged)
{
    double dt = var.dt;
    
    write_info(var, dt);

    char filename[256];
    std::snprintf(filename, 255, "%s.save.%06d", modelname.c_str(), frame);
    BinaryOutput bin(filename);

    bin.write_array(*var.coord, "coordinate", var.coord->size());
    bin.write_array(*var.connectivity, "connectivity", var.connectivity->size());

    bin.write_array(*var.temperature, "temperature", var.temperature->size());

    bin.close();
    std::cout << "  Output # " << frame
              << ", step = " << var.steps
              << ", time = " << var.time / YEAR2SEC << " yr"
              << ", dt = " << dt / YEAR2SEC << " yr.\n";

    frame ++;
}
