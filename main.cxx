#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"
#include "binaryio.hpp"
#include "fields.hpp"
#include "geometry.hpp"
#include "ic.hpp"
#include "input.hpp"
#include "matprops.hpp"
#include "mesh.hpp"
#include "output.hpp"
#include "rheology.hpp"

#ifdef WIN32
#ifdef _MSC_VER
#define snprintf _snprintf
#endif // _MSC_VER
namespace std { using ::snprintf; }
#endif // WIN32

void init_var(const Param& param, Variables& var)
{
    var.time = 0;
    var.steps = 0;
}


void init(const Param& param, Variables& var)
{
    std::cout << "Initializing mesh and field data...\n";

    create_new_mesh(param, var);
    create_boundary_flags(var);
    create_boundary_nodes(var);
    create_boundary_facets(var);
    create_support(var);

    allocate_variables(param, var);

    compute_volume(*var.coord, *var.connectivity, *var.volume);

    // initialize shape functions here
    // compute_shape_fn(...)

    // apply_bcs(param, var, *var.vel);

    // temperature should be init'd before stress and strain
    initial_temperature(param, var, *var.temperature);
}


int main(int argc, const char* argv[])
{
    double start_time = 0;

    //
    // read command line
    //
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " config_file\n";
        std::cout << "       " << argv[0] << " -h or --help\n";
        return -1;
    }

    Param param;
    get_input_parameters(argv[1], param);

    //
    // run simulation
    //
    static Variables var; // declared as static to silence valgrind's memory leak detection
    init_var(param, var);

    Output output(param, start_time,
                  (param.sim.is_restarting) ? param.sim.restarting_from_frame : 0);


    init(param, var);


    var.dt = compute_dt(param, var);
    output.write(var, false);

    double starting_time = var.time; // var.time & var.steps might be set in restart()
    double starting_step = var.steps;
    int next_regular_frame = 1;  // excluding frames due to output_during_remeshing

    std::cout << "Starting simulation...\n";
    do {
        var.steps ++;
        var.time += var.dt;

    
        update_temperature(param, var, *var.temperature, *var.ntmp);

        
        if (
            ( (var.steps - starting_step) == next_regular_frame * param.sim.output_step_interval ) ||
            ( (var.time - starting_time) > next_regular_frame * param.sim.output_time_interval_in_yr * YEAR2SEC )
           )
        {
            output.write(var);

            next_regular_frame ++;
        }
    } while (var.steps < param.sim.max_steps && var.time <= param.sim.max_time_in_yr * YEAR2SEC);

    std::cout << "Ending simulation.\n";
    return 0;
}
