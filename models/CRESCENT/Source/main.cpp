#include <AMReX.H>
#include <AMReX_ParmParse.H>

#include <iostream>
#include <sys/stat.h>

#include "CRESCENT.H"

namespace {

bool inputPathExists(const char* path)
{
    struct stat path_status;
    return ::stat(path, &path_status) == 0;
}

void applyAmrDefaultInputs()
{
    // Manually set defaults for certain AMReX options
    amrex::ParmParse pp("amr");

    pp.add("refine_grid_layout_x", 1);
    pp.add("refine_grid_layout_y", 1);
    pp.add("refine_grid_layout_z", 0);

    pp.add("n_proper", 2);

    int max_grid_size = 2048;
    pp.queryAdd("max_grid_size", max_grid_size);
    
    int blocking_factor = 1;
    pp.queryAdd("blocking_factor", blocking_factor);
    
    int n_error_buf = 0;
    pp.queryAdd("n_error_buf", n_error_buf);
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << " No input file provided!" << std::endl;
        return 1;
    }

    if (!inputPathExists(argv[1])) {
        std::cerr << " Input file '" << std::string(argv[1]) << "' does not exist!" << std::endl;
        return 1;
    }

    amrex::Initialize(argc, argv, true, MPI_COMM_WORLD, applyAmrDefaultInputs);

    {
        // Setup the CRESCENT solver
        CRESCENT crescent;

        // Initialize the data structures
        crescent.Init();

        // Main simulation loop
        crescent.Run();
    }

    amrex::Finalize();
    return 0;
}
