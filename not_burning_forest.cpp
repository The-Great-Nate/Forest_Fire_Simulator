#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cstdlib>
#include <filesystem>
#include <mpi.h>

// Forest Fires are not good

/**
 * Struct that stores whever the forest must continue burning and/or if columns must be communicated between processes
 */
struct ForestState
{
    bool burning;
    bool comm_cols;

    ForestState(bool burning_val, bool comm_cols_val)
        : burning(burning_val), comm_cols(comm_cols_val)
    {
    }
};

/**
 * Random number generator. Only really intended for one process to run and then broadcast.
 * @return Random float between 0 and 1
 */
double rng()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    double random_float = dist(gen);
    return random_float;
}

/**
 * Generates a random string of characters for assigning randomly generated grid id's
 * @return string of length 5
 */
std::string generate_config_id()
{
    const std::string chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
    std::string id = "";

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, chars.size() - 1);

    // Generate 5 random characters/digits
    for (int i = 0; i < 5; ++i)
    {
        int rand_ind = dist(gen);
        id += chars[rand_ind];
    }

    return id;
}

/**
 * Calculates index from 2D array -> 1D equivalent
 * @return int
 */
int get_1D_index(int i, int j, int Nj)
{
    return (i * Nj) + j;
}

/**
 * Outputs grid state held by specific iproc. Helpful funcion for debugging
 * @param[in] Ni - int - No. row's of grid
 * @param[in] Nj - int - No. columns's of grid
 * @param[in] iproc - int - Process ID to output grid it has stored
 * @param[in] nproc - int - No. of processes (Essential for filtering out other processes from running this at the same time)
 * @param[in] grid - const std::vector<int> & - Grid to output. Passed by reference
 */
void display_grid(int Ni, int Nj, int iproc, int nproc, const std::vector<int> &grid)
{
    for (int p = 0; p < nproc; p++)
    {
        if (p == iproc)
        {
            for (int i = 0; i < Ni; ++i)
            {
                for (int j = 0; j < Nj; ++j)
                {
                    int ind = get_1D_index(i, j, Nj); // Calculate 1D index based off i and j
                    std::cout << grid[ind] << " ";    // Output grid
                }
                std::cout << "\n"; // New line foe each row
            }
            std::cout << "\n"; // New line for next output for tidiness
        }
    }
}

/**
 * Outputs grid state of only parts of grid managed by specific iproc. Unmanaged part set to 0
 * Mainly for write_grid() to write to file
 * @param[in] Ni - int - No. row's of grid
 * @param[in] Nj - int - No. columns's of grid
 * @param[in] iproc - int - Process ID to output grid it has stored
 * @param[in] nproc - int - No. of processes (Essential for filtering out other processes from running this at the same time)
 * @param[in] j0 - int - Left boundary index stored by iproc
 * @param[in] j1 - int - Right boundary index stored by iproc
 * @param[in] grid - const std::vector<int> - Grid to output. Passed by reference
 */
std::vector<int> send_partial_grid(int Ni, int Nj, int j0, int j1, std::vector<int> grid)
{
    int partial_columns = j1 - j0;
    std::vector<int> partial_grid(Ni * partial_columns);

    for (int i = 0; i < Ni; ++i)
    {
        for (int j = 0; j < Nj; ++j)
        {
            if (j < j0 || j >= j1)
            {
                int index = get_1D_index(i, j, Nj);
                grid[index] = 0;
            }
        }
    }
    return grid;
}

/**
 * Outputs grid state of only parts of grid managed by specific iproc. Helpful funcion for debugging
 * @param[in] grid_file - const std::vector<int> & - filestream object that stores opened file
 * @param[in] Ni - int - No. row's of grid
 * @param[in] Nj - int - No. columns's of grid
 * @param[in] iproc - int - Process ID
 * @param[in] nproc - int - No. of processes
 * @param[in] j0 - int - Left boundary index stored by iproc
 * @param[in] j1 - int - Right boundary index stored by iproc
 * @param[in] grid - const std::vector<int> & - Grid to output. Passed by reference
 */
void write_grid(std::ofstream &grid_file, int Ni, int Nj, int j0, int j1, std::vector<int> &grid, int iproc, int nproc)
{
    std::vector<int> full_grid(Ni * Nj, 0);
    std::vector<int> partial_grid = send_partial_grid(Ni, Nj, j0, j1, grid);

    MPI_Reduce(partial_grid.data(), full_grid.data(), Ni * Nj, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (iproc == 0)
    {
        grid_file << Ni << " " << Nj << std::endl;
        for (int i = 0; i < Ni; ++i)
        {
            for (int j = 0; j < Nj; ++j)
            {
                int ind = get_1D_index(i, j, Nj);
                grid_file << i << " " << j << " " << full_grid[ind] << std::endl;
            }
        }
    }
}

/**
 * Distributes columns to tasks based on id.
 * @param[in] grid_file - const std::vector<int> & - filestream object that stores opened file
 * @param[in] Nj - int - No. columns's of grid
 * @param[in] iproc - int - Process ID
 * @param[in] nproc - int - No. of processes
 * @param[in] j0 - int & - Left boundary index to be edited by iproc
 * @param[in] j1 - int & - Right boundary index to be edited by iproc
 */
void distribute_grid(int Nj, int iproc, int nproc, int &j0, int &j1)
{
    int ni = Nj / nproc;
    int rem = Nj % nproc;

    if (iproc < rem)
    {
        j0 = iproc * (Nj + 1);
        j1 = j0 + Nj + 1;
    }
    else
    {
        j0 = iproc * Nj + rem; // Remaining processes get only `Nj` columns
        j1 = j0 + Nj;
    }
}

/**
 * Randomly generates grid of trees based on probability p
 * @param[in] grid_file - const std::vector<int> & - filestream object that stores opened file
 * @param[in] Ni - int - No. row's of grid
 * @param[in] Nj - int - No. columns's of grid
 * @param[in] p - double - Probability that tree is generated
 */
std::vector<int> random_grid(int Ni, int Nj, double p)
{
    std::vector<int> grid(Ni * Nj, 0);
    // randomly fill the initial grid with trees, using probability p
    for (int i = 0; i < Ni; ++i)
    {
        for (int j = 0; j < Nj; ++j)
        {
            // generate a random floating number
            double rn = rng();
            int ind = get_1D_index(i, j, Nj);

            // if the normalised random number is less than probability p, fill the site with a tree
            if (rn < p)
            {
                grid[ind] = 1;
            }
            ind += 1;
        }
    }

    // Set the top row (i = 0) to 2 if it was 1 as top row must start on fire (as requested :o )
    for (int j = 0; j < Nj; ++j)
    {
        int ind = get_1D_index(0, j, Nj);
        if (grid[ind] == 1)
        {
            grid[ind] = 2;
        }
    }
    return grid;
}

/**
 * Processes file input to generate grid from file
 * Robustly fills grid with 0s if user inputs too large values for Ni and Nj
 * Cuts off input if user inputs too little, Simulation should still work
 * @param[in] file_name & - const std::vector<int> & - filestream object that stores opened file
 * @param[in] Ni - int & - No. row's of grid
 * @param[in] Nj - int & - No. columns's of grid
 */
std::vector<int> read_forest(std::string &file_name, int &Ni, int &Nj)
{
    std::cout << "Reading in forest data from " << file_name << std::endl;

    // open the file
    std::ifstream forest_file(file_name);
    if (forest_file.is_open())
    {
        // read in the image data to a single 1D vector
        std::vector<int> forest_data;
        std::vector<int> temp_data;

        // read in state and push to temporary data
        int state;
        while (forest_file >> state)
        {
            temp_data.push_back(state);
        }
        forest_data.resize(Ni * Nj, 0);

        // Loop to write the read data into rectanglular/square forest vector of 0's
        // Loop takes into account if user enters very large Ni and Nj
        // Loop can cut off extra data if user specifies Ni and Nj to be too small
        // The input forest is created by the user... The user must have enough inteligence about their own forest
        for (int i = 0; i < temp_data.size() && i < Ni * Nj; i++)
        {
            forest_data[i] = temp_data[i];
        }

        // Ignite top row
        for (int j = 0; j < Nj; j++)
        {
            int ind = get_1D_index(0, j, Nj);
            if (forest_data[ind] == 1)
            {
                forest_data[ind] = 2;
            }
        }

        // close the file
        forest_file.close();
        return forest_data;
    }
    // Return empty vector for error output
    else
    {
        return {};
    }
}

/**
 * Communicates boundary columns to neighbouring processes
 * Allows fire to spread to columns outside of an iproc's bounds
 * Is only ran when needed
 * @param[in] Ni - int - No. row's of grid
 * @param[in] Nj - int - No. columns's of grid
 * @param[in] j0 - int - Left boundary column held by an iproc
 * @param[in] j1 - int - Right boundary column held by an iproc
 * @param[in] iproc - int - Process ID
 * @param[in] nproc - int - Number of processes
 * @param[in] old_forest - std::vector<int> & - old config of forest before fire spread passed by ref
 * @param[in] new_forest - std::vector<int> & - new config of forest before fire spread passed by ref
 * @param[in] communication_time - double - total communication time accumulated by process passed by ref
 */
void Communicate_Cols(int Ni, int Nj, int j0, int j1, int iproc, int nproc,
                      std::vector<int> &old_forest,
                      std::vector<int> &new_forest,
                      double &communication_time)
{
    double start_communication = MPI_Wtime();
    // Sending and recving by odd and even id's
    for (int oe = 0; oe <= 1; oe++)
    {
        // send
        if (iproc > 0 && iproc % 2 == oe)
        {
            // find the index of the first pixel on row i and column j0
            for (int i = 0; i < Ni; i++)
            {
                int ind = get_1D_index(i, j0, Nj);
                MPI_Send(&old_forest[ind], 1, MPI_INT, iproc - 1, i, MPI_COMM_WORLD); //   send the data, using i as a tag
            }
        }
        // receive
        if (iproc < nproc - 1 && iproc % 2 == (oe + 1) % 2)
        {
            // find the index of the first pixel on row i and column j1
            for (int i = 0; i < Ni; i++)
            {
                int ind = get_1D_index(i, j1, Nj);
                MPI_Recv(&old_forest[ind], 1, MPI_INT, iproc + 1, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //   recieve the data, based on tag i
            }
        }
        MPI_Barrier(MPI_COMM_WORLD); // Barrier to make sure all proc are synchronised before next send - recv phase
    }

    // now send in the other direction, again using even/odd split
    for (int oe = 0; oe <= 1; oe++)
    {

        // send
        if (iproc < nproc - 1 && iproc % 2 == oe)
        {
            // find the index of the first pixel on row i and column j1-1
            for (int i = 0; i < Ni; i++)
            {
                int ind = get_1D_index(i, j1 - 1, Nj);
                MPI_Send(&old_forest[ind], 1, MPI_INT, iproc + 1, i, MPI_COMM_WORLD); //   send the data, using i as a tag
            }
        }
        // receive
        if (iproc > 0 && iproc % 2 == (oe + 1) % 2)
        {
            // find the index of the first pixel on row i and column j0-1
            for (int i = 0; i < Ni; i++)
            {
                int ind = get_1D_index(i, j0 - 1, Nj);
                MPI_Recv(&old_forest[ind], 1, MPI_INT, iproc - 1, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //   recieve the data, based on tag i
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double end_communication = MPI_Wtime();
    communication_time += (end_communication - start_communication); // Barrier to make sure all proc are synchronised before next phase of program
}

/**
 * Updates column boundaries from each iproc reciving a column outside their boundary
 * @param[in] Ni - int - No. row's of grid
 * @param[in] Nj - int - No. columns's of grid
 * @param[in] j0 - int - Left boundary column held by an iproc
 * @param[in] j1 - int - Right boundary column held by an iproc
 * @param[in] iproc - int - Process ID
 * @param[in] nproc - int - Number of processes
 * @param[in] old_forest - std::vector<int> & - old config of forest before fire spread passed by ref
 * @param[in] new_forest - std::vector<int> & - new config of forest before fire spread passed by ref
 * @param[in] calculation_time - double - total calculation time accumulated by process passed by ref
 */
void Update_Col_Boundaries(int Ni, int Nj, int j0, int j1, int iproc, int nproc,
                           std::vector<int> &old_forest,
                           std::vector<int> &new_forest,
                           double &calculation_time)
{
    double start_calculation = MPI_Wtime();
    for (int i = 0; i < Ni; i++)
    {
        if (iproc > 0) // Update left boundary without iproc 0 doing so
        {
            int ind = get_1D_index(i, j0 - 1, Nj);
            if (old_forest[ind] == 2) // if cell on left of boundary is on fire
            {
                int ind_local = get_1D_index(i, j0, Nj);
                if (new_forest[ind_local] == 1) // If cell on boundary is tree
                {
                    new_forest[ind_local] = 2; // uh oh...
                }
            }
        }

        // Repeat the same but for right boundary
        if (iproc < nproc - 1) // Update right boundary
        {
            int ind = get_1D_index(i, j1, Nj);
            if (old_forest[ind] == 2) // if cell on left of boundary is on fire
            {
                int ind_local = get_1D_index(i, j1 - 1, Nj);
                if (new_forest[ind_local] == 1) // If cell on boundary is tree
                {
                    new_forest[ind_local] = 2; // uh oh... (again)
                }
            }
        }
    }
    double end_calculation = MPI_Wtime();

    calculation_time += end_calculation - start_calculation; // Cumulate calculation time for each process
}

/**
 * Sets old forest to new forest config. Equiv. with saying old_forest = new_forest if nproc = 1
 * @param[in] Ni - int - No. row's of grid
 * @param[in] Nj - int - No. columns's of grid
 * @param[in] j0 - int - Left boundary column held by an iproc
 * @param[in] j1 - int - Right boundary column held by an iproc
 * @param[in] old_forest - std::vector<int> & - old config of forest before fire spread passed by ref
 * @param[in] new_forest - std::vector<int> & - new config of forest before fire spread passed by ref
 */
void Reset_Old_Forest(int Ni, int Nj, int j0, int j1,
                      std::vector<int> &old_forest,
                      std::vector<int> &new_forest)
{
    for (int i = 0; i < Ni; i++)
    {
        for (int j = j0; j < j1; j++)
        {
            int index = get_1D_index(i, j, Nj);
            old_forest[index] = new_forest[index];
        }
    }
}

/**
 * The main part of the Forest Fire Model. 
 * For a cell on fire (State 2),
 *      if Von-Neumann neighbouring cells are trees (State 1)
 *          neighboruing cells burn
 *      otherwise do nothing and check next cells
 * @param[in] Ni - int - No. row's of grid
 * @param[in] Nj - int - No. columns's of grid
 * @param[in] j0 - int - Left boundary column held by an iproc
 * @param[in] j1 - int - Right boundary column held by an iproc
 * @param[in] iproc - int - Process ID
 * @param[in] nproc - int - Number of processes
 * @param[in] old_forest - std::vector<int> & - old config of forest before fire spread passed by ref
 * @param[in] new_forest - std::vector<int> & - new config of forest before fire spread passed by ref
 * @param[in] reach_bot - bool & - flag that changes to true if the fire has reached the bottom
 * @return ForestState - Struct that holds whever the fire is still burning or columns need to be exchanged.
 */
ForestState update_forest_state(int Ni, int Nj, int j0, int j1, int iproc, int nproc,
                                std::vector<int> &old_forest,
                                std::vector<int> &new_forest,
                                bool &reach_bot)
{
    // Set local flags to false for ForestState struct
    bool burning = false;
    bool comm_cols = false;
    ForestState state(burning, comm_cols);

    // Run Forest Fire Model Calculations
    for (int i = 0; i < Ni; i++)
    {
        for (int j = j0; j < j1; j++)
        {
            int index = get_1D_index(i, j, Nj);
            if (old_forest[index] == 2) // If tree is burning
            {
                // Store indices for Von Neumann Neighbours
                int abv_index = get_1D_index(i - 1, j, Nj);
                int bel_index = get_1D_index(i + 1, j, Nj);
                int lef_index = get_1D_index(i, j - 1, Nj);
                int rig_index = get_1D_index(i, j + 1, Nj);

                // Check neighbours for burnables (trees)
                if (i > 0 && old_forest[abv_index] == 1) // Above
                {
                    new_forest[abv_index] = 2;
                    burning = true;
                }
                if (i + 1 < Ni && old_forest[bel_index] == 1) // Below
                {
                    new_forest[bel_index] = 2;
                    burning = true;
                    if (i + 1 == Ni - 1) // Check if fire has reached bottom
                    {
                        reach_bot = true;
                    }
                }
                if (j > 0 && old_forest[lef_index] == 1) // Left
                {
                    new_forest[lef_index] = 2;
                    burning = true;
                    if (j == j0 && nproc > 1) // Check if fire is on left of left boundary
                    {
                        comm_cols = true;
                    }
                }
                if (j < Nj && old_forest[rig_index] == 1) // Right
                {
                    new_forest[rig_index] = 2;
                    burning = true;
                    if (j + 1 == j1 && nproc > 1) // Check if fire is on left of left boundary
                    {
                        comm_cols = true;
                    }
                }
                new_forest[index] = 3; // tree has deceased :')
            }
        }
    }
    if (nproc == 1)
    {
        // No MPI communication needed for a single process, simply update forest states as on the tin
        old_forest = new_forest;
        state = {burning, comm_cols};
    }
    else
    {
        // set flags for all processes with Allreduce
        bool global_burning = burning;
        bool global_comm_cols = comm_cols;
        MPI_Allreduce(&burning, &global_burning, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
        MPI_Allreduce(&comm_cols, &global_comm_cols, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
        state = {global_burning, global_comm_cols};
    }

    double end_calculation = MPI_Wtime();

    return state;
}

/**
 * Main (im sure nearly every c++ program has one of those)
 * Takes cmd line arguments
 */
int main(int argc, char **argv)
{
    // Check if there's too little or too many cmd line args
    if (argc < 6 || argc > 7)
    {
        std::cerr << "<tasks> <write (true/false)> <Ni> <Nj> <p> <file name [optional]>\n"; // output to user what we want
        return EXIT_FAILURE; // Not me I hope
    }

    // Initiallise MPI
    MPI_Init(&argc, &argv);

    ////////////////////////////////////////////////////////////
    // read in various input parameters from the command line //
    ////////////////////////////////////////////////////////////

    // Setting the no. of processes/tasks using cmd line arg
    int nproc = atoi(argv[1]);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // Get ID
    int iproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

    std::string write_str = argv[2]; // whether to write the grid to file
    bool write = (write_str == "true"); // Set that to true if it is
    int Ni = atoi(argv[3]);   // the size of the grid in the x-dimension
    int Nj = atoi(argv[4]);   // the size of the grid in the y-dimension
    double p = atof(argv[5]); // probability tree is generated

    // Check if file name is provided
    std::string file_name;
    if (argc > 6)
    {
        file_name = argv[6];
    }
    else
    {
        file_name = "no file"; // No file name provided
    }

    // Define grids
    std::vector<int> grid;
    std::vector<int> new_grid;

    // Define grid config id and char array for later broadcast
    std::string config_id;
    char recieved_id[6];

    // Make process 0 generate and initialise grid and new_grid
    if (iproc == 0)
    {
        if (file_name != "no file") // If theres file
        {
            grid = read_forest(file_name, Ni, Nj); // Attempt to read forest
            if (grid.size() == 0) // If forest read attempt returned nothing
            {
                std::cerr << file_name << " Not found!!!" << std::endl; // report file not found 
                return EXIT_FAILURE;
            }
            new_grid = grid;
        }
        // Root process initializes the grid if no file is input/attempted to be input
        else
        {
            grid = random_grid(Ni, Nj, p); // Randomly generate grid
            new_grid = grid; // Set new grid as old grid
        }
        config_id = generate_config_id(); // Generate grid config ID
        std::sprintf(recieved_id, "%s", config_id.c_str()); // Convert  ID string into char array stored in recieved_id
    }

    // initialise grid and new_grid before proc0 broadcasts
    else
    {
        grid.resize(Ni * Nj, 0);
        new_grid.resize(Ni * Nj, 0);
    }

    double start_broadcast = MPI_Wtime(); // Start broadcast time measurement
    if (nproc > 1)
    {
        // Share the grid size
        MPI_Bcast(&Ni, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&Nj, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // Share grid data
        MPI_Bcast(grid.data(), Ni * Nj, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(new_grid.data(), Ni * Nj, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(recieved_id, 6, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
    double end_boadcast = MPI_Wtime(); // End broadcast time measurement

    // Set config ID to be receved id
    config_id = recieved_id;

    // Start distributing columns to processes
    int j0, j1;
    distribute_grid(Nj, iproc, nproc, j0, j1);

    // Handling file writing (if needed)
    std::ofstream grid_file; // Initialise fstream object

    // write to file if write = true
    if (write)
    {
        if (iproc == 0)
        {
            std::string file_name = "output_data/burning_forest_" + std::to_string(Ni) + "x" + std::to_string(Nj) + "_p" + std::to_string(p) + "_" + config_id + ".dat";
            grid_file.open(file_name, std::ios::trunc); // Create file called the file_name generated
        }
    }

    // Initialise burning state, col communication state and step count
    ForestState state(true, false);
    int nsteps = 0;
    int nsteps_bot = -1;
    bool take_nsteps_bot = false;

    // Synchronise processes before the fun (forest simulation) begins
    MPI_Barrier(MPI_COMM_WORLD);
    if (write)
    {
        write_grid(grid_file, Ni, Nj, j0, j1, grid, iproc, nproc);
    }

    double start_burn = MPI_Wtime(); // Start burn time measurement
    double calculation_time = 0.0;   // Initialise calculation time measurement
    double communication_time = 0.0; // Initialise proc communication time measurement

    // Update forest while it is burning (burning = true)
    while (state.burning == true)
    {
        take_nsteps_bot = false;
        double start_calculation = MPI_Wtime();
        state = update_forest_state(Ni, Nj, j0, j1, iproc, nproc, grid, new_grid, take_nsteps_bot);
        nsteps += 1; // Increment step count by 1
        if (take_nsteps_bot && nsteps_bot < nsteps)
        {
            nsteps_bot = nsteps;
        }
        if (state.comm_cols && nproc > 1)
        {
            Communicate_Cols(Ni, Nj, j0, j1, iproc, nproc, grid, new_grid, communication_time);
            Update_Col_Boundaries(Ni, Nj, j0, j1, iproc, nproc, grid, new_grid, calculation_time);
            Reset_Old_Forest(Ni, Nj, j0, j1, grid, new_grid);
        }
        else if (state.comm_cols == false && nproc > 1)
        {
            Reset_Old_Forest(Ni, Nj, j0, j1, grid, new_grid);
        }
        else if (nproc == 1)
        {
            grid = new_grid;
        }
        double end_calculation = MPI_Wtime();
        calculation_time += end_calculation - start_calculation;
        if (write)
        {
            write_grid(grid_file, Ni, Nj, j0, j1, grid, iproc, nproc);
        }
    }
    double end_burn = MPI_Wtime();

    // Synchronise processes before final claculations
    MPI_Barrier(MPI_COMM_WORLD);

    // Calculate all times taken to carry out various steps
    double broadcast_time = end_boadcast - start_broadcast;
    double burning_time = end_burn - start_burn;
    double total_MPI_time = end_burn - start_broadcast;

    int max_n_steps_bot;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&nsteps_bot, &max_n_steps_bot, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max_n_steps_bot, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Output
    std::cout << config_id << "\t" << iproc << "\t" << nproc << "\t" << Ni << "\t" << Nj << "\t" << p << "\t" << total_MPI_time << "\t" << broadcast_time << "\t" << calculation_time << "\t" << communication_time << "\t" << burning_time << "\t" << nsteps << "\t" << max_n_steps_bot << std::endl;

    // Finalise MPI (but in american spelling)
    MPI_Finalize();

    return EXIT_SUCCESS; // eqv. MPI_Gain(sanity back) but thats not what has happened
}