#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cstdlib>
#include <mpi.h>

// Forest Fires are not good

double rng()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    double random_float = dist(gen);
    return random_float;
}

int get_1D_index(int i, int j, int Nj)
{
    return (i * Nj) + j;
}

void display_grid(int Ni, int Nj, const std::vector<int> &grid)
{
    for (int i = 0; i < Ni; ++i)
    {
        for (int j = 0; j < Nj; ++j)
        {
            int ind = get_1D_index(i, j, Nj);
            std::cout << grid[ind] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

void display_grid_proc(int Ni, int Nj, int j0, int j1, const std::vector<int> &grid)
{
    for (int i = 0; i < Ni; ++i)
    {
        for (int j = j0; j < j1; ++j)
        {
            int ind = get_1D_index(i, j, Nj);
            std::cout << grid[ind] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

// function to append the current grid of to an already open file
void write_grid(std::ofstream &grid_file, int Ni, int Nj, std::vector<int> grid)
{
    // write the size of the image
    grid_file << Ni << " " << Nj << std::endl;
    int ind = 0;
    // now write out the grid
    for (int i = 0; i < Ni; ++i)
    {
        for (int j = 0; j < Nj; ++j)
        {
            ind = get_1D_index(i, j, Nj);
            grid_file << i << " " << j << " " << grid[ind] << std::endl;
        }
    }
}

void distribute_grid(int Nj, int iproc, int nproc, int &i0, int &i1)
{
    int ni = Nj / nproc;
    int rem = Nj % nproc;

    if (iproc < rem)
    {
        i0 = iproc * (ni + 1);
        i1 = i0 + ni + 1;
    }
    else
    {
        i0 = iproc * ni + rem; // Remaining processes get only `Nj` columns
        i1 = i0 + ni;
    }
}

// function which initialises a random grid
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

    // Set the top row (i = 0) to 2 if it was 1
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

std::vector<int> read_forest(std::string filename, int &Ni, int &Nj)
{

    std::cout << "Reading in forest data from " << filename << std::endl;

    // open the file
    std::ifstream forest_file(filename);

    // read in the image data to a single 1D vector
    std::vector<int> forest_data;
    int state;
    while (forest_file >> state)
    {
        forest_data.push_back(state);
    }

    std::cout << "\n";

    for (int j = 0; j < Nj; j++)
    {
        int ind = get_1D_index(0, j, Nj);
        if (forest_data[ind] == 1)
        {
            forest_data[ind] = 2;
        }
    }
    std::cout << "\n";
    // close the file
    forest_file.close();
    return forest_data;
}

bool update_forest_MPI(int Ni, int Nj, int j0, int j1, int iproc, int nproc,
                       std::vector<int> &old_forest,
                       std::vector<int> &new_forest,
                       double &calculation_time,
                       double &communication_time,
                       std::vector<int> &recvcounts,
                       std::vector<int> &displs)
{
    bool burning = false;
    double start_calculation = MPI_Wtime();
    if (iproc == 0)
    {
        // display_grid_proc(Ni, Nj, j0, j1, old_forest);
        display_grid(Ni, Nj, old_forest);
    }
    for (int i = 0; i < Ni; i++)
    {
        for (int j = j0; j < j1; j++)
        {
            int index = get_1D_index(i, j, Nj);
            if (old_forest[index] == 2) // If tree is burning
            {
                int abv_index = get_1D_index(i - 1, j, Nj);
                int bel_index = get_1D_index(i + 1, j, Nj);
                int lef_index = get_1D_index(i, j - 1, Nj);
                int rig_index = get_1D_index(i, j + 1, Nj);
                if (iproc == 0 && i == 2)
                {
                    std::cout << "THE THING IN THE MIDDLE " << old_forest[index] << std::endl;
                    std::cout << "THE THING TO THE RIGHT " << old_forest[rig_index] << std::endl;
                }
                if (i > 0 && old_forest[abv_index] == 1) // Above
                {
                    new_forest[abv_index] = 2;
                    burning = true;
                }
                if (i + 1 < Ni && old_forest[bel_index] == 1) // Below
                {
                    new_forest[bel_index] = 2;
                    burning = true;
                }
                if (j > 0 && old_forest[lef_index] == 1) // Left
                {
                    new_forest[lef_index] = 2;
                    burning = true;
                }
                if (j + 1 < Nj && old_forest[rig_index] == 1) // Right
                {
                    std::cout << "I SHOULD GO RIGHT MABE EHEHEHEHE" << std::endl;
                    new_forest[rig_index] = 2;
                    burning = true;
                }
                new_forest[index] = 3;
            }
        }
    }

    // Update old forest after new trees meet their doom...

    double start_communication = MPI_Wtime();

    if (nproc > 1)
    {
        for (int oe = 0; oe <= 1; oe++)
        {

            // send
            if (iproc > 0 && iproc % 2 == oe)
            {
                // find the index of the first pixel on row i and column j0
                for (int i = 0; i < Ni; i++)
                {
                    int ind = get_1D_index(i, j0, Nj);
                    // std::cout << "Sending1 row " << i << " col " << j0 << " (ind=" << ind << ") from task " << iproc << " to task " << iproc - 1 << std::endl;
                    // send the data, using i as a tag
                    MPI_Send(&new_forest[ind], 1, MPI_INT, iproc - 1, i, MPI_COMM_WORLD);
                }
            }
            // receive
            if (iproc < nproc - 1 && iproc % 2 == (oe + 1) % 2)
            {
                // find the index of the first pixel on row i and column j1
                for (int i = 0; i < Ni; i++)
                {
                    int ind = get_1D_index(i, j1, Nj);
                    // receive the data, where each row has Nj * 3 elements, using i1 as the tag
                    // std::cout << "Receiving1 row " << i << " col " << j1 << " (ind=" << ind << ") on task " << iproc << " from task " << iproc + 1 << std::endl;
                    MPI_Recv(&old_forest[ind], 1, MPI_INT, iproc + 1, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }

        // now send in the other direction, again using our even/odd split
        for (int oe = 0; oe <= 1; oe++)
        {

            // send
            if (iproc < nproc - 1 && iproc % 2 == oe)
            {
                // find the index of the first pixel on row i and column j1-1
                for (int i = 0; i < Ni; i++)
                {
                    int ind = get_1D_index(i, j1 - 1, Nj);
                    // send the data, using i as a tag
                     std::cout << "sending val " << new_forest[ind] <<" Sending2 row " << i << " col " << j1 - 1 << " (ind=" << ind << ") from task " << iproc << " to task " << iproc + 1 << std::endl;
                     //std::cout << "it's send: " << new_forest[ind] << std::endl;
                    MPI_Send(&new_forest[ind], 1, MPI_INT, iproc + 1, i, MPI_COMM_WORLD);
                }
            }
            // receive
            if (iproc > 0 && iproc % 2 == (oe + 1) % 2)
            {
                // find the index of the first pixel on row i and column j0-1
                for (int i = 0; i < Ni; i++)
                {
                    int ind = get_1D_index(i, j0 - 1, Nj);
                     //std::cout << "it's before: " << new_forest[ind] << std::endl;
                    MPI_Recv(&old_forest[ind], 1, MPI_INT, iproc - 1, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                     //std::cout << "it's: Now " << old_forest[ind] << std::endl;
                     std::cout << "recieing val " << old_forest[ind] << " Receiving2 row " << i << " col " << j0 - 1 << " (ind=" << ind << ") on task " << iproc << " from task " << iproc - 1 << std::endl;
                }
            }
        } // Maybe try updating the grid again after the send recieves...
        MPI_Barrier(MPI_COMM_WORLD);

        for (int i = 0; i < Ni; i++)
        {
            if (iproc > 0) // Update left boundary
            {
                int ind = get_1D_index(i, j0 - 1, Nj);
                if (old_forest[ind] == 2)
                {
                    int ind_local = get_1D_index(i, j0, Nj);
                    if (new_forest[ind_local] == 1)
                    {
                        new_forest[ind_local] = 2;
                        burning = true;
                    }
                }
            }

            if (iproc < nproc - 1) // Update right boundary
            {
                int ind = get_1D_index(i, j1, Nj);
                if (old_forest[ind] == 2) // If the received cell is burning
                {
                    int ind_local = get_1D_index(i, j1 - 1, Nj);
                    if (new_forest[ind_local] == 1)
                    {
                        new_forest[ind_local] = 2;
                        burning = true;
                    }
                }
            }
        }

        // copy the items belonging to a task from new forest to old forest
        for (int i = 0; i < Ni; i++)
        {
            for (int j = j0; j < j1; j++)
            {
                int ind = get_1D_index(i, j, Nj);
                old_forest[ind] = new_forest[ind];
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        bool local_burning = burning; // Store local value
        // int local_size = Ni * (j1 - j0);
        MPI_Allreduce(&local_burning, &burning, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
        // MPI_Allgatherv(new_forest.data(), local_size, MPI_INT, old_forest.data(), recvcounts.data(), displs.data(), MPI_INT, MPI_COMM_WORLD);
    }
    if (nproc == 1)
    {
        // No MPI communication needed for a single process
        old_forest = new_forest;
    }

    double end_communication = MPI_Wtime();

    calculation_time += start_communication - start_calculation;
    communication_time += end_communication - start_communication;

    return burning;
}

int main(int argc, char **argv)
{
    if (argc < 6 || argc > 7)
    {
        std::cerr << "<tasks> <write (true/false)> <Ni> <Nj> <p> <file name [optional]>\n";
        return EXIT_FAILURE;
    }
    // Initiallise MPI
    MPI_Init(&argc, &argv);

    // read in various input parameters from the command line

    // Setting the no. of processes/tasks using cmd line arg
    int nproc = atoi(argv[1]);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // Get ID
    int iproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

    std::string write_str = argv[2]; // whether to write the grid to file
    bool write = (write_str == "true");
    int Ni = atoi(argv[3]);   // the size of the grid in the x-dimension
    int Nj = atoi(argv[4]);   // the size of the grid in the y-dimension
    double p = atof(argv[5]); // probability tree is generated

    // Check if name is provided
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

    // Make process 0 generate and initialise grid and new_grid
    if (iproc == 0)
    {
        if (file_name != "no file")
        {
            grid = read_forest(file_name, Ni, Nj);
            new_grid = grid;
        }
        // Root process initializes the grid
        else
        {
            grid = random_grid(Ni, Nj, p);
            new_grid = grid;
        }
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
    }
    double end_boadcast = MPI_Wtime(); // End broadcast time measurement

    // Start distributing columns to processes
    int j0, j1;
    distribute_grid(Nj, iproc, nproc, j0, j1);

    // Initialise values needed for MPI_Allgatherv
    std::vector<int> recvcounts(nproc); // Amount of elements each process contributes
    std::vector<int> displs(nproc);     // Position of where data starts in recive buffer
    for (p = 0; p < nproc; p++)
    {
        int temp_j0, temp_j1;
        distribute_grid(Nj, p, nproc, temp_j0, temp_j1);
        recvcounts[p] = Ni * (temp_j1 - temp_j0);
        displs[p] = temp_j0;
    }
    // std::cout << iproc << " " << j0 << " " << j1 << std::endl;
    //  std::cout << iproc << " " << j0 << " " << j1 << std::endl;
    //   Handling file writing (if needed)

    std::ofstream grid_file; // Initialise fstream object

    // write to file if write = true
    if (write)
    {
        if (iproc == 0)
        {
            grid_file.open("burning_forest.dat", std::ios::trunc);
            write_grid(grid_file, Ni, Nj, grid);
        }
    }

    // Initialise burning state and step count
    bool burning = true;
    int nsteps = 0;

    // Synchronise processes before the fun (forest simulation) begins
    MPI_Barrier(MPI_COMM_WORLD);

    double start_burn = MPI_Wtime(); // Start burn time measurement
    double calculation_time = 0.0;   // Initialise calculation time measurement
    double communication_time = 0.0; // Initialise proc communication time measurement

    // Update forest while it is burning (burning = true)
    while (burning)
    {
        // std::cout << nsteps << std::endl;
        if (iproc == 0)
        {
            std::cout << "-------------------- " << nsteps << " ---------------------" << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // display_grid(Ni, Nj, grid);
        burning = update_forest_MPI(Ni, Nj, j0, j1, iproc, nproc, grid, new_grid, calculation_time, communication_time, recvcounts, displs);
        // Carry out if write = true
        if (write)
        {
            if (iproc == 0)
            {
                write_grid(grid_file, Ni, Nj, grid);
            }
        }
        nsteps += 1; // Increment step count by 1
    }
    double end_burn = MPI_Wtime();

    // Synchronise processes before final claculations
    MPI_Barrier(MPI_COMM_WORLD);

    // Calculate all times taken to carry out various steps
    double broadcast_time = end_boadcast - start_broadcast;
    double burning_time = end_burn - start_burn;
    double total_MPI_time = end_burn - start_broadcast;

    // Output via 1 process only.
    if (iproc == 0)
    {
        std::cout << iproc << "\t" << nproc << "\t" << total_MPI_time << "\t" << broadcast_time << "\t" << calculation_time << "\t" << communication_time << "\t" << burning_time << "\t" << nsteps << std::endl;
    }

    // Finalise MPI (but in american spelling)
    MPI_Finalize();

    return EXIT_SUCCESS; // eqv. MPI_Gain(sanity back) but thats not what 99.95% of the world has
}