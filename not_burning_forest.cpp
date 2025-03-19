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

// function to append the current grid of to an already open file
void write_grid(std::ofstream &grid_file, int Ni, int Nj, std::vector<std::vector<int>> grid)
{
    // write the size of the image
    grid_file << Ni << " " << Nj << std::endl;

    // now write out the grid
    for (int i = 0; i < Ni; ++i)
    {
        for (int j = 0; j < Nj; ++j)
        {
            grid_file << i << " " << j << " " << grid[i][j] << std::endl;
        }
    }
}

// function which initialises a random grid
std::vector<int> random_grid(int Ni, int Nj, double p)
{
    std::vector<int> grid(Ni * Nj, 0);
    int ind = 0
        // randomly fill the initial grid with trees, using probability p
        for (int i = 0; i < Ni; ++i)
    {
        for (int j = 0; j < Nj; ++j)
        {
            // generate a random floating number
            double rn = rng();

            // if the normalised random number is less than our probability p, fill the site with a tree
            if (rn < p)
            {
                grid[ind] = 1;
            }
            ind += 1
        }
    }

    // Set the top row (i = 0) to 2 if it was 1
    for (int j = 0; j < Nj; ++j)
    {
        ind = get_1d_ind(0, j, Nj);
        if (grid[ind] == 1)
        {
            grid[ind] = 2;
        }
    }
    return grid;
}

void distribute_grid(int Nj, int iproc, int nproc, int &j0, int &j1)
{

    j0 = 0;
    j1 = Nj;
    if (nproc > 1)
    {
        int nj = Nj / nproc;
        j0 = iproc * nj;
        j1 = j0 + nj;
        // make sure we take care of the fact that the grid might not easily divide into slices
        if (iproc == nproc - 1)
        {
            j1 = Nj;
        }
    }
}

int get_1d_ind(int i, int j, int Nj)
{
    return (j * Nj) + i
}

bool update_forest_MPI(int Ni, int Nj, int j0, int j1, int iproc, int nproc, std::vector<int> &old_forest, std::vector<int> &new_forest, double &calculation_time, double &communication_time)
{
    bool burning = false;
    double start_calculation = MPI_Wtime();
    for (int i = 0; i < Ni; i++)
    {
        for (int j = j0; j < j1; j++)
        {
            int index = get_1d_ind(i, j, Nj);
            if (old_forest[index] == 2) // If tree is burning
            {
                int abv_index = get_1d_ind(i - 1, j, Nj);
                int bel_index = get_1d_ind(i + 1, j, Nj);
                int lef_index = get_1d_ind(i, j - 1, Nj);
                int rig_index = get_1d_ind(i - 1, j + 1, Nj);
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
                    new_forest[rig_index] = 2;
                    burning = true;
                }
                new_forest[index] = 3; // Tree evolves into Ash... or dies, oops
            }
        }
    }

    double start_communication = MPI_Wtime();

    if (nproc > 1)
    {
        for (int oe = 0; oe <= 1; oe++) // Handle sending of even columns, for odd columns to recive
        {
            if (iproc > 0 && iproc % 2 == oe) // send left col to left neighbour
            {
                int ind = get_1d_ind(0, j0, Nj)
                    MPI_Send(&new_forest[ind], Nj, MPI_INT, iproc - 1, j0, MPI_COMM_WORLD);
            }
            if (iproc < nproc - 1 && iproc % 2 == (oe + 1) % 2) // recieve left col from right neighbour -> updating right neighbour right col
            {
                int ind = get_1d_ind(0, j1, Nj)
                    MPI_Recv(&old_forest[ind], Nj, MPI_INT, iproc + 1, j1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        for (int oe = 0; oe <= 1; oe++)
        {
            if (iproc < nproc - 1 && iproc % 2 == oe) // send right col to right neighbour
            {
                int ind = get_1D_index(0, j1 - 1, Nj);
                MPI_Send(&new_forest[ind], Nj, MPI_INT, iproc + 1, j1 - 1, MPI_COMM_WORLD);
            }
            if (iproc > 0 && iproc % 2 == (oe + 1) % 2) // recieve right column from left neighbour -> updating right neighbour left col
            {
                int ind = get_1D_index(0, j0 - 1, Nj);
                MPI_Recv(&old_forest[ind], Nj, MPI_INT, iproc - 1, j0 - 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        int local_size = Nj * (j1 - j0);
        std::vector<int> temp_forest(Ni * Nj); // temp grid to store reassembled forest
        int start_index = get_1d_ind(0, j0, Nj);
        // Assemble and gather columns from all processes and write to central forest
        MPI_Allgatherv(&new_forest[start_index], local_size, MPI_INT, temp_forest.data(), recvcounts.data(), displs.data(), MPI_INT, MPI_COMM_WORLD);
        old_forest = temp_forest;
    }
    else
    {
        // No MPI communication needed for a single process
        old_forest = new_forest;
    }

    double end_communication = MPI_Wtime();

    bool global_burning;
    MPI_Reduce(&burning, &global_burning, 1, MPI_C_BOOL, MPI_LOR, 0, MPI_COMM_WORLD); // Must check if at least 1 process is still burning
    if (iproc == 0)
    {
        burning = global_burning;
    }

    calculation_time += start_communication - start_calculation;
    communication_time += end_communication - start_communication;


    return burning; // Return true if fire is still burning
}

void display_grid(const std::vector<std::vector<int>> &grid)
{
    int Ni = grid.size();
    int Nj = grid[0].size();
    for (int i = 0; i < Ni; ++i)
    {
        for (int j = 0; j < Nj; ++j)
        {
            std::cout << grid[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

int main(int argc, char **argv)
{
    if (argc < 6)
    {
        std::cerr << "<tasks> <write (true/false)> <Ni> <Nj> <p>\n";
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
    std::vector<int> grid;

    if (iproc == 0)
    {
        // Root process initializes the grid
        grid = random_grid(Ni, Nj, p);
    }

    double start_broadcast = MPI_Wtime();
    if (nproc > 1)
    {
        // first share the grid size
        MPI_Bcast(&Ni, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&Nj, 1, MPI_INT, 0, MPI_COMM_WORLD);
        // second share the actual grid
        for (int i = 0; i < Ni * Nj; i++)
        {
            grid.push_back(0);
        }
        MPI_Bcast(grid.data(), Ni * Nj, MPI_INT, 0, MPI_COMM_WORLD);
    }
    double end_boadcast = MPI_Wtime();
    

    std::vector<int> new_grid(Ni * Nj, 0);

    std::vector<int> recvcounts(nproc);
    std::vector<int> displs(nproc);

    int j0, j1;
    distribute_grid(Nj, iproc, nproc, j0, j1);
    for (int p = 0; p < nproc; p++)
    {
        recvcounts[p] = Ni * (j1 - j0); // Find no. of elements each process stores
        displs[p] = j0;
    }

    std::ofstream grid_file("burning_forest.dat");
    if (write)
    {
        write_grid(grid_file, Ni, Nj, grid);
    }

    bool burning = true;
    int nsteps = 0;
    double start_burn = MPI_Wtime();
    double calculation_time = 0.0, communication_time = 0.0;
    while (burning)
    {
        burning = update_forest_MPI(Ni, Nj, j0, j1, iproc, nproc, grid, new_grid, calculation_time, communication_time);
        if (burning)
        {
            if (write)
            {
                write_grid(grid_file, Ni, Nj, grid);
            }
        }
        nsteps += 1
    }
    double end_burn = MPI_Wtime();

    double broadcast_time = end_boadcast - start_broadcast;
    double burning_time = end_burn - start_burn;
    double total_MPI_time = end_time - start_broadcast;
    // Finalise MPI (but in american spelling)
    MPI_Finalize();
    

    std::cout << nproc << "\t" << total_MPI_time << "\t" << broadcast_time << "\t" << calculation_time << "\t" << communication_time << "\t" << burning_time << "\t" << nsteps << std::endl;
}