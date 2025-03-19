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
std::vector<std::vector<int>> random_grid(int Ni, int Nj, double p)
{
    std::vector<std::vector<int>> grid(Ni, std::vector<int>(Nj, 0));
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
                grid[i][j] = 1;
            }
        }
    }

    // Set the top row (i = 0) to 2 if it was 1
    for (int j = 0; j < Nj; ++j)
    {
        if (grid[0][j] == 1)
        {
            grid[0][j] = 2;
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

bool update_forest_serial(int Ni, int Nj, std::vector<std::vector<int>> &old_forest)
{
    bool burning = false;
    std::vector<std::vector<int>> new_forest = old_forest;
    for (int i = 0; i < Ni; i++)
    {
        for (int j = 0; j < Nj; j++)
        {
            if (old_forest[i][j] == 2) // If tree is burning
            {
                if (i > 0 && old_forest[i - 1][j] == 1) // Above
                {
                    new_forest[i - 1][j] = 2;
                    burning = true;
                }
                if (i + 1 < Ni && old_forest[i + 1][j] == 1) // Below
                {
                    new_forest[i + 1][j] = 2;
                    burning = true;
                }
                if (j > 0 && old_forest[i][j - 1] == 1) // Left
                {
                    new_forest[i][j - 1] = 2;
                    burning = true;
                }
                if (j + 1 < Nj && old_forest[i][j + 1] == 1) // Right
                {
                    new_forest[i][j + 1] = 2;
                    burning = true;
                }
                new_forest[i][j] = 3; // Tree evolves into Ash... or dies, oops
            }
        }
    }

    old_forest = new_forest; // Copy updated grid back
    return burning;          // Return true if fire is still burning
}

bool update_forest_MPI(int Ni, int Nj, int j0, int j1, int iproc, int nproc, std::vector<std::vector<int>> &old_forest, std::vector<std::vector<int>> &new_forest, double &calculation_time, double &communication_time)
{
    bool burning = false;
    double start_calculation = MPI_Wtime();

    for (int i = 0; i < Ni; i++)
    {
        for (int j = j0; j < j1; j++)
        {
            if (old_forest[i][j] == 2) // If tree is burning
            {
                if (i > 0 && old_forest[i - 1][j] == 1) // Above
                {
                    new_forest[i - 1][j] = 2;
                    burning = true;
                }
                if (i + 1 < Ni && old_forest[i + 1][j] == 1) // Below
                {
                    new_forest[i + 1][j] = 2;
                    burning = true;
                }
                if (j > 0 && old_forest[i][j - 1] == 1) // Left
                {
                    new_forest[i][j - 1] = 2;
                    burning = true;
                }
                if (j + 1 < Nj && old_forest[i][j + 1] == 1) // Right
                {
                    new_forest[i][j + 1] = 2;
                    burning = true;
                }
                new_forest[i][j] = 3; // Tree evolves into Ash... or dies, oops
            }
        }
    }

    double end_calculation = MPI_Wtime();

    double start_communication = MPI_Wtime();

    if (nproc > 1)
    {
        for (int oe = 0; oe <= 1; oe++) // Handle sending of even columns, for odd columns to recive
        {
            if (iproc > 0 && iproc % 2 == oe) // send left col to left neighbour
            {
                for (int i = i0; i < i1; i++)
                {
                    MPI_Send(&new_forest[i][j0], 1, MPI_INT, iproc - 1, j0, MPI_COMM_WORLD);
                }
            }
            if (iproc < nproc - 1 && iproc % 2 == (oe + 1) % 2) // recieve left col from right neighbour -> updating right neighbour right col
            {
                for (int i = i0; i < i1; i++)
                {
                    MPI_Recv(&old_forest[i][j1], 1, MPI_INT, iproc + 1, j1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }

        for (int oe = 0; oe <= 1; oe++)
        {
            if (iproc < nproc - 1 && iproc % 2 == oe) // send right col to right neighbour
            {
                for (int i = i0; i < i1; i++)
                {
                    MPI_Send(&new_forest[i][j1 - 1], 1, MPI_INT, iproc + 1, j1-1, MPI_COMM_WORLD);
                }
            }
            if (iproc > 0 && iproc % 2 == (oe + 1) % 2) // recieve right column from left neighbour -> updating right neighbour left col
            {
                for (int i = i0; i < i1; i++)
                {
                    MPI_Recv(&old_forest[i][j0 - 1], 1, MPI_INT, iproc - 1, j0-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        
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

    double local_computation_time = end_calculation - start_calculation;
    double local_communication_time = end_communication - start_communication;

    std::cout << "Process\t" << iproc << "\tCalculation Time\t" << local_computation_time << "\tCommunication Time\t" << local_communication_time << std::endl;

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
    std::vector<std::vector<int>> grid;

    if (iproc == 0)
    {
        // Root process initializes the grid
        grid = random_grid(Ni, Nj, p);
    }

    if (nproc > 1)
    {
        // first share the image size
        MPI_Bcast(&Ni, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&Nj, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // flattening the grid for broadcasting
        std::vector<int> flat_grid(Ni * Nj);
        if (iproc == 0)
        {
            int f_ind = 0;
            for (int i = 0; i < Ni; ++i)
            {
                for (int j = 0; j < Nj; ++j)
                {
                    flat_grid[f_ind] = grid[i][j];
                    f_ind += 1;
                }
            }
        }
        MPI_Bcast(flat_grid.data(), Ni * Nj, MPI_INT, 0, MPI_COMM_WORLD);
    }
    
    // reconstruct grid
    if (iproc != 0)
    {
        int f_ind = 0;
        for (int i = 0; i < Ni; ++i)
        {
            for (int j = 0; j < Nj; ++j)
            {
                grid[i][j] = flat_grid[f_ind];
                f_ind += 1
            }
        }
    }

    std::vector < int > new_grid(Ni*Nj, 0);

    int j0, j1;
    distribute_grid(Nj, iproc, nproc, j0, j1);

    std::ofstream grid_file("burning_forest.dat");
    if (write)
    {
        write_grid(grid_file, Ni, Nj, grid);
    }

    bool burning = true;
    int nsteps = 0;
    double start_time = MPI_Wtime();
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
    double end_time = MPI_Wtime();
    // Finalise MPI (but in american spelling)
    MPI_Finalize();

    std::cout << "Process_vs_Steps\t" << nproc << "\t" << nsteps << "\t" << end_time - start_time << std::endl;
}