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

                // if the normalised random number is less than our probability p, fill the site with a tree
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
                }
            }
        }

        // Update old forest after new trees meet their doom...
        for (int i = 0; i < Ni; i++)
        {
            for (int j = j0; j < j1; j++)
            {
                int index = get_1D_index(i, j, Nj);
                if (old_forest[index] == 2)
                {
                    new_forest[index] = 3; // Tree evolves into Ash... or dies, oops
                }
            }
        }
        double start_communication = MPI_Wtime();

        if (nproc > 1)
        {
            for (int oe = 0; oe <= 1; oe++)
            {

                // send
                if (iproc > 0 && iproc % 2 == oe)
                {
                    // find the index of the first pixel on row i0 - this is only guaranteed to work if j0=0, so we could put 0 instead (see comment on BB)
                    for (int i = 0; i < Ni; i++)
                    {
                        int ind = get_1D_index(i, j0, Nj);
                        std::cout << "Sending1 row " << i << " (ind=" << ind << ") from task " << iproc << " to task " << iproc - 1 << std::endl;
                        MPI_Send(&new_forest[ind], Nj, MPI_INT, iproc - 1, i, MPI_COMM_WORLD);
                    }
                }
                MPI_Barrier(MPI_COMM_WORLD);
                // receive
                if (iproc < nproc - 1 && iproc % 2 == (oe + 1) % 2)
                {
                    // find the index of the first pixel on row i1 - this is only guaranteed to work if j0=0, so we could put 0 instead (see comment on BB)
                    for (int i = 0; i < Ni; i++)
                    {
                        int ind = get_1D_index(i, j1, Nj);
                        // receive the data, where each row has Nj * 3 elements, using i1 as the tag
                        std::cout << "Receiving1 row " << i << " (ind=" << ind << ") on task " << iproc << " from task " << iproc + 1 << std::endl;
                        MPI_Recv(&old_forest[ind], Nj, MPI_INT, iproc + 1, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }

            // now send in the other direction, again using our even/odd split
            for (int oe = 0; oe <= 1; oe++)
            {

                // send
                if (iproc < nproc - 1 && iproc % 2 == oe)
                {
                    // find the index of the first pixel on row i1-1 - this is only guaranteed to work if j0=0, so we could put 0 instead (see comment on BB)
                    for (int i = 0; i < Ni; i++)
                    {
                        int ind = get_1D_index(i, j1, Nj);
                        // send the data, where each row has Nj * 3 elements, using i1-1 as the tag
                        std::cout << "Sending2 row " << i << " (ind=" << ind << ") from task " << iproc << " to task " << iproc + 1 << std::endl;
                        MPI_Send(&new_forest[ind], Nj, MPI_INT, iproc + 1, i, MPI_COMM_WORLD);
                    }
                }
                MPI_Barrier(MPI_COMM_WORLD);
                // receive
                if (iproc > 0 && iproc % 2 == (oe + 1) % 2)
                {
                    // find the index of the first pixel on row i0-1 - this is only guaranteed to work if j0=0, so we could put 0 instead (see comment on BB)
                    for (int i = 0; i < Ni; i++)
                    {
                        int ind = get_1D_index(i, j0, Nj);
                        // receive the data, where each row has Nj * 3 elements, using i1 as the tag
                        std::cout << "Receiving2 row " << i << " (ind=" << ind << ") on task " << iproc << " from task " << iproc - 1 << std::endl;
                        MPI_Recv(&old_forest[ind], Nj, MPI_INT, iproc - 1, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
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
            int local_size = Ni * (j1 - j0);
            MPI_Allreduce(&local_burning, &burning, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
            MPI_Allgatherv(new_forest.data(), local_size, MPI_INT, old_forest.data(), recvcounts.data(), displs.data(), MPI_INT, MPI_COMM_WORLD);
        }
        else
        {
            // No MPI communication needed for a single process
            old_forest = new_forest;
        }

        double end_communication = MPI_Wtime();

        calculation_time += start_communication - start_calculation;
        communication_time += end_communication - start_communication;

        return burning;
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

        // Define grids
        std::vector<int> grid;
        std::vector<int> new_grid;

        // Make process 0 generate and initialise grid and new_grid
        if (iproc == 0)
        {
            // Root process initializes the grid
            grid = random_grid(Ni, Nj, p);
            new_grid = grid;
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
        std::cout << iproc << " " << j0 << " " << j1 << std::endl;
        //std::cout << iproc << " " << j0 << " " << j1 << std::endl;
        // Handling file writing (if needed)
        std::ofstream grid_file("burning_forest.dat"); // Initialise fstream object
        // write to file if write = true
        if (write)
        {
            write_grid(grid_file, Ni, Nj, grid);
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
            burning = update_forest_MPI(Ni, Nj, j0, j1, iproc, nproc, grid, new_grid, calculation_time, communication_time, recvcounts, displs);
            // display_grid(Ni, Nj, grid);
            std::cout << nsteps << std::endl;
            // Carry out if write = true
            if (write)
            {
                write_grid(grid_file, Ni, Nj, grid);
            }
            if (nsteps > 3)
            {
                break;
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