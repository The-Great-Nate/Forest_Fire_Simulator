#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cstdlib>

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

bool update_forest_serial(int Ni, int Nj, std::vector<std::vector<int>> &grid)
{
    bool burning = false;
    std::vector<std::vector<int>> new_grid = grid;
    for (int i = 0; i < Ni; i++)
    {
        for (int j = 0; j < Nj; j++)
        {
            if (grid[i][j] == 2) // If tree is burning
            {
                if (i > 0 && grid[i - 1][j] == 1) // Above
                {
                    new_grid[i - 1][j] = 2;
                    burning = true;
                }
                if (i + 1 < Ni && grid[i + 1][j] == 1) // Below
                {
                    new_grid[i + 1][j] = 2;
                    burning = true;
                }
                if (j > 0 && grid[i][j - 1] == 1) // Left
                {
                    new_grid[i][j - 1] = 2;
                    burning = true;
                }
                if (j + 1 < Nj && grid[i][j + 1] == 1) // Right
                {
                    new_grid[i][j + 1] = 2;
                    burning = true;
                }
                new_grid[i][j] = 3; // Tree evolves into Ash... or dies, oops
            }
        }
    }

    grid = new_grid; // Copy updated grid back
    return burning;  // Return true if fire is still burning
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
    if (argc < 5)
    {
        std::cerr << " <write (true/false)> <Ni> <Nj> <p>\n";
        return EXIT_FAILURE;
    }
    // read in various input parameters from the command line
    std::string write_str = argv[1]; // whether to write the grid to file
    bool write = (write_str == "true");
    int Ni = atoi(argv[2]);   // the size of the grid in the x-dimension
    int Nj = atoi(argv[3]);   // the size of the grid in the y-dimension
    double p = atof(argv[4]); // probability tree is generated
    std::vector<std::vector<int>> grid = random_grid(Ni, Nj, p);
    //display_grid(grid);

    std::ofstream grid_file("burning_forest.dat");
    if (write)
    {
        write_grid(grid_file, Ni, Nj, grid);
    }

    bool burning = true;
    int nsteps = 0;
    while (burning)
    {
        burning = update_forest_serial(Ni, Nj, grid);
        if (burning)
        {
            if (write)
            {
                write_grid(grid_file, Ni, Nj, grid);
            }
        }
    }
}   