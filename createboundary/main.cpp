#include <iostream>
#include <fstream>
#include <sstream>


int main (int argc, char* argv[])
{
    if (argc != 6)
    {
        std::cerr << "ERROR: cmdl arguments not matching!" << std::endl;
        std::cout << "use selectcells [inFile] [outFile] [MaxX] [MaxY] [MaxZ]\n" << std::endl;
        return -1;
    }

    std::string inFileName    = argv[1];
    std::string outFileName   = argv[2];
    unsigned long nx        = atof(argv[3]);
    unsigned long ny        = atof(argv[4]);
    unsigned long nz        = atof(argv[5]);

    std::ifstream infile(inFileName);
    std::ofstream outfile (outFileName);

    if (infile.fail() || outfile.fail())
    {
        throw std::string("cannot open file");
    }
    std::string line = "";
    unsigned long linesloaded = 0;

    while(std::getline(infile, line))   // copy lines
    {
        outfile << line << "\n";
    }

    std::cout << "writing border voxels\n";

    for (int i = -1; i !=nx+1; ++i)
    for (int j = -1; j !=ny+1; ++j)
    {
        if (i% 4 != 0) continue;
        if (j% 4 != 0) continue;
        outfile << 0 << " " << i << " " << j << " " << -1 << " " << " 1 \n";
        outfile << 0 << " " << i << " " << j << " " << nz << " " << " 1 \n";
    }
    for (int i = -1; i !=nx+1; ++i)
    for (int k = -1; k !=nz+1; ++k)
    {
        if (i% 4 != 0) continue;
        if (k% 4 != 0) continue;
        outfile << 0 << " " << i << " " << -1 << " " << k << " " << " 1 \n";
        outfile << 0 << " "  << i << " " << ny << " " << k << " " << " 1 \n";
    }
    for (int j = -1; j !=ny+1; ++j)
    for (int k = -1; k !=nz+1; ++k)
    {
        if (j% 4 != 0) continue;
        if (k% 4 != 0) continue;
        outfile << 0 << " " << -1 << " " << j << " " << k << " " << " 1 \n";
        outfile << 0 << " " << nx << " "  << j << " " << k << " " << " 1 \n";
    }

    infile.close();
    outfile.close();
}
