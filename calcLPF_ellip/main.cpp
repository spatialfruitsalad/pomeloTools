#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <cmath>
#include <numeric>
#include <algorithm>

#include "splitstring.hpp"

double cx = 566;
double cy = 566;
double rmax = 300;
double zmin = 200;
double zmax = 800;

bool checkpoint (double x, double y, double z)
{
    if (z < zmin || z > zmax) return false;

    double dx = cx - x;
    double dy = cy - y;
    double r = sqrt(dx*dx + dy*dy);
    if (r > rmax) return false;

    return true;
}


int main (int argc, char* argv[])
{
    if (argc != 10)
    {
        std::cout << "This is calcLPF_ellip, a tool twhich can be used to calculate local packing fractions from pomelo in- and output\nA cylindrical container is assumed.\nThe parameters cx, cy, rmax, zmin, zmax can be used to discard particles which are positoined at the container boundary.\n\n";
        std::cout << "[setVoronoiVolumes.dat -- path to the setVoronoiVolumes.dat file, created by pomelo]" << std::endl;
        std::cout << "[ellipFileName -- path to the ellip file, which was used as an input for pomelo]" << std::endl;
        std::cout << "[outFolder -- path to output folder. The output of this program will be written here. ]" << std::endl;
        std::cout << "[infoName -- tag to identify this analysis]" << std::endl;
        std::cout << "[cx cy rmax zmin zmax -- parameters for discarding particles at the container boundary]" << std::endl;
        std::cerr << "ERROR: cmdl arguments not matching!" << std::endl;
        return -1;
    }
    
    std::string setVoronoiVolumes = argv[1];
    std::string ellipFileName = argv[2];
    std::string outFolder    = argv[3];
    std::string infoName     = argv[4];
    cx = std::stod(argv[5]);
    cy = std::stod(argv[6]);
    rmax = std::stod(argv[7]);
    zmin = std::stod(argv[8]);
    zmax = std::stod(argv[9]);
   

///////////////////////
// Parse ellip data
///////////////////////
std::map<unsigned long, double> ellipVolumes;
{
    std::cout << "# parse ellip data" << std::endl;
    std::ifstream infile;
    infile.open(ellipFileName);
    if (infile.fail())
    {
        throw std::string("cannot open ellip input file");
    }
    std::string line = "";
    unsigned long linesloaded = 0;

    unsigned long erased = 0;
    while(std::getline(infile, line))   // parse lines
    {
        if(line.find('#') != std::string::npos) continue;

        unsigned long label;
        double x, y, z;
        double a1, a2, a3;
        splitstring ss(line.c_str());
        auto ssv = ss.split(' ');
        label = stoi(ssv[0]);
        x = stod(ssv[1]);
        y = stod(ssv[2]);
        z = stod(ssv[3]);
        a1 = stod(ssv[4]);
        a2 = stod(ssv[8]);
        a3 = stod(ssv[12]);

        linesloaded++;
        if (!checkpoint(x,y,z))
        {
            erased++;
            continue;
        }
        
        double v = 4.0/3.0*3.141592*a1*a2*a3;

        ellipVolumes[label] = v;
    }
    std::cout << "outside ellipsoids: " << erased << std::endl;

}

    std::cout << "# parsed tetrahedra: " << ellipVolumes.size() << std::endl;

///////////////////////
// Parse CellVolumes
///////////////////////
    std::map<unsigned long, double> cellVolumes;
{
    std::cout << "# parse cell volume data" << std::endl;
    std::ifstream infile;
    infile.open(setVoronoiVolumes);
    if (infile.fail())
    {
        throw std::string("cannot open input file");
    }
    std::string line = "";
    while(std::getline(infile, line))   // parse lines
    {
        if(line.find('#') != std::string::npos) continue;
        
        unsigned long l = 0;
        double v = 0;


        std::istringstream iss(line);
        if (!(iss >> l >> v))
        {
            std::cerr << "error parsing line" << std::endl;
        }

        if (ellipVolumes.count(l) == 1)
        {
            cellVolumes[l] = v;    
        }
    }
}

std::cout << "done" << std::endl;

///////////////////////
// calc particle volume distribution
///////////////////////
double totalParticleVolume = 0;
double totalCellColume = 0;
{
    std::cout << "# calculating ellipsoid volume distribution ... " << std::flush;
    // particleVolumeList
    std::vector<double> pVL;

    std::ofstream volumeFile ( outFolder + "/" + infoName + "_voldist.dat");

    for (auto it = ellipVolumes.begin(); it != ellipVolumes.end(); ++it)
    {
        pVL.push_back(it->second);
        totalParticleVolume += it->second;
    }

    double minV = *std::min_element(pVL.begin(), pVL.end());
    double maxV = *std::max_element(pVL.begin(), pVL.end());

    const unsigned long bins = ((fabs(maxV-minV)< 0.1)? 1 :  9);
    volumeFile << "# minV= " << minV << " maxV= " << maxV << std::endl;
 
    std::vector<unsigned int> hist(bins+1,0);   
    for (double v : pVL)
    {
        unsigned int index = static_cast<unsigned int>((v - minV)/(maxV-minV)*bins);
        hist[index] += 1;
    }

    {
        for (size_t i = 0; i != hist.size(); ++i)
        {
            volumeFile  << minV +  i*(maxV-minV)/bins << " " << hist[i] << std::endl;
        }
    }

    std::cout << "# done" << std::endl;

}

///////////////////////
// calc mean and std
///////////////////////
 
    
    std::ofstream lpfFile(outFolder + "/" + infoName + "_lpf.dat");
    std::vector<double> labels;
    std::vector<double> lpfs;
    double inversemean = 0;
    double cellVolumeTotal = 0;
    double tetraVolumeTotal = 0;

    for (auto it = ellipVolumes.begin(); it != ellipVolumes.end(); ++it)
    {
        double cv = cellVolumes[it->first];
        totalCellColume +=cv;
        lpfs.push_back(it->second/cv);
        labels.push_back(it->first);
        inversemean += 1.0/(it->second/cv);
        cellVolumeTotal += cv;
        tetraVolumeTotal += it->second;
    }
    inversemean /=labels.size();
    double harmonicmean = 1./inversemean;
    double phiG = tetraVolumeTotal/cellVolumeTotal;

    double min = *std::min_element(lpfs.begin(), lpfs.end());
    double max = *std::max_element(lpfs.begin(), lpfs.end());
    auto size = lpfs.size();
    std::cout << "# number of elements: " << size << "\n#min: " << min << "\n#max: " << max << std::endl;
    lpfFile << "# number of elements: " << size << "\n#min: " << min << "\n#max: " << max << std::endl;
    lpfFile << "# totalParticleVolume: " << totalParticleVolume << std::endl;
    lpfFile << "# totalCellVoume: " << totalCellColume << std::endl;

    lpfFile << "# meanParticleVolume: " << totalParticleVolume/size << std::endl;

    double sum = std::accumulate(lpfs.begin(), lpfs.end(),0.0);
    double mean = sum/size;
    std::cout << "# sum: " << sum << std::endl;
    std::cout << "#mean: " << mean << std::endl;
    std::cout << "#harmonicmean: " << harmonicmean << std::endl;
    std::cout << "#PhiG: " << phiG << std::endl;

    lpfFile << "# sum: " << sum << std::endl;
    lpfFile << "#mean: " << mean << std::endl;
    lpfFile << "#harmonicmean: " << harmonicmean << std::endl;
    lpfFile << "#PhiG: " << phiG << std::endl;

    double sum_sq = std::inner_product(lpfs.begin(), lpfs.end(), lpfs.begin(), 0.0);
    double std = std::sqrt(fabs(sum_sq/size - mean*mean));

    std::cout << "# sum_sq: " << sum_sq << std::endl;
    std::cout << "# sum_sq/size: " << sum_sq/size << std::endl;
    std::cout << "# mean*mean: " << mean*mean << std::endl;
    std::cout << "# std: " << std << std::endl;
    lpfFile << "# sum_sq: " << sum_sq << std::endl;
    lpfFile << "# sum_sq/size: " << sum_sq/size << std::endl;
    lpfFile << "# mean*mean: " << mean*mean << std::endl;
    lpfFile << "# std: " << std << std::endl;
    

///////////////////////
// do output
///////////////////////
    std::cout << "#1_Label #2_LPF #3_CellVolume"  << std::endl;
    for (size_t i = 0; i != lpfs.size(); ++i)
    {
        lpfFile << labels[i] << " " << lpfs[i] << " " << cellVolumes[labels[i]] << std::endl;
    }


    return 0;
}
