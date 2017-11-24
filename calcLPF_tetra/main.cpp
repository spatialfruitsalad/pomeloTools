#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <cmath>
#include <numeric>
#include <algorithm>

#include "splitstring.hpp"


double calculateLength (double x1, double y1, double z1, double x2, double y2, double z2)
{
    double dx = x2-x1;
    double dy = y2-y1;
    double dz = z2-z1;

    double l = sqrt( dx*dx + dy*dy + dz*dz);

    return l;
}

bool checkpoint (double x, double y, double z)
{
    if (z < 200 || z > 800) return false;

    double cx = 640;
    double cy = 640;

    double dx = cx - x;
    double dy = cy - y;
    double r = sqrt(dx*dx + dy*dy);
    if (r > 400) return false;

    return true;
}


int main (int argc, char* argv[])
{
    if (argc != 5)
    {
        std::cout << "[w000FileName ]" << std::endl;
        std::cout << "[dataFileName ]" << std::endl;
        std::cout << "[outFolder    ]" << std::endl;
        std::cout << "[infoName     ]" << std::endl;
        std::cerr << "ERROR: cmdl arguments not matching!" << std::endl;
        return -1;
    }
    
    std::string w000FileName = argv[1];
    std::string dataFileName = argv[2];
    std::string outFolder    = argv[3];
    std::string infoName     = argv[4];

   

///////////////////////
// Parse tetra data
///////////////////////
    double meanLength = 0;
std::map<unsigned long, double> tetraVolumes;
{
    std::cout << "# parse tetra data" << std::endl;
    std::ifstream infile;
    infile.open(dataFileName);
    if (infile.fail())
    {
        throw std::string("cannot open tetra input file");
    }
    std::string line = "";
    unsigned long linesloaded = 0;

    std::ofstream lengthFile ( outFolder + "/" + infoName + "_lengths.dat");

    unsigned long erased = 0;
    while(std::getline(infile, line))   // parse lines
    {
        if(line.find('#') != std::string::npos) continue;

        std::istringstream iss(line);
        unsigned long label;
        double x1, y1, z1;
        double x2, y2, z2;
        double x3, y3, z3;
        double x4, y4, z4;
        if (!(iss >> label >> x1 >> y1 >> z1 >> x2 >> y2 >> z2>> x3 >> y3 >> z3 >> x4 >> y4 >> z4 ))
        {
            std::cerr << "ERROR: Parsing tetra file" << std::endl;

            return -2;
        }
        linesloaded++;
        double cx = (x1 + x2 + x3 + x4)/4.0;
        double cy = (y1 + y2 + y3 + y4)/4.0;
        double cz = (z1 + z2 + z3 + z4)/4.0;
        if (!checkpoint(cx,cy,cz))
        {
            erased++;
            continue;
        }
        
        std::vector<double> lines(6,0);
        lines[0] = calculateLength(x1,y1,z1, x2, y2, z2);
        lines[1] = calculateLength(x1,y1,z1, x3, y3, z3);
        lines[2] = calculateLength(x2,y2,z2, x3, y3, z3);
        lines[3] = calculateLength(x1,y1,z1, x4, y4, z4);
        lines[4] = calculateLength(x2,y2,z2, x4, y4, z4);
        lines[5] = calculateLength(x3,y3,z3, x4, y4, z4);
       
        double a = std::accumulate(lines.begin(), lines.end(), 0)/6.0;
        meanLength += lines[0];
        double v = a*a*a /12.0 * std::sqrt(2.0);

        tetraVolumes[label] = v;

        for (double l : lines)
        {
            lengthFile << l << std::endl;
        }
    }
    std::cout << "outside tetrahedra: " << erased << std::endl;

}
meanLength /= tetraVolumes.size();

    std::cout << "# parsed tetrahedra: " << tetraVolumes.size() << std::endl;

///////////////////////
// Parse CellVolumes
///////////////////////
    std::map<unsigned long, double> cellVolumes;
{
    std::cout << "# parse cell volume data" << std::endl;
    std::ifstream infile;
    infile.open(w000FileName);
    if (infile.fail())
    {
        throw std::string("cannot open input file");
    }
    std::string line = "";
    unsigned long linesloaded = 0;
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

        if (tetraVolumes.count(l) == 1)
        {
            cellVolumes[l] = v;    
        }
    }
}


///////////////////////
// Parse w000 data
///////////////////////
if (false)
{
    std::cout << "# parse w000 data" << std::endl;
    std::ifstream infile;
    infile.open(w000FileName);
    if (infile.fail())
    {
        throw std::string("cannot open tetra input file");
    }
    std::string line = "";
    unsigned long linesloaded = 0;
    unsigned long erased = 0;
    while(std::getline(infile, line))   // parse lines
    {
        splitstring ss(line.c_str());
        std::vector<std::string> split = ss.split(' ');
        if (split.size() <= 2) continue;
        if (split[2] != "w000") continue;

        unsigned long l = std::stoi(split[0]);
        if (tetraVolumes.count(l) == 1)
        {
            //std::cout << split[1] << std::endl;
            if(split[1] == "ERROR") 
            {
                erased++;
                tetraVolumes.erase(l);
                continue;
            }            
            double cv = std::stod(split[1]);

            cellVolumes[l] = cv;    
        }
    }
    std::cout << "incorrect karambola calculations: " << erased << std::endl;
}

std::cout << "done" << std::endl;

///////////////////////
// calc particle volume distribution
///////////////////////
double totalParticleVolume = 0;
double totalCellColume = 0;
{
    std::cout << "# calculating tetrahedra volume distribution ... " << std::flush;
    // particleVolumeList
    std::vector<double> pVL;

    std::ofstream volumeFile ( outFolder + "/" + infoName + "_voldist.dat");

    for (auto it = tetraVolumes.begin(); it != tetraVolumes.end(); ++it)
    {
        pVL.push_back(it->second);
        totalParticleVolume += it->second;
    }

    double minV = *std::min_element(pVL.begin(), pVL.end());
    double maxV = *std::max_element(pVL.begin(), pVL.end());

    const unsigned long bins = ((fabs(maxV-minV)< 0.1)? 1 :  9);
    volumeFile << "# minV= " << minV << " maxV= " << maxV << std::endl;
    //double dV = (maxV - minV)/bins;
 
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

    for (auto it = tetraVolumes.begin(); it != tetraVolumes.end(); ++it)
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
    lpfFile << "# meanParticleLength: " << meanLength << std::endl;

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
