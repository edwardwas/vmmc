/*
  Copyright (c) 2015 Lester Hedges <lester.hedges+vmmc@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifdef ISOTROPIC
#error patchy_disc.cpp cannot be linked to isotropic VMMC library!
#endif

#include "Demo.h"
#include "VMMC.h"

void writeStartOptions( const std::string fileName
        , const double boxSize
        , const double energy
        , const double patchDiameter
        , const int numParts
        , const int writeRate 
        , const int numSweeps){
    std::ofstream fileStream;
    fileStream.open(fileName,std::ios::trunc);
    fileStream << "{\n";

    fileStream << "theta = ["
        << 0 << ","
        << 0.5*M_PI << ","
        << M_PI << ","
        << 1.5*M_PI << "]\n";
    fileStream << "boxSize = " << boxSize << "\n";
    fileStream << "sameEnergy = " << energy << "\n";
    fileStream << "differentEnergy = " << energy << "\n";
    fileStream << "patchDiameter = " << patchDiameter << "\n";
    fileStream << "particleDiameter = " << 1 << "\n";
    fileStream << "numParts = " << numParts << "\n";
    fileStream << "fileName = \"" << fileName << "\"\n";
    fileStream << "writeRate = " << writeRate << "\n";
    fileStream << "numSweeps = " << numSweeps << "\n";
    fileStream << "singleMoveMonteCarlo = " << false << "\n";
    fileStream << "}\n";
    fileStream.close();
}

void writeFrame(const std::string fileName, const std::vector<Particle> parts){
    std::ofstream fileStream;
    fileStream.open(fileName,std::ios::app);

    fileStream << parts.size() << "\n\n";
    for (int i = 0; i < parts.size(); i++){
        fileStream << parts[i].position[0] << " ";
        fileStream << parts[i].position[1] << " ";
        fileStream << parts[i].orientation[0] << " ";
        fileStream << parts[i].orientation[1] << " ";
        fileStream << "Red Blue Red Blue\n";
    }
    fileStream << "\n";
}

int main(int argc, char** argv)
{
    int crystalSize = 9;
    // Simulation parameters.
    unsigned int dimension = 2;                     // dimension of simulation box
    unsigned int nParticles = crystalSize*crystalSize;                 // number of particles
    double interactionEnergy = 3.0;                 // pair interaction energy scale (in units of kBT)
    double interactionRange = 0.15;                  // diameter of patch (in units of particle diameter)
    double baseLength = 21;                              // base length of simulation box
    unsigned int maxInteractions = 4;               // maximum number of interactions per particle (number of patches)

    const int numSweeps = 5e4;
    const int writeRate = numSweeps / 100;

    // Data structures.
    std::vector<Particle> particles(nParticles);    // particle container
    CellList cells;                                 // cell list
    bool isIsotropic[nParticles];                   // whether the potential of each particle is isotropic

    // Resize particle container.
    particles.resize(nParticles);

    // Work out base length of simulation box (particle diameter is one).
    //if (dimension == 2) baseLength = std::pow((nParticles*M_PI)/(2.0*density), 1.0/2.0);
    //else baseLength = std::pow((nParticles*M_PI)/(6.0*density), 1.0/3.0);

    std::vector<double> boxSize;
    for (unsigned int i=0;i<dimension;i++)
        boxSize.push_back(baseLength);

    // Initialise simulation box object.
    Box box(boxSize);

    // Initialise cell list.
    cells.setDimension(dimension);
    cells.initialise(box.boxSize, 1 + 0.5*interactionRange);

    // Initialise the patchy disc model.
    PatchyDisc patchyDisc(box, particles, cells,
        maxInteractions, interactionEnergy, interactionRange);

    // Initialise random number generator.
    MersenneTwister rng;

    // Initialise particle initialisation object.
    Initialise initialise;

    // Generate a random particle configuration.
    initialise.random(particles, cells, box, rng, false);
    cells.reset();

    int i = 0;
    double deltaX = 1 + 0.5*interactionRange;
    double startX = 0.5 * (baseLength - crystalSize * deltaX);
    for (int x = 0; x < crystalSize; x ++){
        for (int y = 0; y < crystalSize; y++){
            std::cout << particles[i].cell << "-";
            particles[i].position[0] = startX + x * deltaX;
            particles[i].position[1] = startX + y * deltaX;
            particles[i].orientation[0] = 1;
            particles[i].orientation[1] = 0;
            particles[i].cell = cells.getCell(particles[i]);
            cells.initCell(particles[i].cell, particles[i]);
            std::cout << particles[i].cell << "\n";
            i++;
        }
    }

    // Initialise data structures needed by the VMMC class.
    double coordinates[dimension*nParticles];
    double orientations[dimension*nParticles];

    // Copy particle coordinates and orientations into C-style arrays.
    for (unsigned int i=0;i<nParticles;i++)
    {
        for (unsigned int j=0;j<dimension;j++)
        {
            coordinates[dimension*i + j] = particles[i].position[j];
            orientations[dimension*i + j] = particles[i].orientation[j];
        }

        // Set all particles as anisotropic.
        isIsotropic[i] = false;
    }

    // Initialise the VMMC callback functions.
    using namespace std::placeholders;
    vmmc::CallbackFunctions callbacks;
#ifndef ISOTROPIC
    callbacks.energyCallback =
        std::bind(&PatchyDisc::computeEnergy, patchyDisc, _1, _2, _3);
    callbacks.pairEnergyCallback =
        std::bind(&PatchyDisc::computePairEnergy, patchyDisc, _1, _2, _3, _4, _5, _6);
    callbacks.interactionsCallback =
        std::bind(&PatchyDisc::computeInteractions, patchyDisc, _1, _2, _3, _4);
    callbacks.postMoveCallback =
        std::bind(&PatchyDisc::applyPostMoveUpdates, patchyDisc, _1, _2, _3);
#else
    callbacks.energyCallback =
        std::bind(&PatchyDisc::computeEnergy, patchyDisc, _1, _2);
    callbacks.pairEnergyCallback =
        std::bind(&PatchyDisc::computePairEnergy, patchyDisc, _1, _2, _3, _4);
    callbacks.interactionsCallback =
        std::bind(&PatchyDisc::computeInteractions, patchyDisc, _1, _2, _3);
    callbacks.postMoveCallback =
        std::bind(&PatchyDisc::applyPostMoveUpdates, patchyDisc, _1, _2);
#endif

    // Initialise VMMC object.
    vmmc::VMMC vmmc(nParticles, dimension, coordinates, orientations,
        0.15, 0.2, 0.5, 0.5, maxInteractions, &boxSize[0], isIsotropic, false, callbacks);

    std::cout << "Writing start options\n";
    writeStartOptions("out.txt",baseLength,interactionEnergy,interactionRange,nParticles,writeRate,numSweeps);
    std::cout << "Start options writen\n";


    // Execute the simulation.
    for (unsigned int i=0;i<numSweeps/writeRate;i++)
    {
        vmmc += writeRate*nParticles;
        writeFrame("out.txt",particles);

        std::cout << "After " << (100.0 * double(i*writeRate) / double(numSweeps)) << "% complete, the energy is " << patchyDisc.getEnergy() << "\n";
    }

    std::cout << "\nComplete!\n";

    // We're done!
    return (EXIT_SUCCESS);
}
