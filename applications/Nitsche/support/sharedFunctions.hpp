#ifndef __femus_sharedFunctions_hpp__
#define __femus_sharedFunctions_hpp__

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "PetscMatrix.hpp"
#include "slepceps.h"

#include "LinearImplicitSystem.hpp"
#include "Marker.hpp"
#include "Line.hpp"

#include "Fluid.hpp"
#include "Solid.hpp"
#include "Parameter.hpp"

using namespace femus;

void GetInterfaceElementEigenvalues(MultiLevelSolution& mlSol, Line* line3, Line* lineI, const double &deps);
void GetInterfaceElementEigenvaluesAD(MultiLevelSolution& mlSol, Line* line3, Line* lineI, const double &deps);
void GetParticleWeights(MultiLevelSolution & mlSol, Line* line3);

#endif
 
