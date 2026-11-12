#ifndef FEMUS_EX45_ALGOIM_QUADRATURE_HPP
#define FEMUS_EX45_ALGOIM_QUADRATURE_HPP

#include <vector>

namespace DirectAlgoim {

// Points use standard FEM reference coordinates; weights integrate dVref and dSref with half measure on exact element-boundary interfaces, and unit normals point toward positive G.
void BuildQuadrature(const std::vector<double>& nodalValues, std::vector<std::vector<double>>& positivePoints, std::vector<double>& positiveWeights, std::vector<std::vector<double>>& negativePoints, std::vector<double>& negativeWeights, std::vector<std::vector<double>>& surfacePoints, std::vector<double>& surfaceWeights, std::vector<std::vector<double>>& surfaceNormals, unsigned quadratureOrder = 20, unsigned subdivisions = 0);

}

#endif
