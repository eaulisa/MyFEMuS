#ifndef FEMUS_ALGOIM_HPP
#define FEMUS_ALGOIM_HPP

#include <vector>

namespace FemusAlgoim {

// Surface modes use the reference gradient, a constant inverse Jacobian, or geometryCoordinates[physical][node] including the volume determinant.
void IntegrateNodalMoments(const std::vector<double>& nodalValues, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder = 20, unsigned subdivisions = 0, bool referenceSurface = false, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<double>>* geometryCoordinates = nullptr);
void IntegrateQuad9Moments(const std::vector<std::vector<double>>& coefficient, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder = 20, unsigned subdivisions = 0, bool referenceSurface = false, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<double>>* geometryCoordinates = nullptr);
void IntegrateTri7Moments(const std::vector<double>& coefficient, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder = 20, unsigned subdivisions = 0, bool referenceSurface = false, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<double>>* geometryCoordinates = nullptr);
void IntegrateHex27Moments(const std::vector<std::vector<std::vector<double>>>& coefficient, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder = 20, unsigned subdivisions = 0, bool referenceSurface = false, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<double>>* geometryCoordinates = nullptr);

// TET15 moments use cumulative coordinates x=xi+eta+zeta, y=eta+zeta, z=zeta with 0<=z<=y<=x<=1.
void IntegrateTet15NodalMoments(const std::vector<double>& nodalValues, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder = 20, unsigned subdivisions = 0, bool referenceSurface = false, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<double>>* geometryCoordinates = nullptr);

}

#endif

