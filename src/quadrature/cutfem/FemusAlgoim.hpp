#ifndef FEMUS_ALGOIM_HPP
#define FEMUS_ALGOIM_HPP

#include <vector>

namespace FemusAlgoim {

// Return paired canonical moments in the total-degree ordering used by CutFemWeight.
void IntegrateNodalMoments(const std::vector<double>& nodalValues, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder = 20, unsigned subdivisions = 0);
void IntegrateQuad9Moments(const std::vector<std::vector<double>>& coefficient, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder = 20, unsigned subdivisions = 0);
void IntegrateTri7Moments(const std::vector<double>& coefficient, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder = 20, unsigned subdivisions = 0);
void IntegrateHex27Moments(const std::vector<std::vector<std::vector<double>>>& coefficient, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, unsigned quadratureOrder = 20, unsigned subdivisions = 0);

}

#endif
