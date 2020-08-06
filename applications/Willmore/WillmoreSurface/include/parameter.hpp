
class Parameter {
  public:
    Parameter(const std::string name1, const unsigned simulation1, const bool surface1,
              const bool constraintIsOn1, const unsigned numberOfUniformLevels1,
              const unsigned numberOfSmoothingSteps1,
              const bool finalSmoothIsOn1, const unsigned numberOfNonLinearSteps1,
              const unsigned numberOfIterations1, const double beltramiNorm1 = 0.
             ):
      name(name1),
      simulation(simulation1),
      surface(surface1),
      constraintIsOn(constraintIsOn1),
      numberOfUniformLevels(numberOfUniformLevels1),
      numberOfSmoothingSteps(numberOfSmoothingSteps1),
      finalSmoothIsOn(finalSmoothIsOn1),
      numberOfNonLinearSteps(numberOfNonLinearSteps1),
      numberOfIterations(numberOfIterations1),
      beltramiNorm(beltramiNorm1) {
    }
    Parameter() {}

    std::string name;
    unsigned simulation;
    bool surface;
    bool constraintIsOn;
    unsigned numberOfUniformLevels;
    unsigned numberOfSmoothingSteps;
    bool finalSmoothIsOn;
    unsigned numberOfNonLinearSteps;
    unsigned numberOfIterations;
    double beltramiNorm;

    void print() {
      std::cout << "\nBenchmarck values for " << name << std::endl;
      std::cout << "simultation number = " << simulation << std::endl;
      std::cout << "surface = " << ((surface) ? "true" : "false") << std::endl;
      std::cout << "contraint is on = " << ((constraintIsOn) ? "true" : "false") << std::endl;
      std::cout << "number of uniform refinements = " << numberOfUniformLevels << std::endl;
      std::cout << "number of smoothing steps = " << numberOfSmoothingSteps << std::endl;
      std::cout << "final smooth is on = " << ((finalSmoothIsOn) ? "true" : "false") << std::endl;
      std::cout << "number of non linear steps = " << numberOfNonLinearSteps << std::endl;
      std::cout << "number of time iterations = " << numberOfIterations << std::endl;
      std::cout << "Beltrami norm = " << beltramiNorm << std::endl;
    }
};
