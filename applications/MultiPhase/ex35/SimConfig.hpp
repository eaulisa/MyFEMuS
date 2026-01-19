 
#include <string>
#include <vector>
#include <stdexcept>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>

struct SimConfig {
  RungeKutta::VelKind velocityType = RungeKutta::VelKind::Vortex;

  MeshSeedFactory::Type    meshSeed  = MeshSeedFactory::Type::SquareQuad9;
  MeshSeedFactory::Reshape meshShape = MeshSeedFactory::Reshape::Identity;

  unsigned dim = 2; // derived from meshSeed

  std::vector<double> shift; // size dim
  std::vector<double> xc;    // bubble center size dim
  double r = 0.15;

  double period = 8.0;       // derived from dim+velocityType
  unsigned uniformRefinementLevel = 1;

  unsigned print_step = 10;
  unsigned reinit_step = std::numeric_limits<unsigned>::max();
  unsigned levelNstart = 8;  // default 2D=8, 3D=6
  unsigned delta_depth = 5;  // default 2D=5, 3D=4

  unsigned nSteps = 320;     // default, user can shorten
};

static inline unsigned dimFromSeed(MeshSeedFactory::Type t) {
  switch (t) {
    case MeshSeedFactory::Type::SquareQuad9:
    case MeshSeedFactory::Type::SquareTri7:  return 2;
    case MeshSeedFactory::Type::CubeHex27:
    case MeshSeedFactory::Type::CubeWedge21:
    case MeshSeedFactory::Type::CubeTet15:   return 3;
  }
  throw std::runtime_error("Unknown mesh seed type");
}


static inline void applyDefaults(SimConfig& c) {
  c.dim = dimFromSeed(c.meshSeed);

  // defaults by dimension
  if (c.dim == 2) {
    if (c.levelNstart == 0) c.levelNstart = 8;
    if (c.delta_depth == 0) c.delta_depth = 5;
    if (c.shift.empty()) c.shift = {0.0, 0.0};
    if (c.xc.empty())    c.xc    = {0.0, 0.25};
  } else {
    if (c.levelNstart == 0) c.levelNstart = 6;
    if (c.delta_depth == 0) c.delta_depth = 4;
    if (c.shift.empty()) c.shift = {0.0, 0.0, 0.0};
    if (c.xc.empty())    c.xc    = {0.0, 0.0, 0.25};
  }

  // period rules
  if (c.velocityType == RungeKutta::VelKind::Rotation) {
    c.period = 2.0 * M_PI;
  } else {
    c.period = (c.dim == 2 ? 8.0 : 4.0);
  }

  // rising bubble shift rule
  if (c.velocityType == RungeKutta::VelKind::RisingBubble) {
    if (c.dim == 2) {
      // needs +0.5 in y (unless user already overrode shift)
      if (c.shift.size() != 2) c.shift = {0.0, 0.0};
      c.shift[1] = 0.5;
    } else {
      if (c.shift.size() != 3) c.shift = {0.0, 0.0, 0.0};
      c.shift[2] = 0.5;
    }
  } else {
    // vortex/rotation default box => shift 0
    // (only if user didn't override shift explicitly)
    // If you want "always force zero", remove this guard.
  }

  // uniform refinement default tweaks you mentioned
  if (c.meshSeed == MeshSeedFactory::Type::CubeTet15) {
    if (c.uniformRefinementLevel < 2) c.uniformRefinementLevel = 2;
  }

  // final sanity
  if (c.shift.size() != c.dim) throw std::runtime_error("shift has wrong size");
  if (c.xc.size()    != c.dim) throw std::runtime_error("xc has wrong size");
}


static inline bool starts_with(const char* s, const char* pref) {
  return std::strncmp(s, pref, std::strlen(pref)) == 0;
}

static inline std::string toLower(std::string s) {
  for (auto& ch : s) ch = char(std::tolower((unsigned char)ch));
  return s;
}

static inline RungeKutta::VelKind parseVel(const std::string& v) {
  const std::string s = toLower(v);
  if (s == "vortex") return RungeKutta::VelKind::Vortex;
  if (s == "rotation") return RungeKutta::VelKind::Rotation;
  if (s == "rising" || s == "risingbubble" || s == "bubble") return RungeKutta::VelKind::RisingBubble;
  throw std::runtime_error("Unknown --vel");
}

static inline MeshSeedFactory::Reshape parseShape(const std::string& v) {
  const std::string s = toLower(v);
  if (s == "identity" || s == "box") return MeshSeedFactory::Reshape::Identity;
  if (s == "ball") return MeshSeedFactory::Reshape::Ball;
  if (s == "funnel") return MeshSeedFactory::Reshape::Funnel;
  throw std::runtime_error("Unknown --shape");
}

static inline MeshSeedFactory::Type parseSeed(const std::string& v) {
  const std::string s = toLower(v);
  if (s == "tri7")  return MeshSeedFactory::Type::SquareTri7;
  if (s == "quad9") return MeshSeedFactory::Type::SquareQuad9;
  if (s == "hex27") return MeshSeedFactory::Type::CubeHex27;
  if (s == "wedge21") return MeshSeedFactory::Type::CubeWedge21;
  if (s == "tet15") return MeshSeedFactory::Type::CubeTet15;
  throw std::runtime_error("Unknown --mesh");
}

static inline void parseVec(std::vector<double>& out, const char* s) {
  // expects comma-separated "a,b" or "a,b,c"
  out.clear();
  std::string str(s);
  std::size_t pos = 0;
  while (true) {
    std::size_t comma = str.find(',', pos);
    std::string tok = (comma == std::string::npos) ? str.substr(pos) : str.substr(pos, comma - pos);
    out.push_back(std::stod(tok));
    if (comma == std::string::npos) break;
    pos = comma + 1;
  }
}



static inline void printHelp(const char* prog) {
  std::cout
  << "\n"
  << "Usage:\n"
  << "  " << prog << " [options]\n"
  << "\n"
  << "Goal:\n"
  << "  Run a 2D or 3D simulation choosing velocity test, mesh seed, and domain reshape.\n"
  << "  Most defaults depend on dimension (2D vs 3D) and velocity test.\n"
  << "\n"
  << "Options (all optional):\n"
  << "  --help\n"
  << "      Print this help and exit.\n"
  << "\n"
  << "  --vel <vortex|rotation|rising|risingBubble|bubble>\n"
  << "      Velocity test.\n"
  << "      Default: vortex\n"
  << "\n"
  << "  --mesh <quad9|tri7|hex27|wedge21|tet15>\n"
  << "      Mesh seed type.\n"
  << "      2D: quad9, tri7\n"
  << "      3D: hex27, wedge21, tet15\n"
  << "      Default: quad9 (2D)\n"
  << "\n"
  << "  --shape <box|identity|ball|funnel>\n"
  << "      Domain reshape applied to the default box [-0.5,0.5]^dim.\n"
  << "      Default: box (identity)\n"
  << "\n"
  << "  --steps <N>\n"
  << "      Number of time steps.\n"
  << "      Default: 320\n"
  << "\n"
  << "  --print <k>\n"
  << "      Output frequency (print every k steps).\n"
  << "      Default: 10\n"
  << "\n"
  << "  --uniformRef <L>\n"
  << "      Uniform refinement level before adaptive refinement.\n"
  << "      Default: 1 (but tet15 is forced to at least 2)\n"
  << "\n"
  << "  --levelNstart <L>\n"
  << "      Starting AMR level (depends on dimension).\n"
  << "      Default: 8 in 2D, 6 in 3D\n"
  << "\n"
  << "  --deltaDepth <D>\n"
  << "      AMR depth increment (depends on dimension).\n"
  << "      Default: 5 in 2D, 4 in 3D\n"
  << "\n"
  << "  --shift a,b[,c]\n"
  << "      Shift applied to the domain/mesh seed.\n"
  << "      Default depends on velocity test:\n"
  << "        - vortex/rotation: (0,0) in 2D or (0,0,0) in 3D\n"
  << "        - rising bubble:   forced to (0,0.5) in 2D or (0,0,0.5) in 3D\n"
  << "      Note: shift length must match dimension.\n"
  << "\n"
  << "  --xc a,b[,c]\n"
  << "      Bubble center location.\n"
  << "      Default: (0,0.25) in 2D, (0,0,0.25) in 3D\n"
  << "      Note: xc length must match dimension.\n"
  << "\n"
  << "  --r <radius>\n"
  << "      Bubble radius.\n"
  << "      Default: 0.15\n"
  << "\n"
  << "  --period <T>\n"
  << "      Override the default period.\n"
  << "      Default rule:\n"
  << "        - rotation: 2*pi\n"
  << "        - vortex or rising bubble: 8 (2D) or 4 (3D)\n"
  << "\n"
  << "  --reinit <n>\n"
  << "      Enable field reinitialization each <n> iterations.\n"
  << "      Default rule: off"
  << "\n"
  << "Defaults summary (if you do not set anything):\n"
  << "  --vel vortex\n"
  << "  --mesh quad9  (=> 2D)\n"
  << "  --shape box\n"
  << "  --steps 320\n"
  << "  --print 10\n"
  << "  --uniformRef 1\n"
  << "  2D defaults: levelNstart=8, deltaDepth=5, xc=(0,0.25), shift=(0,0)\n"
  << "  3D defaults: levelNstart=6, deltaDepth=4, xc=(0,0,0.25), shift=(0,0,0)\n"
  << "  Rising bubble forces shift to +0.5 in the vertical direction (y in 2D, z in 3D).\n"
  << "\n"
  << "Examples:\n"
  << "  # 2D vortex on tri7 in a ball domain, short run\n"
  << "  " << prog << " --vel vortex --mesh tri7 --shape ball --steps 80\n"
  << "\n"
  << "  # 2D rising bubble: shift is forced to (0,0.5)\n"
  << "  " << prog << " --vel rising --mesh quad9 --steps 160\n"
  << "\n"
  << "  # 3D rotation on hex27 in funnel domain\n"
  << "  " << prog << " --vel rotation --mesh hex27 --shape funnel\n"
  << "\n"
  << "  # Override bubble center and radius explicitly\n"
  << "  " << prog << " --vel rising --mesh wedge21 --xc 0,0,0.3 --r 0.12\n"
  << "\n";
}





static inline SimConfig parseArgs(int argc, char** argv) {
  SimConfig c;
  c.levelNstart = 0;
  c.delta_depth = 0;

  auto readValue = [&](int& i) -> std::string {
    // expects argv[i] is either --key=value OR --key value
    std::string a = argv[i];
    auto eq = a.find('=');
    if (eq != std::string::npos) {
      return a.substr(eq + 1);
    }
    // no '=', so value must be next token
    if (i + 1 >= argc) throw std::runtime_error("Missing value for " + a);
    return std::string(argv[++i]);
  };

  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];

    if (a == "--help" || a == "-h") {
      printHelp(argv[0]);
      std::exit(0);
    }

    // Handle --key=value by stripping to key
    std::string key = a;
    auto eq = key.find('=');
    if (eq != std::string::npos) key = key.substr(0, eq);

    if (key == "--vel") {
      const std::string v = readValue(i);
      c.velocityType = parseVel(v);
    }
    else if (key == "--mesh") {
      const std::string v = readValue(i);
      c.meshSeed = parseSeed(v);
    }
    else if (key == "--shape") {
      const std::string v = readValue(i);
      c.meshShape = parseShape(v);
    }
    else if (key == "--shift") {
      const std::string v = readValue(i);
      parseVec(c.shift, v.c_str());
    }
    else if (key == "--xc") {
      const std::string v = readValue(i);
      parseVec(c.xc, v.c_str());
    }
    else if (key == "--r") {
      const std::string v = readValue(i);
      c.r = std::stod(v);
    }
    else if (key == "--period") {
      const std::string v = readValue(i);
      c.period = std::stod(v);
    }
    else if (key == "--steps") {
      const std::string v = readValue(i);
      c.nSteps = (unsigned)std::stoul(v);
    }
    else if (key == "--print") {
      const std::string v = readValue(i);
      c.print_step = (unsigned)std::stoul(v);
    }
    else if (key == "--levelNstart") {
      const std::string v = readValue(i);
      c.levelNstart = (unsigned)std::stoul(v);
    }
    else if (key == "--deltaDepth") {
      const std::string v = readValue(i);
      c.delta_depth = (unsigned)std::stoul(v);
    }
    else if (key == "--uniformRef") {
      const std::string v = readValue(i);
      c.uniformRefinementLevel = (unsigned)std::stoul(v);
    }
    else if (key == "--reinitFixed") {
      const std::string v = readValue(i);
      c.reinit_step = (unsigned)std::stoul(v);
    }
    else {
      throw std::runtime_error("Unknown option: " + a);
    }
  }

  applyDefaults(c);
  return c;
}

