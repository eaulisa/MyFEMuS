// Loibl et al. (2026), Comparative Study of Different Quadrature Methods
// for Cut Elements: PAPER-MATCHED 2D/3D BENCHMARK DRIVER for FEMuS/rifat/ex44_test.
//
// Reproduces the 2D input settings of the companion repository, v2.2.0:
//   - Sect. 3.2: three area sweeps, nQP,set=1..4, 1x1 mesh.
//   - Sect. 3: area h-refinement, ne=1,2,4,8,16,32; nQP,set=1..4.
//   - Sect. 3: monomial sweeps q=1..6, nQP,set=1..p*(q+1)+1.
//   - Sect. 4.1.4: quadratic case, q=6, h-refinement, nQP,set=3,6.
//   - Sect. 4.1.1: quintic case, 2x2 mesh, nQP,set=3; point export.
//   - Sect. 4.1.2: all 26 UniBw 2D geometries, ne=1 and ne=8;
//                  nQP,set=ceil((p+1)/2) (FEMuS/Algoim).
//   - Sect. 4.1.3: 1000 moving-parabola positions, ne=8, nQP,set=2.
//   - Sect. 4.2.1 / Fig. 42 / Table 4: moving ellipsoid, 1000 shifts,
//     8x8x8 elements, nQP,set=2, exact quadratic level set.
//   - Sect. 4.2.1 / Fig. 42: optional equivalent moving-background-mesh
//     ellipsoid study (authors' QuESo/STL configuration).
//   - Sect. 4.2.2 / Fig. 44: torus volume, 1..32 elements per direction,
//     nQP,set=2, exact quartic Algoim geometry and HEX27-Q2 FEMuS geometry.
//   - Sect. 4.2.3 / Fig. 46: elliptic cylinder f_q = X^q Y^q Z^q,
//     q=3,5, 8x8x8 mesh, nQP,set=1..4.
//   - Sect. 4.2.3 / Fig. 47: elliptic cylinder q=5, nQP,set=3,
//     1..32 elements per direction.
// NOTE: only FEMuS moment fitting and direct Algoim are implemented. The
// original paper includes other tools; this is not their implementation.
//
// Geometry and monomials are copied/adapted from:
//   B2-M/CutElementIntegration, tag v2.2.0
//   utilities/examples/getInterfaceCase2D.m
//   utilities/examples/getTestCase2D.m
//   publications/Paper_CutElemComp/CutElemComp_example_*.m
//   doi:10.1007/s11831-026-10619-2
//
// IMPORTANT SCIENTIFIC LIMITATIONS:
// - FEMuS uses a QUAD9 quadratic interpolant for Psi: degree 3-5
//   curves, a circle, and piecewise boundaries are therefore APPROXIMATE.
// - FEMuS's FE Gauss table supports up to index 37; the current mapping
//   momentDegree=2*nQP,set-1 needs index=4*nQP,set-2, so FEMuS
//   cannot run nQP,set>9 without changing the FEMuS library. Such rows
//   are explicitly marked SKIPPED; Algoim still runs the paper's sweep.
// - Algoim uses exact polynomial geometry, including subelement splits
//   for suite case 6 and the piecewise curve in case 8.
// - The moving lens is an intersection of TWO parabolic half-domains.
//   Algoim is applied via signed difference of two single-boundary
//   polynomial quadrature rules, NOT as the paper's native multi-interface
//   Algoim construction. Its node counts are thus NOT the paper's counts.
// - 'femus-moment' is the in-house method, NOT BoSSS/mlhp/QuESo.
// - Run on ONE MPI rank: failures are recorded and the sweep continues.
// - Exact polynomial 3D geometry is used by Algoim; FEMuS integrates the
//   HEX27-Q2 interpolant, so the quartic torus has a geometry error.
//   For coarse torus meshes a Q2 cell can miss a fully interior component.
// - 3D polynomial references use an analytical disk-moment formula.
// - These tests are NOT verified by compilation in an actual FEMuS build.
//
// Build in the existing ex44_test CMake target or add a separate target.
// Run: ./MultiPhase_review_2d_3d_benchmark --study area --mode compare
//      ./MultiPhase_review_2d_3d_benchmark --study torus-href --mode compare
//      ./MultiPhase_review_2d_3d_benchmark --study cylinder-monomial --mode compare
//      ./MultiPhase_review_2d_3d_benchmark --study 3d --mode compare

#include "FemusInit.hpp"
#include "MultiLevelMesh.hpp"
#include "MultiLevelSolution.hpp"
#include "ElemTypeEnum.hpp"
#include "Fem.hpp"
#include "IntegrateCut.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <utility>

using namespace femus;

namespace {

struct Options {
  std::string study="area", mode="compare", csv="review_2d.csv";
  std::string pointsCsv="review_2d_points.csv";
  unsigned caseId=0, suiteId=0, qMin=1, qMax=6;
  unsigned nqpMin=1, nqpMax=0, levels=5, hrefNqp=0;
  unsigned algoimSubdivisions=0, steps=1000, stepStart=0, stepEnd=0;
  bool quiet=false;
};
unsigned AsUnsigned(const std::string& s) {
  if(s.empty() || s[0]=='-') throw std::invalid_argument("Expected nonnegative integer: "+s);
  std::size_t pos=0;
  const unsigned long x=std::stoul(s,&pos);
  if(pos!=s.size() || x>std::numeric_limits<unsigned>::max())
    throw std::invalid_argument("Invalid unsigned integer: "+s);
  return static_cast<unsigned>(x);
}
void Help() {
  std::cout <<
    "PAPER-MATCHED Loibl et al. (2026) FEMuS/Algoim 2D + 3D tests\n"
    " --study area|area-href|monomial|monomial-href|layout|suite26|moving\n"
    "         ellipsoid-moving|ellipsoid-mesh-moving|torus-href\n"
    "         cylinder-monomial|cylinder-href\n"
    "         core|3d|all\n"
    "   core: four central 2D studies; 3d: all five 3D studies\n"
    "   all: all twelve 2D+3D study families (potentially very expensive)\n"
    " --case 1|2|3|all       (three main Fig. 3 geometries)\n"
    " --suite-case 1..26|all (Fig. 35 geometry selection)\n"
    " --mode moment|algoim|compare\n"
    " --q-min 1 --q-max 6  (2D powers q=1..6; 3D cylinder q=3,5)\n"
    " --nqp-min 1 --nqp-max N (0/omitted => PAPER range for study)\n"
    " --levels 5           (ne=1,2,4,8,16,32)\n"
    " --href-nqp N        (override paper nQP in h-refinement)\n"
    " --steps 1000 --step-start 0 --step-end N (0 = final step)\n"
    "   ellipsoid shifts d=0..0.55; 2D moving shifts d=0..0.25\n"
    "   ellipsoid-mesh-moving models QuESo's mesh-shift variant (nQP,set=1)\n"
    "   torus 3D: ne=1..32, nQP,set=2\n"
    "   cylinder 3D: ne=8, nQP,set=1..4, q=3,5\n"
    "   cylinder 3D href: ne=1..32, nQP,set=3, q=5\n"
    " --algoim-subdivisions N (default 0: automatic)\n"
    " --csv FILE --points-csv FILE --quiet --help\n"
    " nQP,set is a REQUESTED parameter; quadrature_nodes_cut is ACTUAL nQP.\n"
    " FEMuS moment nQP,set>9 cannot run with the unmodified library and\n"
    " is recorded explicitly as SKIPPED, not silently lowered.\n";
}
Options Parse(int argc,char** argv) {
  Options o;
  for(int i=1;i<argc;++i) {
    std::string k=argv[i];
    if(k=="--help" || k=="-h") {Help(); std::exit(0);}
    if(k=="--quiet") {o.quiet=true;continue;}
    if(i+1>=argc) throw std::invalid_argument("Missing value for "+k);
    std::string v=argv[++i];
    if(k=="--study")o.study=v;
    else if(k=="--mode")o.mode=v;
    else if(k=="--case")o.caseId=(v=="all"?0:AsUnsigned(v));
    else if(k=="--suite-case")o.suiteId=(v=="all"?0:AsUnsigned(v));
    else if(k=="--q-min")o.qMin=AsUnsigned(v);
    else if(k=="--q-max")o.qMax=AsUnsigned(v);
    else if(k=="--nqp-min")o.nqpMin=AsUnsigned(v);
    else if(k=="--nqp-max")o.nqpMax=AsUnsigned(v);
    else if(k=="--levels")o.levels=AsUnsigned(v);
    else if(k=="--href-nqp")o.hrefNqp=AsUnsigned(v);
    else if(k=="--steps")o.steps=AsUnsigned(v);
    else if(k=="--step-start")o.stepStart=AsUnsigned(v);
    else if(k=="--step-end")o.stepEnd=AsUnsigned(v);
    else if(k=="--algoim-subdivisions")o.algoimSubdivisions=AsUnsigned(v);
    else if(k=="--csv")o.csv=v;
    else if(k=="--points-csv")o.pointsCsv=v;
    else throw std::invalid_argument("Unknown option: "+k);
  }
  if(o.study!="area" && o.study!="area-href" && o.study!="monomial" &&
     o.study!="monomial-href" && o.study!="layout" &&
     o.study!="suite26" && o.study!="moving" &&
     o.study!="ellipsoid-moving" && o.study!="ellipsoid-mesh-moving" &&
     o.study!="torus-href" &&
     o.study!="cylinder-monomial" && o.study!="cylinder-href" &&
     o.study!="core" && o.study!="3d" && o.study!="all")
    throw std::invalid_argument("Invalid --study");
  if(o.mode!="moment" && o.mode!="algoim" && o.mode!="compare")
    throw std::invalid_argument("Invalid --mode");
  if(o.caseId>3 || o.suiteId>26 || o.nqpMin<1 ||
     (o.nqpMax && o.nqpMax<o.nqpMin) || o.nqpMax>100 ||
     o.hrefNqp>100 || o.levels>6 || o.qMin<1 || o.qMax>6 ||
     o.qMin>o.qMax || o.steps<2 || o.steps>100000 ||
     o.stepStart>=o.steps || (o.stepEnd && (o.stepEnd>=o.steps || o.stepEnd<o.stepStart)) ||
     o.algoimSubdivisions>32 || o.csv.empty() || o.pointsCsv.empty())
    throw std::invalid_argument("Invalid numerical settings");
  return o;
}

struct PolyPart {
  long double xLow=0, xHigh=1;
  std::vector<long double> xPoly, yPoly;
  int sign=1; // +1 for volume integrals; -1 for signed Algoim subtraction
};
struct Geometry {
  unsigned id=0, curveDegree=1, mainCase=0;
  std::string name;
  double referenceArea=0;
  std::vector<PolyPart> parts;
  bool moving=false;
  bool intersection=false;
};
PolyPart Part(std::vector<long double> x, std::vector<long double> y,
              long double xLow=0, long double xHigh=1, int sign=1) {
  PolyPart p; p.xLow=xLow;p.xHigh=xHigh;
  p.xPoly=std::move(x);p.yPoly=std::move(y);p.sign=sign;
  return p;
}

// Suite numbers 1..26 are those of paper Fig. 35, corresponding to the
// authors' getTestCase2D.m indices 9..34 and getInterfaceCase2D.m IDs 10..35.
// All polynomial coefficients below represent Psi>0 INSIDE (negative of
// the MATLAB phi in most cases). The references are authors' sizeInsideDomain.
Geometry Suite(unsigned id) {
  Geometry g;g.id=id;g.name="suite-"+std::to_string(id);
  switch(id) {
    case 1: g.curveDegree=1;g.mainCase=1;g.referenceArea=.7;
      g.parts={Part({-.2L,1.L},{0.L,-.2L})};break;
    case 2: g.curveDegree=2;g.mainCase=2;g.referenceArea=.75;
      g.parts={Part({.8L,-1.L},{0.L,.3L,-.6L})};break;
    case 3: g.curveDegree=3;g.referenceArea=.341666666666667;
      g.parts={Part({.9L,-1.L,-3.1L,3.9L},{0.L,-1.L})};break;
    case 4: g.curveDegree=4;g.referenceArea=.6;
      g.parts={Part({-.7L,-.3L,1.5L,3.4L,-4.5L},{0.L,1.L})};break;
    case 5: g.curveDegree=5;g.referenceArea=.392578532959995;
      g.parts={Part({.5L,.7L,-1.5L,3.4L,-4.5L,-.8L},{0.L,-1.L})};break;
    case 6: g.curveDegree=2;g.referenceArea=.1125;
      // MATLAB also includes the vertical restriction x>=0.25.
      g.parts={Part({.05L,-.4L,.8L},{0.L,-1.L},.25L,1.L)};break;
    case 7: g.curveDegree=2;g.referenceArea=.286475828547308;
      g.parts={Part({-1.5L,1.L},{0.L,2.6L,-1.7L})};break;
    case 8: g.curveDegree=2;g.referenceArea=.208457587049037;
      // Piecewise parabola, change at x=0.3; two exact polynomial slices.
      g.parts={Part({-.4L,-1.2L,.4L},{0.L,1.L},0.L,.3L),
               Part({-.658L,-.1L,-.4L},{0.L,1.L},.3L,1.L)};break;
    case 9: g.curveDegree=3;g.referenceArea=.195397387787143;
      g.parts={Part({-.2L,-1.L},{0.L,1.L,2.1L,-3.4L})};break;
    case 10: g.curveDegree=4;g.referenceArea=.055443038657972;
      g.parts={Part({-1.2L,1.L},{0.L,1.6L,-2.1L,.7L,-.4L})};break;
    case 11: g.curveDegree=5;g.referenceArea=.088229807652097;
      g.parts={Part({-1.3L,.7L,1.2L,-.7L,-.5L,-1.1L},{0.L,1.L})};break;
    case 12: g.curveDegree=5;g.referenceArea=.088229807652097;
      g.parts={Part({-.3L,.7L,1.2L,-.7L,-.5L,-1.1L},{0.L,-1.L})};break;
    case 13: g.curveDegree=5;g.mainCase=3;g.referenceArea=.607421467040005;
      g.parts={Part({-.5L,-.7L,1.5L,-3.4L,4.5L,.8L},{0.L,1.L})};break;
    case 14: g.curveDegree=3;g.referenceArea=.763596079675832;
      g.parts={Part({.7L,-1.L,-3.1L,3.L},{0.L,1.L})};break;
    case 15: g.curveDegree=2;g.referenceArea=.424193452759087;
      g.parts={Part({1.7L,-4.5L,3.1L},{0.L,-1.L})};break;
    case 16: g.curveDegree=4;g.referenceArea=.680254721840308;
      g.parts={Part({.7L,-1.3L,3.1L,-2.4L,1.L},{0.L,-1.L})};break;
    case 17: g.curveDegree=5;g.referenceArea=.598099177308039;
      g.parts={Part({-.3L,1.L},{0.L,-1.3L,2.1L,-2.7L,1.L,2.5L})};break;
    case 18: g.curveDegree=2;g.referenceArea=.8843362938564083;
      // Authors' rational NURBS arc is exactly the circle of center
      // (1.1,0.1), radius 0.4. Integrate its exterior within [0,1]^2.
      g.parts={Part({1.21L,-2.2L,1.L},{-.15L,-.2L,1.L})};break;
    case 19: g.curveDegree=2;g.referenceArea=.713524171452692;
      g.parts={Part({1.5L,-1.L},{0.L,-2.6L,1.7L})};break;
    case 20: g.curveDegree=3;g.referenceArea=.588688397937261;
      g.parts={Part({.1L,1.L},{0.L,-1.7L,.5L,.7L})};break;
    case 21: g.curveDegree=3;g.referenceArea=.804603612212857;
      g.parts={Part({.2L,1.L},{0.L,-1.L,-2.1L,3.4L})};break;
    case 22: g.curveDegree=4;g.referenceArea=.944556961342028;
      g.parts={Part({1.2L,-1.L},{0.L,-1.6L,2.1L,-.7L,.4L})};break;
    case 23: g.curveDegree=5;g.referenceArea=.911770192347903;
      g.parts={Part({1.3L,-.7L,-1.2L,.7L,.5L,1.1L},{0.L,-1.L})};break;
    case 24: g.curveDegree=5;g.referenceArea=.911770192347903;
      g.parts={Part({.3L,-.7L,-1.2L,.7L,.5L,1.1L},{0.L,1.L})};break;
    case 25: g.curveDegree=2;g.referenceArea=.0057329597491502;
      g.parts={Part({-2.34L,3.23L,-1.9L},{0.L,1.L})};break;
    case 26: g.curveDegree=2;g.referenceArea=.9942670402508498;
      g.parts={Part({2.34L,-3.23L,1.9L},{0.L,-1.L})};break;
    default:throw std::invalid_argument("Suite geometry must be 1..26");
  }
  return g;
}
Geometry MainCase(unsigned id) {
  unsigned suiteId=0;
  const char* name="";
  if(id==1) {suiteId=1;name="linear";}
  else if(id==2) {suiteId=2;name="quadratic";}
  else if(id==3) {suiteId=13;name="quintic";}
  else throw std::invalid_argument("Main case must be 1..3");
  Geometry g=Suite(suiteId);
  g.id=id; // Main study IDs are 1,2,3; suite ID 13 denotes quintic.
  g.name=name;
  return g;
}
long double Eval(const std::vector<long double>& a,long double x) {
  long double y=0;
  for(auto i=a.rbegin();i!=a.rend();++i)y=y*x+*i;
  return y;
}
Geometry Moving(long double shift) {
  Geometry g;g.id=0;g.curveDegree=2;g.name="moving-parabolic-lens";
  g.referenceArea=2.0*1.155625408802561/16.0;
  g.moving=true;g.intersection=true;
  // x-range of the parabola intersections; the lens translates in y.
  const long double a=(26.L-std::sqrt(52.L))/41.6L;
  const long double b=(26.L+std::sqrt(52.L))/41.6L;
  // U(x)=-10.4*x^2+13*x-3.25 + shift, L(x) = +10.4*x^2-13*x+4.25 + shift.
  // Original script: shift = y_shift - d, y_shift=0.0625, d=0..0.25.
  g.parts={Part({-3.25L+shift,13.L,-10.4L},{0.L,-1.L},a,b,+1),
           Part({ 4.25L+shift,-13.L,10.4L},{0.L,-1.L},a,b,-1)};
  return g;
}
double Psi(const Geometry& g,double x,double y) {
  if(g.intersection) {
    const auto& upper=g.parts.at(0);const auto& lower=g.parts.at(1);
    const long double U=Eval(upper.xPoly,x)+Eval(upper.yPoly,y);
    const long double L=Eval(lower.xPoly,x)+Eval(lower.yPoly,y);
    return static_cast<double>(std::min(U,-L));
  }
  for(const auto& part:g.parts) {
    if(x>=part.xLow-1.e-14 && x<=part.xHigh+1.e-14) {
      long double value=Eval(part.xPoly,x)+Eval(part.yPoly,y);
      if(g.id==6) value=std::min(value,static_cast<long double>(x)-.25L);
      return static_cast<double>(value);
    }
  }
  // Outside a clipped domain: negative level set.
  return -1.;
}

// Main-case scaled monomials (Eq. 34) and an independent reference integral.
double Function(const Geometry& g,unsigned q,double x,double y) {
  if(q==0)return 1.;
  const double normalizedX=g.mainCase==1 ? (x-.2)/.8 : x/.95;
  return std::pow(normalizedX,static_cast<int>(q))*
         std::pow(y,static_cast<int>(q));
}
struct Node {long double x,w;};
std::vector<Node> ReferenceGauss() {
  const unsigned n=64;
  const long double pi=std::acos(-1.L);
  std::vector<Node> nodes(n);
  for(unsigned i=0;i<n/2;++i) {
    long double z=std::cos(pi*(i+.75L)/(n+.5L));
    for(unsigned k=0;k<100;++k) {
      long double p0=1,p1=z;
      for(unsigned j=2;j<=n;++j) {
        const long double p2=((2*j-1)*z*p1-(j-1)*p0)/j;
        p0=p1;p1=p2;
      }
      const long double dp=n*(z*p1-p0)/(z*z-1.L);
      const long double dz=p1/dp;
      z-=dz;
      if(std::abs(dz)<8*std::numeric_limits<long double>::epsilon())break;
    }
    long double p0=1,p1=z;
    for(unsigned j=2;j<=n;++j) {
      const long double p2=((2*j-1)*z*p1-(j-1)*p0)/j;
      p0=p1;p1=p2;
    }
    const long double dp=n*(z*p1-p0)/(z*z-1.L);
    const long double weight=1.L/((1.L-z*z)*dp*dp);
    nodes[i]={.5L*(1.L-z),weight};
    nodes[n-1-i]={.5L*(1.L+z),weight};
  }
  return nodes;
}
long double Reference(const Geometry& g,unsigned q) {
  if(q==0)return g.referenceArea;
  if(g.mainCase<1 || g.mainCase>3)
    throw std::invalid_argument("Polynomial integrands are only defined for the 3 main cases");
  static const auto nodes=ReferenceGauss();
  const long double qq=static_cast<long double>(q+1);
  const long double cutoff=.750417596651957L;
  long double integral=0;
  for(const auto& node:nodes) {
    const long double t=node.x;
    long double value=0;
    if(g.mainCase==1) {
      value=std::pow(t,static_cast<int>(q))*.8L/qq*
        (1.L-std::pow(.25L*t,static_cast<int>(q+1)));
    }else if(g.mainCase==2) {
      const long double xmax=.8L+.3L*t-.6L*t*t;
      value=std::pow(t,static_cast<int>(q))*.95L/qq*
        std::pow(xmax/.95L,static_cast<int>(q+1));
    }else {
      const long double x0=cutoff*t,x1=cutoff+(1.L-cutoff)*t;
      const auto& p=g.parts[0];
      const long double h0=-Eval(p.xPoly,x0); // boundary y=-px(x)
      value=cutoff*std::pow(x0/.95L,static_cast<int>(q))/qq*
              (1.L-std::pow(h0,static_cast<int>(q+1)))+
            (1.L-cutoff)*std::pow(x1/.95L,static_cast<int>(q))/qq;
    }
    integral+=node.w*value;
  }
  return integral;
}
void SetPsi(MultiLevelSolution& ml,const Geometry& g) {
  Mesh& mesh=*ml.GetMultilevelMesh()->GetLevel(0);
  Solution& sol=*ml.GetLevel(0);
  const unsigned rank=mesh.processor_id();
  const unsigned index=ml.GetIndex("Psi");
  const unsigned type=ml.GetSolutionType("Psi");
  auto* field=sol._Sol[index];field->zero();
  for(unsigned e=mesh._elementOffset[rank];e<mesh._elementOffset[rank+1];++e) {
    for(unsigned j=0;j<mesh.GetElementDofNumber(e,type);++j) {
      const unsigned coordDof=mesh.GetSolutionDof(j,e,2);
      const double x=(*mesh._topology->_Sol[0])(coordDof);
      const double y=(*mesh._topology->_Sol[1])(coordDof);
      field->set(mesh.GetSolutionDof(j,e,type),Psi(g,x,y));
    }
  }
  field->close();
}
long double Binomial(unsigned n,unsigned k) {
  long double r=1.;for(unsigned j=1;j<=k;++j)r=r*(n+1-j)/j;return r;
}
// Power coefficients of Psi(x0+hx*u,y0+hy*v), with u,v in [0,1].
std::vector<long double> ExactCoefficients(const PolyPart& part,
                         long double x0,long double y0,long double hx,long double hy,
                         unsigned& nx,unsigned& ny) {
  nx=static_cast<unsigned>(part.xPoly.size());
  ny=static_cast<unsigned>(part.yPoly.size());
  std::vector<long double> c(nx*ny,0.L);
  for(unsigned m=0;m<nx;++m)for(unsigned k=0;k<=m;++k)
    c[k*ny]+=part.xPoly[m]*Binomial(m,k)*
      std::pow(x0,static_cast<int>(m-k))*std::pow(hx,static_cast<int>(k));
  for(unsigned m=0;m<ny;++m)for(unsigned k=0;k<=m;++k)
    c[k]+=part.yPoly[m]*Binomial(m,k)*
      std::pow(y0,static_cast<int>(m-k))*std::pow(hy,static_cast<int>(k));
  return c;
}
struct Counts {
  unsigned long long all=0,cut=0,full=0,qpCut=0,qpFull=0,qpSurface=0;
  unsigned long long total()const{return qpCut+qpFull;}
};
void Reduce(Counts& x) {
  unsigned long long a[]={x.all,x.cut,x.full,x.qpCut,x.qpFull,x.qpSurface};
  MPI_Allreduce(MPI_IN_PLACE,a,6,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
  x.all=a[0];x.cut=a[1];x.full=a[2];x.qpCut=a[3];x.qpFull=a[4];x.qpSurface=a[5];
}
struct Result {double value=0;Counts counts;unsigned nominal=0;};
struct Bounds {double x0,y0,hx,hy;};
Bounds ElementBounds(Mesh& mesh,unsigned e) {
  const unsigned v0=mesh.GetSolutionDof(0,e,2);
  const unsigned v2=mesh.GetSolutionDof(2,e,2);
  const double x0=(*mesh._topology->_Sol[0])(v0);
  const double y0=(*mesh._topology->_Sol[1])(v0);
  const double hx=(*mesh._topology->_Sol[0])(v2)-x0;
  const double hy=(*mesh._topology->_Sol[1])(v2)-y0;
  if(!(hx>0 && hy>0))
    throw std::runtime_error("Requires axis-aligned QUAD9 background cells");
  return {x0,y0,hx,hy};
}
void ReadPsi(MultiLevelSolution& ml,unsigned e,std::vector<double>& psi) {
  Mesh& mesh=*ml.GetMultilevelMesh()->GetLevel(0);
  Solution& sol=*ml.GetLevel(0);
  const unsigned idx=ml.GetIndex("Psi");
  psi.resize(mesh.GetElementDofNumber(e,2));
  for(unsigned j=0;j<psi.size();++j)
    psi[j]=(*sol._Sol[idx])(mesh.GetSolutionDof(j,e,2));
}
void WritePoint(std::ostream* out,const std::string& method,const Geometry& g,
  unsigned ne,unsigned nqp,unsigned elem,double x,double y,double w,bool cut) {
  if(!out)return;
  *out<<method<<','<<g.id<<','<<ne<<','<<nqp<<','<<elem<<','
      <<std::setprecision(17)<<x<<','<<y<<','<<w<<','<<(cut?1:0)<<'\n';
}
Result RunMoment(MultiLevelSolution& ml,const Geometry& g,unsigned power,unsigned nqp,
                 std::ostream* points=nullptr) {
  IntegrationOptions opts;
  opts.mode=IntegrationMode::MomentFitting;
  opts.momentDegree=2*nqp-1; // Paper: kq_requested=2*nQP,set-1.
  Result r;
  r.value=Integrate(ml,"Psi",[&](const std::vector<double>& xy){
    return Function(g,power,xy[0],xy[1]);},opts).omega1;
  Mesh& mesh=*ml.GetMultilevelMesh()->GetLevel(0);
  const unsigned rank=mesh.processor_id();
  Fem fem(2*opts.momentDegree,2);
  for(unsigned e=mesh._elementOffset[rank];e<mesh._elementOffset[rank+1];++e) {
    ++r.counts.all;
    const unsigned geom=mesh.GetElementType(e);
    const unsigned ngauss=fem.GetFiniteElement(geom,2)->GetGaussPointNumber();
    r.nominal=ngauss;
    std::vector<double> psi;ReadPsi(ml,e,psi);
    const bool cut=IsCutElementQuadratic(psi,geom);
    if(cut) {++r.counts.cut;r.counts.qpCut+=ngauss;}
    else if(psi[0]>0.) {++r.counts.full;r.counts.qpFull+=ngauss;}
    if(!points || (!cut && psi[0]<=0.))continue;
    const Bounds b=ElementBounds(mesh,e);
    auto& rule=ex44_integration_detail::GetMomentRule(geom,opts.momentDegree);
    std::vector<double> volume,surface;
    if(cut) {
      std::vector<std::vector<double>> xyz(2,std::vector<double>(psi.size()));
      for(unsigned j=0;j<psi.size();++j) {
        const unsigned n=mesh.GetSolutionDof(j,e,2);
        xyz[0][j]=(*mesh._topology->_Sol[0])(n);
        xyz[1][j]=(*mesh._topology->_Sol[1])(n);
      }
      rule.GetPhysicalSurfaceWeightsFromNodalGeometry(psi,xyz,volume,surface,
        opts.absoluteTolerance,opts.relativeTolerance);
    }
    const double* gw=rule.GetGaussWeightPointer();
    const double* gx=rule.GetGaussCoordinatePointer(0);
    const double* gy=rule.GetGaussCoordinatePointer(1);
    for(unsigned j=0;j<ngauss;++j) {
      const double x=b.x0+.5*b.hx*(1+gx[j]);
      const double y=b.y0+.5*b.hy*(1+gy[j]);
      const double v=(cut?volume[j]:1.)*gw[j]*b.hx*b.hy*.25;
      WritePoint(points,"femus-moment",g,0,nqp,e,x,y,v,cut);
    }
  }
  Reduce(r.counts);
  return r;
}
Result RunAlgoim(MultiLevelSolution& ml,const Geometry& g,unsigned power,unsigned nqp,
                unsigned subdiv,std::ostream* points=nullptr) {
  Mesh& mesh=*ml.GetMultilevelMesh()->GetLevel(0);
  const unsigned rank=mesh.processor_id();
  Result r;r.nominal=nqp*nqp;
  long double integral=0;
  for(unsigned e=mesh._elementOffset[rank];e<mesh._elementOffset[rank+1];++e) {
    ++r.counts.all;
    const Bounds b=ElementBounds(mesh,e);
    const long double xa=b.x0,xb=b.x0+b.hx;
    unsigned long long volumeCount=0,surfaceCount=0;
    bool partCut=false;
    long double elementIntegral=0;
    long double elementArea=0;
    for(const auto& part:g.parts) {
      const long double a=std::max(xa,part.xLow);
      const long double z=std::min(xb,part.xHigh);
      if(z<=a)continue;
      const long double width=z-a;
      if(width < b.hx-1.e-14L)partCut=true;
      unsigned nx=0,ny=0;
      const auto powerCoeffs=ExactCoefficients(part,a,b.y0,width,b.hy,nx,ny);
      FemusAlgoimDirect::Rule rule;
      FemusAlgoimDirect::detail::BuildTensor<2>(powerCoeffs,
        algoim::uvector<int,2>(static_cast<int>(nx),static_cast<int>(ny)),
        0,0,0,nqp,subdiv,rule);
      volumeCount+=rule.volume.size();
      surfaceCount+=rule.surfaceDelta.size();
      if(!rule.surfaceDelta.empty())partCut=true;
      for(const auto& qp:rule.volume) {
        const double x=static_cast<double>(a+.5L*width*(qp.xi[0]+1.L));
        const double y=b.y0+.5*b.hy*(qp.xi[1]+1.);
        const double physicalW=static_cast<double>(part.sign*width*b.hy*.25L)*qp.weight;
        elementIntegral+=static_cast<long double>(physicalW)*Function(g,power,x,y);
        elementArea+=physicalW;
        WritePoint(points,g.intersection?"algoim-signed-difference":"algoim-exact",
          g,0,nqp,e,x,y,physicalW,partCut || g.intersection);
      }
    }
    // For the signed-difference representation, a cell may have two full
    // but cancelling rules; all evaluations are still counted.
    integral+=elementIntegral;
    if(volumeCount) {
      const long double cellArea=static_cast<long double>(b.hx)*b.hy;
      const bool genuinelyCut=partCut ||
          (elementArea>1.e-12L*cellArea && elementArea<(1.L-1.e-12L)*cellArea);
      if(genuinelyCut) {++r.counts.cut;r.counts.qpCut+=volumeCount;}
      else if(elementArea>1.e-12L*cellArea) {++r.counts.full;r.counts.qpFull+=volumeCount;}
      // Some signed-difference rules cancel on inactive cells: their
      // function evaluations are counted separately only if nonzero.
    }
    r.counts.qpSurface+=surfaceCount;
  }
  const double local=static_cast<double>(integral);
  MPI_Allreduce(&local,&r.value,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  Reduce(r.counts);
  return r;
}
std::string MethodName(const Geometry& g,bool isAlgoim) {
  if(!isAlgoim)return "femus-moment-Q2";
  if(g.intersection)return "algoim-signed-difference";
  if(g.parts.size()>1 || g.id==6)return "algoim-subcell-polynomial";
  return "algoim-exact";
}
std::string GeometryName(const Geometry& g,bool isAlgoim) {
  if(isAlgoim) {
    if(g.intersection)return "two-exact-polynomials-difference";
    if(g.parts.size()>1 || g.id==6)return "exact-polynomial-subcells";
    return "exact-implicit-polynomial";
  }
  return g.curveDegree<=2 && !g.intersection && g.parts.size()==1 && g.id!=6?
    "Q2-exact-levelset":"Q2-interpolated-levelset";
}
std::string SafeCSV(const std::string& text) {
  std::string out=text;
  for(char& c:out)if(c==',' || c=='\n' || c=='\r' || c=='\t')c='_';
  return out;
}
void Row(std::ostream& f,const Options& o,const std::string& study,
  const std::string& figures,const Geometry& g,unsigned monomial,unsigned ne,unsigned nqp,
  int step,double shift,const std::string& method,const std::string& geometry,
  const std::string& status,const std::string& detail,const Result* r,double exact) {
  const double nan=std::numeric_limits<double>::quiet_NaN();
  const double val=r?r->value:nan;
  const double absErr=r?std::abs(val-exact):nan;
  const double relErr=r?(exact?absErr/std::abs(exact):nan):nan;
  f<<study<<','<<figures<<','<<g.id<<','<<g.name<<','<<g.curveDegree<<','
   <<monomial<<','<<(2*monomial)<<','<<ne<<','<<ne*ne<<','
   <<nqp<<','<<(2*nqp-1)<<','<<step<<','
   <<std::setprecision(17)<<shift<<','<<method<<','<<geometry<<','<<status<<','<<detail<<','
   <<val<<','<<exact<<','<<absErr<<','<<relErr<<','
   <<(r?r->counts.all:0)<<','<<(r?r->counts.cut:0)<<','<<(r?r->counts.full:0)<<','
   <<(r?r->nominal:0)<<','<<(r?r->counts.qpCut:0)<<','
   <<(r?r->counts.qpFull:0)<<','<<(r?r->counts.total():0)<<','
   <<(r?r->counts.qpSurface:0)<<",2"<<'\n';
  if(!o.quiet) {
    std::cout<<study<<" geometry="<<g.id<<" q="<<monomial<<" ne="<<ne
     <<" nQP,set="<<nqp<<" method="<<method<<" status="<<status;
    if(r)std::cout<<" error="<<std::scientific<<relErr<<" cutQP="<<r->counts.qpCut;
    std::cout<<'\n';
  }
}
void RunOne(std::ostream& csv,const Options& o,const std::string& study,
            const std::string& figures,const Geometry& g,unsigned q,unsigned ne,
            unsigned nqp,int rank,int step=-1,double shift=0,
            std::ostream* points=nullptr) {
  // nQP,set is still the paper's requested order, even for skipped methods.
  if(nqp<1 || nqp>100)throw std::invalid_argument("Algoim supports nQP,set=1..100");
  const double exact=static_cast<double>(Reference(g,q));
  const bool wantMoment=o.mode!="algoim";
  const bool wantAlgoim=o.mode!="moment";
  const bool momentSupported=(4*nqp-2)<=37; // Fem allows gaussOrder<=37.
  const bool runMoment=wantMoment && momentSupported;
  MultiLevelMesh mesh;
  mesh.GenerateCoarseBoxMesh(ne,ne,0,0.,1.,0.,1.,0.,0.,QUAD9,"fifth");
  MultiLevelSolution ml(&mesh);
  ml.AddSolution("Psi",LAGRANGE,SECOND,0,false);
  ml.Initialize("All");
  SetPsi(ml,g);
  if(runMoment) {
    try {
      const Result r=RunMoment(ml,g,q,nqp,points);
      if(rank==0)Row(csv,o,study,figures,g,q,ne,nqp,step,shift,
        MethodName(g,false),GeometryName(g,false),"computed","",&r,exact);
    }catch(const std::exception& e) {
      // The benchmark is intended to run serially: a failed quadrature
      // configuration must not prevent later geometries from being tested.
      if(rank==0)Row(csv,o,study,figures,g,q,ne,nqp,step,shift,
        MethodName(g,false),GeometryName(g,false),"failed",
        SafeCSV(e.what()),nullptr,exact);
    }
  }else if(wantMoment && rank==0) {
    Row(csv,o,study,figures,g,q,ne,nqp,step,shift,
      MethodName(g,false),GeometryName(g,false),"skipped",
      "FEMuS_Fem_gaussOrder_limit_37",nullptr,exact);
  }
  if(wantAlgoim) {
    try {
      const Result r=RunAlgoim(ml,g,q,nqp,o.algoimSubdivisions,points);
      if(rank==0)Row(csv,o,study,figures,g,q,ne,nqp,step,shift,
        MethodName(g,true),GeometryName(g,true),"computed","",&r,exact);
    }catch(const std::exception& e) {
      if(rank==0)Row(csv,o,study,figures,g,q,ne,nqp,step,shift,
        MethodName(g,true),GeometryName(g,true),"failed",
        SafeCSV(e.what()),nullptr,exact);
    }
  }
}
std::vector<unsigned> SelectedCases(const Options& o) {
  if(o.caseId)return {o.caseId};
  return {1,2,3};
}
bool WantN(const Options& o,unsigned n,unsigned defaultMax) {
  return n>=o.nqpMin && n<=(o.nqpMax?o.nqpMax:defaultMax);
}
void Area(std::ostream& csv,const Options& o,int rank) {
  for(unsigned k:SelectedCases(o)) {
    const auto g=MainCase(k);
    for(unsigned n=1;n<=4;++n)
      if(WantN(o,n,4))RunOne(csv,o,"area","5_9_18_20_23_28_31",g,0,1,n,rank);
  }
}
void AreaHRef(std::ostream& csv,const Options& o,int rank) {
  for(unsigned k:SelectedCases(o)) {
    const auto g=MainCase(k);
    std::vector<unsigned> choices;
    if(o.hrefNqp) choices.push_back(o.hrefNqp);
    else choices={1,2,3,4};
    for(const unsigned n:choices) {
      if(!WantN(o,n,4))continue;
      for(unsigned l=0;l<=o.levels;++l)
        RunOne(csv,o,"area-href","6_10_13_15_24",g,0,1u<<l,n,rank);
    }
  }
}
void Monomial(std::ostream& csv,const Options& o,int rank) {
  for(unsigned k:SelectedCases(o)) {
    const auto g=MainCase(k);
    for(unsigned q=o.qMin;q<=o.qMax;++q) {
      // Companion scripts use p*(q+1)+1, INCLUDING the final extra point.
      const unsigned paperMax=g.curveDegree*(q+1)+1;
      for(unsigned n=o.nqpMin;n<=std::min(paperMax,o.nqpMax?o.nqpMax:paperMax);++n)
        RunOne(csv,o,"monomial","7_11_16_19_21_25_29_32_33_38",g,q,1,n,rank);
    }
  }
}
void MonomialHRef(std::ostream& csv,const Options& o,int rank) {
  // Fig. 39 and Fig. 40 are specifically geometry 2 and q=6.
  // The companion repository also has separate q=6 h-ref scripts for
  // geometry 1 (nQP=4) and geometry 3 (nQP=6).
  std::vector<unsigned> c=o.caseId?std::vector<unsigned>{o.caseId}:std::vector<unsigned>{2};
  for(unsigned k:c) {
    const Geometry g=MainCase(k);
    const std::vector<unsigned> choices=o.hrefNqp?
      std::vector<unsigned>{o.hrefNqp}:
      (k==2?std::vector<unsigned>{3,6}:std::vector<unsigned>{k==1?4u:6u});
    for(unsigned n:choices) {
      if(!WantN(o,n,100))continue;
      for(unsigned l=0;l<=o.levels;++l)
        RunOne(csv,o,"monomial-href","39_40",g,6,1u<<l,n,rank);
    }
  }
}
void Layout(std::ostream& csv,std::ostream& points,const Options& o,int rank) {
  const Geometry g=MainCase(3); // Fig. 34: 2x2 and nQP,set=3.
  if(WantN(o,3,3))RunOne(csv,o,"layout","34",g,0,2,3,rank,-1,0.,&points);
}
void Suite26(std::ostream& csv,const Options& o,int rank) {
  for(unsigned id=1;id<=26;++id) {
    if(o.suiteId && id!=o.suiteId)continue;
    const Geometry g=Suite(id);
    const unsigned n=static_cast<unsigned>((g.curveDegree+2)/2); // ceil((p+1)/2)
    if(!WantN(o,n,100))continue;
    // Paper Fig. 35: ne=1 and 8 per direction. FEMuS/Algoim are NOT
    // the paper's BoSSS/mlhp/QuESo implementations.
    for(const unsigned ne: {1u,8u})
      RunOne(csv,o,"suite26","35",g,0,ne,n,rank);
  }
}
void MovingTest(std::ostream& csv,const Options& o,int rank) {
  const unsigned n=o.hrefNqp?o.hrefNqp:2; // Fig.37, Table3
  if(!WantN(o,n,100))return;
  const unsigned last=o.stepEnd?o.stepEnd:o.steps-1;
  for(unsigned step=o.stepStart;step<=last;++step) {
    const double d=.25*static_cast<double>(step)/(o.steps-1);
    const double offset=.0625-d;
    const Geometry g=Moving(offset);
    RunOne(csv,o,"moving","37_Table3",g,0,8,n,rank,
           static_cast<int>(step),d);
  }
}

// ----- Three-dimensional studies from Sect. 4.2 of Loibl et al. (2026) -----

enum class Geometry3DKind { Ellipsoid, Torus, Cylinder };
struct Geometry3D {
  Geometry3DKind kind;
  std::string name;
  unsigned degree;
  double lower[3], upper[3];
  double shift=0; // paper displacement d, used by the moving ellipsoid
  unsigned caseId=0;
  double exactVolume=0;
};
Geometry3D Make3D(Geometry3DKind kind,double shift=0) {
  const double pi=std::acos(-1.0);
  if(kind==Geometry3DKind::Ellipsoid) {
    return {kind,"shifted-ellipsoid",2,{-1.1,-1.1,-1.1},
      {1.1,1.1,1.1},shift,1,4.0*pi*(1.0*0.5*(1.0/3.0))/3.0};
  }
  if(kind==Geometry3DKind::Torus) {
    return {kind,"torus",4,{-1.4,-1.4,-1.4},
      {1.4,1.4,1.4},0.0,2,2.0*pi*pi*1.0*(0.2*0.2)};
  }
  return {kind,"elliptic-cylinder",2,{-0.6,-0.6,0.0},
    {0.65,0.65,1.0},0.0,3,pi*0.5*0.25*1.0};
}
// Return positive level set inside the original paper geometry.
double Psi3D(const Geometry3D& g,double x,double y,double z) {
  if(g.kind==Geometry3DKind::Ellipsoid) {
    const double yc=-2.2/8.0+g.shift;
    return 1.0 - x*x - 4.0*(y-yc)*(y-yc) - 9.0*z*z;
  }
  if(g.kind==Geometry3DKind::Cylinder)
    return 1.0 - 4.0*x*x - 16.0*y*y;
  // Equivalent exact quartic implicit polynomial for a torus of
  // major radius R=1, minor radius r=0.2.
  const double rho2=x*x+y*y;
  const double t=rho2+z*z+0.96; // R^2-r^2
  return 4.0*rho2-t*t;
}
double Function3D(unsigned q,double x,double y,double z) {
  if(q==0)return 1.0;
  // Exact bounding box of the NURBS cylinder: x=[-0.5,0.5],
  // y=[-0.25,0.25], z=[0,1]. The paper uses these normalizations.
  return std::pow(x+0.5,static_cast<int>(q)) *
         std::pow(2.0*y+0.5,static_cast<int>(q)) *
         std::pow(z,static_cast<int>(q));
}
long double Reference3D(const Geometry3D& g,unsigned q) {
  if(q==0)return g.exactVolume;
  if(g.kind!=Geometry3DKind::Cylinder)
    throw std::invalid_argument("3D monomials are defined only for the cylinder");
  // Set x=0.5*u, y=0.25*v, where u^2+v^2<=1.
  // The disk moments of u^(2i)v^(2j) are
  // Gamma(i+1/2)*Gamma(j+1/2)/Gamma(i+j+2).
  long double sum=0;
  for(unsigned i=0;2*i<=q;++i)
    for(unsigned j=0;2*j<=q;++j)
      sum+=Binomial(q,2*i)*Binomial(q,2*j)*
        std::tgamma(static_cast<long double>(i)+.5L)*
        std::tgamma(static_cast<long double>(j)+.5L)/
        std::tgamma(static_cast<long double>(i+j)+2.0L);
  return (0.5L*0.25L)*sum /
    (std::pow(4.0L,static_cast<int>(q))*static_cast<long double>(q+1));
}
void SetPsi3D(MultiLevelSolution& ml,const Geometry3D& g) {
  Mesh& mesh=*ml.GetMultilevelMesh()->GetLevel(0);
  Solution& sol=*ml.GetLevel(0);
  const unsigned rank=mesh.processor_id();
  const unsigned idx=ml.GetIndex("Psi"),type=ml.GetSolutionType("Psi");
  auto* field=sol._Sol[idx];field->zero();
  for(unsigned e=mesh._elementOffset[rank];e<mesh._elementOffset[rank+1];++e)
    for(unsigned j=0;j<mesh.GetElementDofNumber(e,type);++j) {
      const unsigned d=mesh.GetSolutionDof(j,e,2);
      const double x=(*mesh._topology->_Sol[0])(d);
      const double y=(*mesh._topology->_Sol[1])(d);
      const double z=(*mesh._topology->_Sol[2])(d);
      field->set(mesh.GetSolutionDof(j,e,type),Psi3D(g,x,y,z));
    }
  field->close();
}
struct Bounds3D {double lo[3],hi[3],length[3];};
Bounds3D ElementBounds3D(Mesh& mesh,unsigned e) {
  Bounds3D b;
  const unsigned nd=mesh.GetElementDofNumber(e,2);
  for(unsigned d=0;d<3;++d){b.lo[d]=std::numeric_limits<double>::infinity();b.hi[d]=-b.lo[d];}
  for(unsigned j=0;j<nd;++j) {
    const unsigned dof=mesh.GetSolutionDof(j,e,2);
    for(unsigned d=0;d<3;++d) {
      const double x=(*mesh._topology->_Sol[d])(dof);
      b.lo[d]=std::min(b.lo[d],x);
      b.hi[d]=std::max(b.hi[d],x);
    }
  }
  for(unsigned d=0;d<3;++d) {
    b.length[d]=b.hi[d]-b.lo[d];
    if(!(b.length[d]>0))throw std::runtime_error("Requires axis-aligned HEX27 cells");
  }
  return b;
}
// Conservative interval tests for empty/full boxes to avoid expensive
// polynomial quadrature on cells far from the boundary.
// Returns -1 outside, +1 full, 0 cut or inconclusive.
static double NearSquare(double lo,double hi,double center) {
  return (lo>center)?std::pow(lo-center,2):(hi<center)?std::pow(hi-center,2):0.0;
}
static double FarSquare(double lo,double hi,double center) {
  return std::max(std::pow(lo-center,2),std::pow(hi-center,2));
}
int ClassifyExact3D(const Geometry3D& g,const Bounds3D& b) {
  double mn=0,mx=0,threshold=1;
  if(g.kind==Geometry3DKind::Ellipsoid) {
    const double yc=-2.2/8.0+g.shift;
    const double centers[3]={0,yc,0},coef[3]={1,4,9};
    for(unsigned d=0;d<3;++d) {
      mn+=coef[d]*NearSquare(b.lo[d],b.hi[d],centers[d]);
      mx+=coef[d]*FarSquare(b.lo[d],b.hi[d],centers[d]);
    }
  } else if(g.kind==Geometry3DKind::Cylinder) {
    mn=4*NearSquare(b.lo[0],b.hi[0],0)+16*NearSquare(b.lo[1],b.hi[1],0);
    mx=4*FarSquare(b.lo[0],b.hi[0],0)+16*FarSquare(b.lo[1],b.hi[1],0);
  } else {
    const double rmin=std::sqrt(NearSquare(b.lo[0],b.hi[0],0)+
                                NearSquare(b.lo[1],b.hi[1],0));
    const double rmax=std::sqrt(FarSquare(b.lo[0],b.hi[0],0)+
                                FarSquare(b.lo[1],b.hi[1],0));
    const double dmin=(1<rmin)?rmin-1:(1>rmax)?1-rmax:0;
    const double dmax=std::max(std::abs(rmin-1),std::abs(rmax-1));
    mn=dmin*dmin+NearSquare(b.lo[2],b.hi[2],0);
    mx=dmax*dmax+FarSquare(b.lo[2],b.hi[2],0);
    threshold=0.04; // tube radius squared
  }
  if(mn>threshold+1.e-12)return -1;
  if(mx<threshold-1.e-12)return 1;
  return 0;
}
struct Term3D {unsigned i,j,k;long double c;};
std::vector<Term3D> PowerTerms3D(const Geometry3D& g) {
  if(g.kind==Geometry3DKind::Cylinder)
    return {{0,0,0,1},{2,0,0,-4},{0,2,0,-16}};
  if(g.kind==Geometry3DKind::Ellipsoid) {
    const long double yc=-2.2L/8.0L+g.shift;
    return {{0,0,0,1.0L-4.0L*yc*yc},{2,0,0,-1},
      {0,1,0,8.0L*yc},{0,2,0,-4},{0,0,2,-9}};
  }
  // 4(x^2+y^2) - (x^2+y^2+z^2+0.96)^2.
  return {{0,0,0,-.9216L},{2,0,0,2.08L},{0,2,0,2.08L},
    {0,0,2,-1.92L},{4,0,0,-1},{0,4,0,-1},{0,0,4,-1},
    {2,2,0,-2},{2,0,2,-2},{0,2,2,-2}};
}
std::vector<long double> ExactCoefficients3D(const Geometry3D& g,
                          const Bounds3D& b,unsigned& nx,unsigned& ny,unsigned& nz) {
  nx=ny=nz=(g.kind==Geometry3DKind::Torus)?5u:3u;
  // One extra zero linear z coefficient avoids degeneracy in the
  // general 3D Bernstein representation of a z-independent cylinder.
  auto result=std::vector<long double>(nx*ny*nz,0.L);
  for(const auto& term:PowerTerms3D(g)) {
    for(unsigned i=0;i<=term.i;++i)
      for(unsigned j=0;j<=term.j;++j)
        for(unsigned k=0;k<=term.k;++k) {
          const long double c=term.c*Binomial(term.i,i)*Binomial(term.j,j)*
            Binomial(term.k,k)*
            std::pow(static_cast<long double>(b.lo[0]),static_cast<int>(term.i-i))*
            std::pow(static_cast<long double>(b.lo[1]),static_cast<int>(term.j-j))*
            std::pow(static_cast<long double>(b.lo[2]),static_cast<int>(term.k-k))*
            std::pow(static_cast<long double>(b.length[0]),static_cast<int>(i))*
            std::pow(static_cast<long double>(b.length[1]),static_cast<int>(j))*
            std::pow(static_cast<long double>(b.length[2]),static_cast<int>(k));
          result[(i*ny+j)*nz+k]+=c;
        }
  }
  return result;
}
Result RunMoment3D(MultiLevelSolution& ml,unsigned power,unsigned nqp) {
  IntegrationOptions opts;
  opts.mode=IntegrationMode::MomentFitting;
  opts.momentDegree=2*nqp-1;
  Result r;
  r.value=Integrate(ml,"Psi",[&](const std::vector<double>& xyz) {
    return Function3D(power,xyz[0],xyz[1],xyz[2]);},opts).omega1;
  Mesh& mesh=*ml.GetMultilevelMesh()->GetLevel(0);
  const unsigned rank=mesh.processor_id();
  Fem fem(2*opts.momentDegree,3);
  for(unsigned e=mesh._elementOffset[rank];e<mesh._elementOffset[rank+1];++e) {
    ++r.counts.all;
    const unsigned geom=mesh.GetElementType(e);
    const unsigned ng=fem.GetFiniteElement(geom,2)->GetGaussPointNumber();
    r.nominal=ng;
    std::vector<double> psi;ReadPsi(ml,e,psi);
    if(IsCutElementQuadratic(psi,geom)) {++r.counts.cut;r.counts.qpCut+=ng;}
    else if(psi[0]>0) {++r.counts.full;r.counts.qpFull+=ng;}
  }
  Reduce(r.counts);
  return r;
}
Result RunAlgoim3D(MultiLevelSolution& ml,const Geometry3D& g,unsigned power,
                   unsigned nqp,unsigned subdiv) {
  Mesh& mesh=*ml.GetMultilevelMesh()->GetLevel(0);
  const unsigned rank=mesh.processor_id();
  Result r;r.nominal=nqp*nqp*nqp;
  long double sum=0;
  for(unsigned e=mesh._elementOffset[rank];e<mesh._elementOffset[rank+1];++e) {
    ++r.counts.all;
    const Bounds3D b=ElementBounds3D(mesh,e);
    const int classification=ClassifyExact3D(g,b);
    if(classification<0)continue;
    FemusAlgoimDirect::Rule rule;
    if(classification>0) {
      FemusAlgoimDirect::detail::AddFullBox<3>(
        algoim::uvector<double,3>(0.0),algoim::uvector<double,3>(1.0),
        nqp,0,rule);
    } else {
      unsigned nx,ny,nz;
      const auto c=ExactCoefficients3D(g,b,nx,ny,nz);
      FemusAlgoimDirect::detail::BuildTensor<3>(c,
        algoim::uvector<int,3>(static_cast<int>(nx),static_cast<int>(ny),
                               static_cast<int>(nz)),0,0,0,nqp,subdiv,rule);
    }
    if(classification==0) {++r.counts.cut;r.counts.qpCut+=rule.volume.size();}
    else {++r.counts.full;r.counts.qpFull+=rule.volume.size();}
    r.counts.qpSurface+=rule.surfaceDelta.size();
    const long double jac=static_cast<long double>(b.length[0])*
        b.length[1]*b.length[2]/8.0L;
    for(const auto& qp:rule.volume) {
      const double x=b.lo[0]+.5*b.length[0]*(qp.xi[0]+1);
      const double y=b.lo[1]+.5*b.length[1]*(qp.xi[1]+1);
      const double z=b.lo[2]+.5*b.length[2]*(qp.xi[2]+1);
      sum+=jac*qp.weight*Function3D(power,x,y,z);
    }
  }
  const double local=static_cast<double>(sum);
  MPI_Allreduce(&local,&r.value,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  Reduce(r.counts);
  return r;
}
void Row3D(std::ostream& csv,const Options& o,const std::string& study,
    const std::string& figure,const Geometry3D& g,unsigned q,unsigned ne,
    unsigned nqp,int step,double displacement,
    const std::string& method,const std::string& model,
    const std::string& status,const std::string& detail,
    const Result* r,long double reference) {
  const double nan=std::numeric_limits<double>::quiet_NaN();
  const double exact=static_cast<double>(reference);
  const double numerical=r?r->value:nan;
  const double err=r?std::abs(numerical-exact):nan;
  const double relative=r?(exact?err/std::abs(exact):nan):nan;
  csv<<study<<','<<figure<<','<<g.caseId<<','<<g.name<<','<<g.degree<<','
     <<q<<','<<3*q<<','<<ne<<','<<(static_cast<unsigned long long>(ne)*ne*ne)<<','
     <<nqp<<','<<(2*nqp-1)<<','<<step<<','<<std::setprecision(17)
     <<displacement<<','<<method<<','<<model<<','<<status<<','<<SafeCSV(detail)<<','
     <<numerical<<','<<exact<<','<<err<<','<<relative<<','
     <<(r?r->counts.all:0)<<','<<(r?r->counts.cut:0)<<','<<(r?r->counts.full:0)<<','
     <<(r?r->nominal:0)<<','<<(r?r->counts.qpCut:0)<<','
     <<(r?r->counts.qpFull:0)<<','<<(r?r->counts.total():0)<<','
     <<(r?r->counts.qpSurface:0)<<",3\n";
  if(!o.quiet) {
    std::cout<<study<<" q="<<q<<" ne="<<ne<<" nQP,set="<<nqp
      <<" method="<<method<<" status="<<status;
    if(r)std::cout<<" relerr="<<std::scientific<<relative
                  <<" cutQP="<<r->counts.qpCut;
    std::cout<<'\n';
  }
}
void RunOne3D(std::ostream& csv,const Options& o,const std::string& study,
    const std::string& figure,const Geometry3D& g,unsigned q,unsigned ne,
    unsigned nqp,int rank,int step=-1,double shift=0) {
  if(nqp<1 || nqp>100)throw std::invalid_argument("3D nQP must be 1..100");
  const long double exact=Reference3D(g,q);
  const bool moment=o.mode!="algoim",algoimOn=o.mode!="moment";
  const bool supported=4*nqp-2<=37;
  MultiLevelMesh mesh;
  mesh.GenerateCoarseBoxMesh(ne,ne,ne,g.lower[0],g.upper[0],
    g.lower[1],g.upper[1],g.lower[2],g.upper[2],HEX27,"fifth");
  MultiLevelSolution ml(&mesh);
  ml.AddSolution("Psi",LAGRANGE,SECOND,0,false);
  ml.Initialize("All");
  SetPsi3D(ml,g);
  if(moment && supported) {
    try {
      const Result r=RunMoment3D(ml,q,nqp);
      if(rank==0)Row3D(csv,o,study,figure,g,q,ne,nqp,step,shift,
        "femus-moment-HEX27",
        g.kind==Geometry3DKind::Torus?"Q2-interpolated-quartic-torus":"Q2-exact-quadratic",
        "computed","",&r,exact);
    } catch(const std::exception& e) {
      if(rank==0)Row3D(csv,o,study,figure,g,q,ne,nqp,step,shift,
        "femus-moment-HEX27","HEX27-Q2","failed",e.what(),nullptr,exact);
    }
  }else if(moment && rank==0) {
    Row3D(csv,o,study,figure,g,q,ne,nqp,step,shift,
      "femus-moment-HEX27","HEX27-Q2","skipped","FEMuS_Gauss_order_limit_37",nullptr,exact);
  }
  if(algoimOn) {
    try {
      const Result r=RunAlgoim3D(ml,g,q,nqp,o.algoimSubdivisions);
      if(rank==0)Row3D(csv,o,study,figure,g,q,ne,nqp,step,shift,
        "algoim-exact-3D","exact-implicit-polynomial","computed","",&r,exact);
    }catch(const std::exception& e) {
      if(rank==0)Row3D(csv,o,study,figure,g,q,ne,nqp,step,shift,
        "algoim-exact-3D","exact-implicit-polynomial","failed",e.what(),nullptr,exact);
    }
  }
}
void EllipsoidMoving3D(std::ostream& csv,const Options& o,int rank) {
  const unsigned n=o.hrefNqp?o.hrefNqp:2;
  if(!WantN(o,n,100))return;
  const unsigned last=o.stepEnd?o.stepEnd:o.steps-1;
  for(unsigned i=o.stepStart;i<=last;++i) {
    const double d=(2.0*2.2/8.0)*static_cast<double>(i)/(o.steps-1);
    RunOne3D(csv,o,"ellipsoid-moving","42_Table4",
      Make3D(Geometry3DKind::Ellipsoid,d),0,8,n,rank,
      static_cast<int>(i),d);
  }
}
void EllipsoidMeshMoving3D(std::ostream& csv,const Options& o,int rank) {
  // Paper's STL/QuESo variant: keep the ellipsoid fixed at y=-0.275,
  // translate its bounding-box mesh by -d (equivalent relative position).
  const unsigned n=o.hrefNqp?o.hrefNqp:1; // QuESo uses nQP,set=1.
  if(!WantN(o,n,100))return;
  const unsigned last=o.stepEnd?o.stepEnd:o.steps-1;
  for(unsigned i=o.stepStart;i<=last;++i) {
    const double d=(2.0*2.2/8.0)*static_cast<double>(i)/(o.steps-1);
    Geometry3D g=Make3D(Geometry3DKind::Ellipsoid);
    g.name="fixed-ellipsoid-moving-mesh";
    g.caseId=7;
    g.lower[1]-=d;
    g.upper[1]-=d;
    RunOne3D(csv,o,"ellipsoid-mesh-moving","42_Table4_QuESo_mesh_variant",
      g,0,8,n,rank,static_cast<int>(i),d);
  }
}
void TorusHRef3D(std::ostream& csv,const Options& o,int rank) {
  const unsigned n=o.hrefNqp?o.hrefNqp:2;
  if(!WantN(o,n,100))return;
  const auto g=Make3D(Geometry3DKind::Torus);
  for(unsigned l=0;l<=o.levels;++l)
    RunOne3D(csv,o,"torus-href","44",g,0,1u<<l,n,rank);
}
void CylinderMonomial3D(std::ostream& csv,const Options& o,int rank) {
  const auto g=Make3D(Geometry3DKind::Cylinder);
  // The published q=3 script sets nQP,set=1..4. Figure 46 additionally
  // presents q=5; apply the same requested order range for a consistent test.
  for(unsigned q:{3u,5u}) {
    if(q<o.qMin || q>o.qMax)continue;
    const unsigned nmax=o.nqpMax?std::min(o.nqpMax,4u):4u;
    for(unsigned n=o.nqpMin;n<=nmax;++n)
      RunOne3D(csv,o,"cylinder-monomial","46",g,q,8,n,rank);
  }
}
void CylinderHRef3D(std::ostream& csv,const Options& o,int rank) {
  if(o.qMin>5 || o.qMax<5)return;
  const unsigned n=o.hrefNqp?o.hrefNqp:3;
  if(!WantN(o,n,100))return;
  const auto g=Make3D(Geometry3DKind::Cylinder);
  for(unsigned l=0;l<=o.levels;++l)
    RunOne3D(csv,o,"cylinder-href","47",g,5,1u<<l,n,rank);
}
} // namespace

int main(int argc,char** argv) {
  FemusInit init(argc,argv,MPI_COMM_WORLD);
  int rank=0,commSize=1;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&commSize);
  try {
    const Options o=Parse(argc,argv);
    if(commSize!=1)
      throw std::runtime_error("Paper-matched serial sweeps require one MPI rank; use ./exe or mpirun -np 1");
    std::ofstream out,pts;
    if(rank==0) {
      out.open(o.csv.c_str(),std::ios::trunc);
      if(!out)throw std::runtime_error("Cannot open "+o.csv);
      out<<"study,paper_figures,geometry_id,geometry_name,boundary_degree,"
           "monomial_q,monomial_total_degree,elements_per_direction,mesh_elements,"
           "nqp_set,kq_requested,shift_step,shift_distance,method,geometry_model,"
           "status,detail,numerical,reference,absolute_error,relative_error,"
           "all_elements,cut_elements,uncut_active_elements,"
           "nominal_gauss_nodes_per_element,quadrature_nodes_cut,"
           "quadrature_nodes_uncut,quadrature_nodes_total,surface_nodes_cut,spatial_dimension\n";
      std::cout<<"Review 2D/3D paper-matched benchmark -> "<<o.csv<<"\n"
        <<"Original FEMuS moment construction retained (NO library changes).\n"
        <<"FEMuS Q2 geometry approximates curves of degree >2.\n"
        <<"Paper nQP sweeps exceeding FEMuS FE gaussOrder 37 are marked SKIPPED.\n";
      if(o.study=="layout" || o.study=="all") {
        pts.open(o.pointsCsv.c_str(),std::ios::trunc);
        if(!pts)throw std::runtime_error("Cannot open "+o.pointsCsv);
        pts<<"method,geometry_id,elements_per_direction,nqp_set,element_id,"
               "x,y,physical_volume_weight,cut_element\n";
      }
    }
    const bool everything=o.study=="all";
    const bool core=everything || o.study=="core";
    if(core || o.study=="area")Area(out,o,rank);
    if(core || o.study=="area-href")AreaHRef(out,o,rank);
    if(core || o.study=="monomial")Monomial(out,o,rank);
    if(core || o.study=="monomial-href")MonomialHRef(out,o,rank);
    if(everything || o.study=="layout")Layout(out,pts,o,rank);
    if(everything || o.study=="suite26")Suite26(out,o,rank);
    if(everything || o.study=="moving")MovingTest(out,o,rank);
    const bool all3d=everything || o.study=="3d";
    if(all3d || o.study=="ellipsoid-moving")EllipsoidMoving3D(out,o,rank);
    if(all3d || o.study=="ellipsoid-mesh-moving")EllipsoidMeshMoving3D(out,o,rank);
    if(all3d || o.study=="torus-href")TorusHRef3D(out,o,rank);
    if(all3d || o.study=="cylinder-monomial")CylinderMonomial3D(out,o,rank);
    if(all3d || o.study=="cylinder-href")CylinderHRef3D(out,o,rank);
  }catch(const std::exception& e) {
    if(rank==0)std::cerr<<"2D/3D paper benchmark error: "<<e.what()<<'\n';
    MPI_Abort(MPI_COMM_WORLD,1);
    return 1;
  }
  return 0;
}
