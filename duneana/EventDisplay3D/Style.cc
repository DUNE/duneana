////////////////////////////////////////////////////////////////////////
/// \file Style.cxx
//
/// \author messier@indiana.edu
////////////////////////////////////////////////////////////////////////
#include "Style.h"
#include "TLine.h"

namespace dune{

  /// Convert PDG code to a latex string (root-style)
  const char* Style::LatexName(int pdgcode)
  {
    switch (pdgcode) {
    case  22:   return "#gamma";
    case -11:   return "e^{+}";
    case  11:   return "e^{-}";
    case  13:   return "#mu";
    case -15:   return "#bar{#tau}";
    case  15:   return "#tau";
    case -13:   return "#bar{#mu}";
    case  12:   return "#nu_{e}";
    case  14:   return "#nu_{#mu}";
    case  16:   return "#nu_{#tau}";
    case -12:   return "#bar{#nu}_{e}";
    case -14:   return "#bar{#nu}_{#mu}";
    case -16:   return "#bar{#nu}_{#tau}";
    case  111:  return "#pi^{0}";
    case  211:  return "#pi^{+}";
    case -211:  return "#pi^{-}";
    case  321:  return "K^{+}";
    case -321:  return "K^{-}";
    case  130:  return "K^{0}_{L}";
    case  310:  return "K^{0}_{S}";
    case  2112: return "n";
    case  2212: return "p";
    case -2112: return "#bar{n}";
    case -2212: return "#bar{p}";
    case 2224:  return "#Delta^{++}";
    case 1000060120: return "^{12}C";
    case 1000170350: return "^{35}Cl";
    case 1000260560: return "^{56}Fe";
    case 1000220480: return "^{48}Ti";
    case 1000080160: return "^{16}O";
    case 1000070140: return "^{14}N";
    case 1000110230: return "^{23}Na";
    case 1000130270: return "^{27}Al";
    case 1000140280: return "^{28}Si";
    case 1000200400: return "^{40}Ca";
    case 1000561370: return "^{137}Ba";
    case 1000180400: return "^{40}Ar";
    case 1000180390: return "^{39}Ar";
    case 2000000001: return "GENIE_{1}";
    case 2000000002: return "GENIE_{2}";
    default:
      static char buff[256];
      sprintf(buff,"X_{%d}",pdgcode);
      return buff;
    }
    return 0;
  }

  //......................................................................

  int Style::ColorFromPDG(int pdgcode) {
    switch (pdgcode) {
      case 11:
      case -11:
      case 12:
      case -12:
        return kYellow;
      case 13:
      case -13:
      case 14:
      case -14:
        return kCyan;
      case 15:
      case -15:
      case 16:
      case -16:
        return kYellow-1;
      case 22:
        return kGreen;
      case 111:
      case 211:
      case -211:
      case  321:
      case -321:
      case  130:
      case  310:
        return kRed;
      case 2112:
      case 2212:
        return kAzure-2;
      default:
        return kGray;
    }
  }

  //............................................................

  int Style::LineWidthFromPDG(int pdgcode) {
    if (pdgcode == 2112 || pdgcode == 2212) return 3;
    return 2;
  }

  //............................................................

  int Style::LineStyleFromPDG(int pdgcode) {
    switch (pdgcode) {
    case 11:
    case -11:
    case 13:
    case -13:
    case 15:
    case -15:
    case 211:
    case -211:
    case 2212:
      return kSolid;
    case 12:
    case -12:
    case 14:
    case -14:
    case 16:
    case -16:
    case 22:
    case 111:
      return kDashed;
    case 2112:
      return kDotted;
    }
    return 0;
  }

  //............................................................

  void Style::FromPDG(TLine& line, int pdgcode)
  {
    // Many cases handled here for most common particles. Extend list as
    // needed
    int kSolid=1, kDashed=2, kDotted=3 /* kDashDot=4 */;
    int c = kGray;
    int s = kDotted;
    int w = 1;

    switch (pdgcode) {
    case  11:  c=kRed;       s=kSolid;  w=2; break; // e-
    case -11:  c=kRed;       s=kSolid;  w=2; break; // e+
    case  12:  c=kRed;       s=kDotted; w=2; break; // nue
    case -12:  c=kRed;       s=kDotted; w=2; break; // nue-bar
    case  13:  c=kBlue;      s=kSolid;  w=2; break; // mu+
    case -13:  c=kBlue;      s=kSolid;  w=2; break; // mu-
    case  14:  c=kBlue;      s=kDotted; w=2; break; // numu
    case -14:  c=kBlue;      s=kDotted; w=2; break; // numu-bar
    case  15:  c=kGreen;     s=kSolid;  w=2; break; // tau+
    case -15:  c=kGreen;     s=kSolid;  w=2; break; // tau-
    case  16:  c=kGreen;     s=kDotted; w=2; break; // nutau
    case -17:  c=kGreen;     s=kDotted; w=2; break; // nutau-bar
    case  22:  c=kYellow-1;  s=kDotted; w=2; break; // gamma
    case  111: c=kMagenta-3; s=kDashed; w=3; break; // pi0
    case  211: c=kMagenta-3; s=kSolid;  w=3; break; // pi+
    case -211: c=kMagenta-3; s=kSolid;  w=3; break; // pi-
    case 2212: c=kMagenta+3; s=kSolid;  w=4; break; // proton
    case 2112: c=kMagenta+3; s=kDotted; w=4; break; // neutron
    default: break;
    };
    line.SetLineColor(c);
    line.SetLineStyle(s);
    line.SetLineWidth(w);
  }
}//namespace
////////////////////////////////////////////////////////////////////////
