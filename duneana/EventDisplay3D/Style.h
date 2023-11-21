////////////////////////////////////////////////////////////////////////
/// \file Style.h
//
/// \author messier@indiana.edu
////////////////////////////////////////////////////////////////////////
#ifndef EventDisplay3D_LINESTYLE_H
#define EventDisplay3D_LINESTYLE_H
class TLine;

namespace dune {
  /// Parameters for drawing options. Allow a consistent style for
  /// drawing particle tracks
  class Style {
  public:
    static const char* LatexName(int pdgcode);
    static void        FromPDG(TLine& line, int pdgcode);
    static int         ColorFromPDG(int pdgcode);
    static int         LineStyleFromPDG(int pdgcode);
    static int         LineWidthFromPDG(int pdgcode);
  };
}
#endif
