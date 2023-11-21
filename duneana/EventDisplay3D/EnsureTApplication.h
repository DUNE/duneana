#ifndef EventDisplay3D_EnsureTApplication_h
#define EventDisplay3D_EnsureTApplication_h
//
// If an art module wishes to open a TCanvas on the screen
// then the interactive ROOT environment must be created.  However
// it must only be created once per art job.
//
// The purpose of this class is to check to see if that environment
// has been created.  If it has, do nothing.  If it has not, create it.
//

namespace dune {

  class EnsureTApplication {

  public:
    explicit EnsureTApplication(int argc = 0, char** argv = nullptr);

  private:
  };

} // end namespace dune

#endif /* Analyzers_EnsureTApplication_h */
