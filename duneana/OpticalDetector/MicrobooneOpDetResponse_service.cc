////////////////////////////////////////////////////////////////////////
//
//  \file MicrobooneOpDetResponse_service.cc
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

#include "larana/OpticalDetector/OpDetResponseInterface.h"

#include "lardataobj/Simulation/SimPhotons.h"

namespace opdet {
  class MicrobooneOpDetResponse : public OpDetResponseInterface {
  public:
    MicrobooneOpDetResponse(fhicl::ParameterSet const& pset);

  private:
    void doReconfigure(fhicl::ParameterSet const& p) override;
    bool doDetected(int OpChannel, const sim::OnePhoton& Phot, int& newOpChannel) const override;
    bool doDetectedLite(int OpChannel, int& newOpChannel) const override;

    float fQE; // Quantum efficiency of tube

    float fWavelengthCutLow;  // Sensitive wavelength range
    float fWavelengthCutHigh; //
  };                          // class MicrobooneOpDetResponse
}

DECLARE_ART_SERVICE_INTERFACE_IMPL(opdet::MicrobooneOpDetResponse,
                                   opdet::OpDetResponseInterface,
                                   LEGACY)

namespace opdet {
  //--------------------------------------------------------------------
  MicrobooneOpDetResponse::MicrobooneOpDetResponse(fhicl::ParameterSet const& pset)
  {
    this->doReconfigure(pset);
  }

  //--------------------------------------------------------------------
  void MicrobooneOpDetResponse::doReconfigure(fhicl::ParameterSet const& pset)
  {
    fQE = pset.get<double>("QuantumEfficiency");
    //double tempfQE=            pset.get<double>("QuantumEfficiency");
    fWavelengthCutLow = pset.get<double>("WavelengthCutLow");
    fWavelengthCutHigh = pset.get<double>("WavelengthCutHigh");

    /**
     * Don't apply QE here.  It is applied in the uboone
     * electronics simulation.
     **
     // Correct out the prescaling applied during simulation
     auto const* LarProp = lar::providerFrom<detinfo::LArPropertiesService>();
     fQE = tempfQE / LarProp->ScintPreScale();

     if (fQE > 1.0001 ) {
         mf::LogWarning("MicrobooneOpDetResponse_service") << "Quantum efficiency set in OpDetResponse_service, " << tempfQE
                                                           << " is too large.  It is larger than the prescaling applied during simulation, "
                                                           << LarProp->ScintPreScale()
                                                           << ".  Final QE must be equalt to or smaller than the QE applied at simulation time.";
         assert(false);
     }
     **/
  }

  //--------------------------------------------------------------------
  bool MicrobooneOpDetResponse::doDetected(int OpChannel,
                                           const sim::OnePhoton& Phot,
                                           int& newOpChannel) const
  {

    newOpChannel = OpChannel;

    /**
     * Don't apply QE here.  It is applied in the uboone
     * electronics simulation.
     **
     // Check QE
     if ( CLHEP::RandFlat::shoot(1.0) > fQE ) return false;
     **/

    double wavel = wavelength(Phot.Energy);
    // Check wavelength acceptance
    if (wavel < fWavelengthCutLow) return false;
    if (wavel > fWavelengthCutHigh) return false;

    return true;
  }

  //--------------------------------------------------------------------
  bool MicrobooneOpDetResponse::doDetectedLite(int OpChannel, int& newOpChannel) const
  {
    newOpChannel = OpChannel;

    /**
     * Don't apply QE here.  It is applied in the uboone
     * electronics simulation.
     **
     // Check QE
     if ( CLHEP::RandFlat::shoot(1.0) > fQE ) return false;
     **/

    return true;
  }

} // namespace

DEFINE_ART_SERVICE_INTERFACE_IMPL(opdet::MicrobooneOpDetResponse, opdet::OpDetResponseInterface)
