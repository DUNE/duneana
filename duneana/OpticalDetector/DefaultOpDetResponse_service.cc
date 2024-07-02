////////////////////////////////////////////////////////////////////////
//
//  \file DefaultOpDetResponse_service.cc
//
////////////////////////////////////////////////////////////////////////

#include "larana/OpticalDetector/OpDetResponseInterface.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <cassert>

namespace opdet {
  class DefaultOpDetResponse : public OpDetResponseInterface {
  public:
    DefaultOpDetResponse(fhicl::ParameterSet const& pset);

  private:
    void doReconfigure(fhicl::ParameterSet const& p) override;
    bool doDetected(int OpChannel, const sim::OnePhoton& Phot, int& newOpChannel) const override;
    bool doDetectedLite(int OpChannel, int& newOpChannel) const override;

  }; // class DefaultOpDetResponse

}

DECLARE_ART_SERVICE_INTERFACE_IMPL(opdet::DefaultOpDetResponse,
                                   opdet::OpDetResponseInterface,
                                   LEGACY)

namespace opdet {
  //--------------------------------------------------------------------
  DefaultOpDetResponse::DefaultOpDetResponse(fhicl::ParameterSet const& pset)
  {
    this->doReconfigure(pset);
  }

  //--------------------------------------------------------------------
  void DefaultOpDetResponse::doReconfigure(fhicl::ParameterSet const& pset)
  {
    auto const* LarProp = lar::providerFrom<detinfo::LArPropertiesService>();

    if (LarProp->ScintPreScale() < 1) {
      mf::LogWarning("DefaultOpDetResponse_service")
        << "A prescale of " << LarProp->ScintPreScale()
        << " has been applied during optical MC production, "
        << "but DefaultOpDetResponse does not include any QE so this effect is not being corrected "
           "out.";
      assert(false);
    }
  }

  //--------------------------------------------------------------------
  bool DefaultOpDetResponse::doDetected(int OpChannel,
                                        const sim::OnePhoton& /*Phot*/,
                                        int& newOpChannel) const
  {
    newOpChannel = OpChannel;
    return true;
  }

  //--------------------------------------------------------------------
  bool DefaultOpDetResponse::doDetectedLite(int OpChannel, int& newOpChannel) const
  {
    newOpChannel = OpChannel;
    return true;
  }

} // namespace

DEFINE_ART_SERVICE_INTERFACE_IMPL(opdet::DefaultOpDetResponse, opdet::OpDetResponseInterface)
