/**
 * \file PMTAna_module.cc
 *
 * \ingroup PMTAna
 *
 * \brief Class definition file of PMTAna
 *
 * @author Kazu - Nevis 2013
 */

/** \addtogroup PMTAna

@{*/

// ART includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/fwd.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft
#include "lardataobj/RawData/OpDetWaveform.h"

// STL
#include <functional>
#include <numeric>
#include <string>

// ROOT
#include <TTree.h>

// My modules
#include "OpHitFinder/AlgoFixedWindow.h"
#include "OpHitFinder/AlgoThreshold.h"
#include "OpHitFinder/PedAlgoEdges.h"
#include "OpHitFinder/PulseRecoManager.h"

namespace pmtana {

  /**
     \class PMTAna
     PMTAna module to copy LArSoft data contents into LArLight data formatted file
  */
  class PMTAna : public art::EDAnalyzer {

  public:
    /// Constructor
    PMTAna(const fhicl::ParameterSet&);

    /// Function to be called per event
    void analyze(const art::Event&);

  private:
    std::string _fifo_mod_name; ///< Input FIFOChannel producer name
    TTree* _tree;               ///< output data holder TTree

    PulseRecoManager _preco_man;
    AlgoThreshold _th_algo;
    AlgoFixedWindow _fw_algo;
    PedAlgoEdges _ped_algo;
  };

}

namespace pmtana {
  DEFINE_ART_MODULE(PMTAna)
}

namespace pmtana {

  //#######################################################################################################
  PMTAna::PMTAna(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset), _preco_man(), _th_algo(), _fw_algo(), _ped_algo()
  //#######################################################################################################
  {

    // Obtain module names for input data
    _fifo_mod_name = pset.get<std::string>("fModName_FIFOChannel");

    // Next we make storage data class objects for those data types specified in fcl files.
    art::ServiceHandle<art::TFileService const> fileService;

    // Create TTree
    _tree = fileService->make<TTree>("pmt_tree", "Analysis Tree");

    //
    // Demonstration purpose ...
    //
    _preco_man.AddRecoAlgo(&_th_algo);
    _preco_man.AddRecoAlgo(&_fw_algo);
    _preco_man.SetDefaultPedAlgo(&_ped_algo);
  }

  //#######################################################################################################
  void PMTAna::analyze(const art::Event& evt)
  //#######################################################################################################
  {

    //data_ptr->set_event(evt.id().event(), evt.run(), evt.subRun());

    //    std::vector<const optdata::FIFOChannel*> pmtArray;
    std::vector<const raw::OpDetWaveform*> pmtArray;
    try {

      evt.getView(_fifo_mod_name, pmtArray);
    }
    catch (art::Exception const& e) {

      if (e.categoryCode() != art::errors::ProductNotFound) throw;
    }

    for (size_t i = 0; i < pmtArray.size(); ++i) {

      //      const optdata::FIFOChannel* fifo_ptr(pmtArray.at(i));
      const raw::OpDetWaveform* fifo_ptr(pmtArray.at(i));

      _preco_man.Reconstruct(*fifo_ptr);

      //
      // here I add code to store reco-ed pulse w/ channel number
      // OR I may make a singleton storage manager...

      /*
      data_ptr->add_pmtfifo(fifo_ptr->ChannelNumber(),
			    fifo_ptr->Category(),
			    fifo_ptr->Frame(),
			    fifo_ptr->TimeSlice(),
			    *fifo_ptr);
			    */

      //
      //
      //
    }
  }

}

/** @}*/ // end of PMTAna group
