////////////////////////////////////////////////////////////////////////
//
//  PMTPedestalBase source
//
////////////////////////////////////////////////////////////////////////

#include "PMTPedestalBase.h"
#include "OpticalRecoException.h"
#include <sstream>
namespace pmtana {

  //**************************************************************
  PMTPedestalBase::PMTPedestalBase(std::string name) : _name(name), _mean_v(), _sigma_v()
  //**************************************************************
  {}

  //*********************************
  PMTPedestalBase::~PMTPedestalBase()
  //*********************************
  {}

  //**********************************************
  const std::string& PMTPedestalBase::Name() const
  //**********************************************
  {
    return _name;
  }

  //************************************************************
  bool PMTPedestalBase::Evaluate(const ::pmtana::Waveform_t& wf)
  //************************************************************
  {
    _mean_v.resize(wf.size(), 0);
    _sigma_v.resize(wf.size(), 0);

    for (size_t i = 0; i < wf.size(); ++i)
      _mean_v[i] = _sigma_v[i] = 0;

    const bool res = ComputePedestal(wf, _mean_v, _sigma_v);

    if (wf.size() != _mean_v.size())
      throw OpticalRecoException("Internal error: computed pedestal mean array length changed!");
    if (wf.size() != _sigma_v.size())
      throw OpticalRecoException("Internal error: computed pedestal sigma array length changed!");

    return res;
  }

  //*******************************************
  double PMTPedestalBase::Mean(size_t i) const
  //*******************************************
  {
    if (i > _mean_v.size()) {
      std::stringstream ss;
      ss << "Invalid index: no pedestal mean exist @ " << i;
      throw OpticalRecoException(ss.str());
    }
    return _mean_v[i];
  }

  //*******************************************
  double PMTPedestalBase::Sigma(size_t i) const
  //*******************************************
  {
    if (i > _sigma_v.size()) {
      std::stringstream ss;
      ss << "Invalid index: no pedestal sigma exist @ " << i;
      throw OpticalRecoException(ss.str());
    }
    return _sigma_v[i];
  }

  //*************************************************
  const PedestalMean_t& PMTPedestalBase::Mean() const
  //*************************************************
  {
    return _mean_v;
  }

  //***************************************************
  const PedestalSigma_t& PMTPedestalBase::Sigma() const
  //***************************************************
  {
    return _sigma_v;
  }
}
