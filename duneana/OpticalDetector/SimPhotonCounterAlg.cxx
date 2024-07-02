/*!
 * Title:   SimPhotonCounterALG Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: Alg class that counts up sim photons, leading towards
 *              comparisons with flashes and flash hypotheses.
*/

#include "SimPhotonCounterAlg.h"
#include "OpDetResponseInterface.h"
#include "OpDigiProperties.h"

#include "fhiclcpp/ParameterSet.h"

opdet::SimPhotonCounterAlg::SimPhotonCounterAlg(fhicl::ParameterSet const& p)
{
  FillAllRanges(p.get<std::vector<fhicl::ParameterSet>>("SimPhotonCounterParams"));
}

void opdet::SimPhotonCounterAlg::FillAllRanges(std::vector<fhicl::ParameterSet> const& pv)
{
  fTimeRanges.clear();
  fWavelengthRanges.clear();

  fTimeRanges.reserve(pv.size());
  fWavelengthRanges.reserve(pv.size());

  for (auto const& p : pv)
    FillRanges(p);
}

void opdet::SimPhotonCounterAlg::FillRanges(fhicl::ParameterSet const& p)
{
  std::vector<float> time_range(4);
  time_range[0] = p.get<float>("MinPromptTime");
  time_range[1] = p.get<float>("MaxPromptTime");
  time_range[2] = p.get<float>("MinLateTime");
  time_range[3] = p.get<float>("MaxLateTime");

  if (time_range[0] > time_range[1] || time_range[2] > time_range[3] ||
      time_range[1] > time_range[2])
    throw std::runtime_error("ERROR in SimPhotonCounterAlg: Bad time range.");

  fTimeRanges.push_back(time_range);

  std::vector<float> wavelength_range(2);
  wavelength_range[0] = p.get<float>("MinWavelength");
  wavelength_range[1] = p.get<float>("MaxWavelength");

  if (wavelength_range[0] >= wavelength_range[1])
    throw std::runtime_error("ERROR in SimPhotonCounterAlg: Bad wavelength range.");

  fWavelengthRanges.push_back(wavelength_range);
}

void opdet::SimPhotonCounterAlg::InitializeCounters(geo::GeometryCore const& geo,
                                                    opdet::OpDigiProperties const& opdigip)
{
  fCounters.resize(fTimeRanges.size());
  art::ServiceHandle<opdet::OpDetResponseInterface const> odresponse;
  for (size_t i = 0; i < fCounters.size(); i++)
    fCounters[i] = SimPhotonCounter(fTimeRanges[i][0],
                                    fTimeRanges[i][1],
                                    fTimeRanges[i][2],
                                    fTimeRanges[i][3],
                                    fWavelengthRanges[i][0],
                                    fWavelengthRanges[i][1],
                                    std::vector<float>(odresponse->NOpChannels(), opdigip.QE()));
}

void opdet::SimPhotonCounterAlg::AddSimPhotonCollection(sim::SimPhotonsCollection const& ph_col)
{
  if (ph_col.size() != fCounters.at(0).GetVectorSize())
    throw std::runtime_error(
      "ERROR in SimPhotonCounterAlg: Photon collection size and OpDet size not equal.");

  for (auto const& photons : ph_col)
    for (auto& counter : fCounters)
      counter.AddSimPhotons(photons.second);
}

void opdet::SimPhotonCounterAlg::AddSimPhotonsVector(std::vector<sim::SimPhotons> const& spv)
{
  for (auto const& photons : spv)
    for (auto& counter : fCounters)
      counter.AddSimPhotons(photons);
}

void opdet::SimPhotonCounterAlg::ClearCounters()
{
  for (auto& counter : fCounters)
    counter.ClearVectors();
}

std::vector<float> const& opdet::SimPhotonCounterAlg::PromptPhotonVector(size_t i)
{
  return fCounters.at(i).PromptPhotonVector();
}

std::vector<float> const& opdet::SimPhotonCounterAlg::LatePhotonVector(size_t i)
{
  return fCounters.at(i).LatePhotonVector();
}

opdet::SimPhotonCounter const& opdet::SimPhotonCounterAlg::GetSimPhotonCounter(size_t i)
{
  return fCounters.at(i);
}
