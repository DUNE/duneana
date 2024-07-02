#ifndef SIMPHOTONCOUNTER_H
#define SIMPHOTONCOUNTER_H

/*!
 * Title:   SimPhotonCounter Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: Class that counts up sim photons, leading towards
 *              comparisons with flashes and flash hypotheses.
*/

namespace sim {
  struct OnePhoton;
  class SimPhotons;
}

#include <cstddef>
#include <numeric>
#include <vector>

namespace opdet {

  class SimPhotonCounter {

  public:
    SimPhotonCounter() {}
    SimPhotonCounter(size_t s,
                     float t_p1,
                     float t_p2,
                     float t_l1,
                     float t_l2,
                     float min_w = 0,
                     float max_w = 1e6,
                     float e = 1.0);

    SimPhotonCounter(float t_p1,
                     float t_p2,
                     float t_l1,
                     float t_l2,
                     float min_w,
                     float max_w,
                     const std::vector<float>& eV);

    void SetVectorSize(size_t s)
    {
      _photonVector_prompt.resize(s);
      _photonVector_late.resize(s);
      _qeVector.resize(s);
    }
    size_t GetVectorSize() const { return _photonVector_prompt.size(); }

    void SetWavelengthRanges(float min_w, float max_w);

    float MinWavelength() const { return _min_wavelength; }
    float MaxWavelength() const { return _max_wavelength; }

    void SetTimeRanges(float t_p1, float t_p2, float t_l1, float t_l2);

    float MinPromptTime() const { return _min_prompt_time; }
    float MaxPromptTime() const { return _max_prompt_time; }
    float MinLateTime() const { return _min_late_time; }
    float MaxLateTime() const { return _max_late_time; }

    void SetQE(size_t i, float e) { _qeVector.at(i) = e; }
    float QE(size_t i) const { return _qeVector.at(i); }

    void SetQEVector(const std::vector<float>& eV)
    {
      SetVectorSize(eV.size());
      _qeVector = eV;
    }
    std::vector<float> const& QEVector() const { return _qeVector; }

    void AddOnePhoton(size_t i_opdet, const sim::OnePhoton& photon);
    void AddSimPhotons(const sim::SimPhotons& photons);

    void ClearVectors();
    const std::vector<float>& PromptPhotonVector() const { return _photonVector_prompt; }
    const std::vector<float>& LatePhotonVector() const { return _photonVector_late; }
    float PromptPhotonVector(size_t i) const { return _photonVector_prompt.at(i); }
    float LatePhotonVector(size_t i) const { return _photonVector_late.at(i); }

    std::vector<float> TotalPhotonVector() const;
    float TotalPhotonVector(size_t i) const
    {
      return (PromptPhotonVector(i) + LatePhotonVector(i));
    }

    float PromptPhotonTotal() const
    {
      return std::accumulate(_photonVector_prompt.begin(), _photonVector_prompt.end(), 0.0);
    }
    float LatePhotonTotal() const
    {
      return std::accumulate(_photonVector_late.begin(), _photonVector_late.end(), 0.0);
    }
    float PhotonTotal() const { return (PromptPhotonTotal() + LatePhotonTotal()); }

    void Print();

  private:
    std::vector<float> _photonVector_prompt;
    std::vector<float> _photonVector_late;

    float _min_prompt_time; //in ns
    float _max_prompt_time;
    float _min_late_time;
    float _max_late_time;
    std::vector<float> _qeVector;

    float _min_wavelength; //in nm
    float _max_wavelength;

    float Wavelength(const sim::OnePhoton& ph);
  };

}

#endif
