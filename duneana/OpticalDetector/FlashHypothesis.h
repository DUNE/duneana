#ifndef FLASHHYPOTHESIS_H
#define FLASHHYPOTHESIS_H

/*!
 * Title:   FlashHypothesis Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: Class that represents a flash hypothesis (PEs per opdet) and its error
*/

#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace opdet {

  class FlashHypothesis {

  public:
    FlashHypothesis() {}
    FlashHypothesis(size_t s)
    {
      _NPEs_Vector = std::vector<float>(s, 0.0);
      _NPEs_ErrorVector = std::vector<float>(s, 0.0);
    }
    FlashHypothesis(std::vector<float> const& vector,
                    std::vector<float> const& vector_error = std::vector<float>())
    {
      SetHypothesisVectorAndErrorVector(vector, vector_error);
    }

    std::vector<float> const& GetHypothesisVector() const { return _NPEs_Vector; }
    std::vector<float> const& GetHypothesisErrorVector() const { return _NPEs_ErrorVector; }
    void SetHypothesisVector(std::vector<float> v)
    {
      _NPEs_Vector = v;
      _NPEs_ErrorVector.resize(v.size());
    }
    void SetHypothesisErrorVector(std::vector<float> v)
    {
      _NPEs_ErrorVector = v;
      _NPEs_Vector.resize(v.size());
    }
    void SetHypothesisVectorAndErrorVector(std::vector<float> v,
                                           std::vector<float> err = std::vector<float>(0));

    float const& GetHypothesis(size_t i_opdet) const { return _NPEs_Vector.at(i_opdet); }
    float const& GetHypothesisError(size_t i_opdet) const { return _NPEs_ErrorVector.at(i_opdet); }
    void SetHypothesis(size_t i_opdet, float pe) { _NPEs_Vector.at(i_opdet) = pe; }
    void SetHypothesisError(size_t i_opdet, float err) { _NPEs_ErrorVector.at(i_opdet) = err; }

    void SetHypothesisAndError(size_t i_opdet, float pe, float err = -999);

    float GetTotalPEs() const
    {
      return std::accumulate(_NPEs_Vector.begin(), _NPEs_Vector.end(), 0.0);
    }
    float GetTotalPEsError() const
    {
      return std::sqrt(std::inner_product(
        _NPEs_ErrorVector.begin(), _NPEs_ErrorVector.end(), _NPEs_ErrorVector.begin(), 0.0));
    }

    size_t GetVectorSize() const { return _NPEs_Vector.size(); }

    void Normalize(float const& totalPE_target);

    void Print();

    FlashHypothesis operator+(const FlashHypothesis& fh)
    {

      if (_NPEs_Vector.size() != fh.GetVectorSize())
        throw std::runtime_error(
          "ERROR in FlashHypothesisAddition: Cannot add hypothesis of different size");

      FlashHypothesis flashhyp(_NPEs_Vector.size());
      for (size_t i = 0; i < _NPEs_Vector.size(); i++) {
        flashhyp._NPEs_Vector[i] = _NPEs_Vector[i] + fh._NPEs_Vector[i];
        flashhyp._NPEs_ErrorVector[i] =
          std::sqrt(this->_NPEs_ErrorVector[i] * this->_NPEs_ErrorVector[i] +
                    fh._NPEs_ErrorVector[i] * fh._NPEs_ErrorVector[i]);
      }
      return flashhyp;
    }

  private:
    std::vector<float> _NPEs_Vector;
    std::vector<float> _NPEs_ErrorVector;
  };

  class FlashHypothesisCollection {

  public:
    FlashHypothesisCollection() {}
    FlashHypothesisCollection(size_t s)
    {
      _prompt_hyp = FlashHypothesis(s);
      _late_hyp = FlashHypothesis(s);
      UpdateTotalHyp();
    }
    FlashHypothesisCollection(const FlashHypothesis& prompt, const FlashHypothesis& late)
    {
      SetPromptAndLateHyp(prompt, late);
    }

    void SetPromptAndLateHyp(const FlashHypothesis& prompt, const FlashHypothesis& late)
    {
      _prompt_hyp = prompt;
      _late_hyp = late;
      UpdateTotalHyp();
    }
    void SetTotalHypAndPromptFraction(const FlashHypothesis& total, float frac);
    void SetPromptHypAndPromptFraction(const FlashHypothesis& prompt, float frac);

    size_t GetVectorSize() const { return _prompt_hyp.GetVectorSize(); }

    float GetPromptFraction() const { return _prompt_frac; }
    float GetLateFraction() const { return 1. - _prompt_frac; }

    const FlashHypothesis& GetPromptHypothesis() const { return _prompt_hyp; }
    const FlashHypothesis& GetLateHypothesis() const { return _late_hyp; }
    const FlashHypothesis& GetTotalHypothesis() const { return _total_hyp; }

    void Normalize(float totalPEs);

    void Print();

    FlashHypothesisCollection operator+(const FlashHypothesisCollection& fhc)
    {

      if (this->GetVectorSize() != fhc.GetVectorSize())
        throw std::runtime_error(
          "ERROR in FlashHypothesisCollectionAddition: Cannot add hypothesis of different size");

      FlashHypothesis ph = this->GetPromptHypothesis();
      ph = ph + fhc.GetPromptHypothesis();
      FlashHypothesis lh = this->GetLateHypothesis();
      lh = lh + fhc.GetLateHypothesis();

      return (FlashHypothesisCollection(ph, lh));
    }

  private:
    FlashHypothesis _prompt_hyp;
    FlashHypothesis _late_hyp;
    FlashHypothesis _total_hyp;
    float _prompt_frac;

    void CheckFrac(float f);
    void UpdateTotalHyp();
  };

}

#endif
