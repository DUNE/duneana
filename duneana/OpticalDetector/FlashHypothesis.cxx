#include "FlashHypothesis.h"

#include <iostream>
#include <limits>

void opdet::FlashHypothesis::SetHypothesisVectorAndErrorVector(std::vector<float> v,
                                                               std::vector<float> err)
{
  if (err.size() != 0 && err.size() != v.size())
    throw std::runtime_error("ERROR in FlashHypothesisVectorSetter: Vector sizes not equal");

  _NPEs_Vector = v;
  if (err.size() == 0) {
    _NPEs_Vector = v;
    _NPEs_ErrorVector.resize(v.size());
    for (size_t i = 0; i < _NPEs_Vector.size(); i++)
      _NPEs_ErrorVector[i] = std::sqrt(_NPEs_Vector[i]);
  }
  else
    _NPEs_ErrorVector = err;
}

void opdet::FlashHypothesis::SetHypothesisAndError(size_t i_opdet, float pe, float err)
{
  SetHypothesis(i_opdet, pe);
  if (err > 0)
    SetHypothesisError(i_opdet, err);
  else
    SetHypothesisError(i_opdet, std::sqrt(pe));
}

void opdet::FlashHypothesis::Normalize(float const& totalPE_target)
{

  if (GetTotalPEs() < std::numeric_limits<float>::epsilon()) return;

  const float PE_ratio = totalPE_target / GetTotalPEs();
  for (size_t i_opdet = 0; i_opdet < _NPEs_Vector.size(); i_opdet++) {
    _NPEs_Vector[i_opdet] *= PE_ratio;
    _NPEs_ErrorVector[i_opdet] = std::sqrt(_NPEs_Vector[i_opdet]);
  }
}

void opdet::FlashHypothesis::Print()
{
  std::cout << "TotalPEs: " << GetTotalPEs() << " +/- " << GetTotalPEsError() << std::endl;
  std::cout << "Vector size: " << GetVectorSize() << std::endl;
  for (size_t i = 0; i < GetVectorSize(); i++)
    std::cout << "\t" << i << ": " << GetHypothesis(i) << " +/- " << GetHypothesisError(i)
              << std::endl;
}

void opdet::FlashHypothesisCollection::SetTotalHypAndPromptFraction(const FlashHypothesis& total,
                                                                    float frac)
{
  CheckFrac(frac);
  _total_hyp = total;
  _prompt_frac = frac;
  _prompt_hyp = total;
  _prompt_hyp.Normalize(frac * total.GetTotalPEs());
  _late_hyp = total;
  _late_hyp.Normalize((1. - frac) * total.GetTotalPEs());
}

void opdet::FlashHypothesisCollection::SetPromptHypAndPromptFraction(const FlashHypothesis& prompt,
                                                                     float frac)
{
  CheckFrac(frac);
  _prompt_hyp = prompt;
  _prompt_frac = frac;
  _late_hyp = prompt;
  _late_hyp.Normalize((1 / frac - 1.) * prompt.GetTotalPEs());
  _total_hyp = _prompt_hyp + _late_hyp;
}

void opdet::FlashHypothesisCollection::Normalize(float totalPE_target)
{
  _prompt_hyp.Normalize(totalPE_target * _prompt_frac);
  _late_hyp.Normalize(totalPE_target * (1. - _prompt_frac));
  UpdateTotalHyp();
}

void opdet::FlashHypothesisCollection::CheckFrac(float f)
{
  if (std::abs(f - 0.0) < std::numeric_limits<float>::epsilon() ||
      std::abs(f - 1.0) < std::numeric_limits<float>::epsilon() || (f > 0.0 && f < 1.0))
    return;

  throw std::runtime_error(
    "ERROR in FlashHypothesisCollection: Input fraction is not in valid range.");
}

void opdet::FlashHypothesisCollection::UpdateTotalHyp()
{
  _total_hyp = _prompt_hyp + _late_hyp;
  const float total_pe = _total_hyp.GetTotalPEs();
  if (total_pe > std::numeric_limits<float>::epsilon())
    _prompt_frac = _prompt_hyp.GetTotalPEs() / total_pe;
  else
    _prompt_frac = 1.;
}

void opdet::FlashHypothesisCollection::Print()
{
  std::cout << "PromptHyp:" << std::endl;
  _prompt_hyp.Print();
  std::cout << "LateHyp:" << std::endl;
  _late_hyp.Print();
  std::cout << "TotalHyp:" << std::endl;
  _total_hyp.Print();
  std::cout << "PromptFraction: " << _prompt_frac << std::endl;
}
