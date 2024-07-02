/*!
 * Title:   FlashUtilities Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description:
 * Class that contains utility functions for flash and flash hypotheses:
 * --- compare a flash hypothesis to a truth or reco vector
 * --- get an extent of a flash (central point, width)
 * These classes should operate using simple objects, and will need other
 * classes/functions to fill those vectors properly.
*/

#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

#include "FlashHypothesis.h"
#include "FlashUtilities.h"

float opdet::FlashUtilities::CompareByError(const FlashHypothesis& fh,
                                            const std::vector<float>& compare_vector,
                                            std::vector<float>& result_vector)
{
  if (compare_vector.size() != fh.GetVectorSize())
    throw std::runtime_error("ERROR in FlashUtilities Compare: Mismatching vector sizes.");

  result_vector.resize(fh.GetVectorSize());

  float result_total = 0;

  std::vector<float> const& NPEs_Vector(fh.GetHypothesisVector());
  std::vector<float> const& NPEs_ErrorVector(fh.GetHypothesisErrorVector());

  for (size_t i = 0; i < fh.GetVectorSize(); i++) {
    result_total += compare_vector[i];
    float diff = NPEs_Vector[i] - compare_vector[i];
    if (std::abs(diff) < std::numeric_limits<float>::epsilon())
      result_vector[i] = 0;
    else if (NPEs_ErrorVector[i] < std::numeric_limits<float>::epsilon())
      result_vector[i] = diff / std::numeric_limits<float>::epsilon();
    else
      result_vector[i] = diff / NPEs_ErrorVector[i];
  }

  float total_error = fh.GetTotalPEsError();
  float total_diff = fh.GetTotalPEs() - result_total;
  if (std::abs(total_diff) < std::numeric_limits<float>::epsilon())
    result_total = 0;
  else if (total_error < std::numeric_limits<float>::epsilon())
    result_total = total_diff / std::numeric_limits<float>::epsilon();
  else
    result_total = total_diff / total_error;

  return result_total;
}

float opdet::FlashUtilities::CompareByFraction(const FlashHypothesis& fh,
                                               const std::vector<float>& compare_vector,
                                               std::vector<float>& result_vector)
{
  return CompareByFraction(fh.GetHypothesisVector(), compare_vector, result_vector);
}

float opdet::FlashUtilities::CompareByFraction(const std::vector<float>& NPEs_Vector,
                                               const std::vector<float>& compare_vector,
                                               std::vector<float>& result_vector)
{
  if (compare_vector.size() != NPEs_Vector.size())
    throw std::runtime_error("ERROR in FlashUtilities Compare: Mismatching vector sizes.");

  result_vector.resize(NPEs_Vector.size());

  float total_comp = 0;
  float total_hyp = 0;
  for (size_t i = 0; i < NPEs_Vector.size(); i++) {
    total_comp += compare_vector[i];
    total_hyp += NPEs_Vector[i];
    float diff = NPEs_Vector[i] - compare_vector[i];
    if (std::abs(diff) < std::numeric_limits<float>::epsilon())
      result_vector[i] = 0;
    else if (compare_vector[i] < std::numeric_limits<float>::epsilon())
      result_vector[i] = diff / std::numeric_limits<float>::epsilon();
    else
      result_vector[i] = diff / compare_vector[i];
  }

  float result_total = 0.0;
  float total_diff = total_hyp - total_comp;
  if (std::abs(total_diff) < std::numeric_limits<float>::epsilon())
    result_total = 0;
  else if (total_comp < std::numeric_limits<float>::epsilon())
    result_total = total_diff / std::numeric_limits<float>::epsilon();
  else
    result_total = total_diff / total_comp;

  return result_total;
}

void opdet::FlashUtilities::GetPosition(const std::vector<float>& pe_vector,
                                        const std::vector<float>& pos_vector,
                                        float& mean,
                                        float& rms)
{
  if (pe_vector.size() != pos_vector.size())
    throw std::runtime_error("ERROR in FlashUtilities GetPosition: Mismatchin vector sizes.");

  float sum = std::accumulate(pe_vector.begin(), pe_vector.end(), 0.0);

  if (sum < std::numeric_limits<float>::epsilon()) {
    mean = 0;
    rms = 0;
    return;
  }

  mean = std::inner_product(pe_vector.begin(), pe_vector.end(), pos_vector.begin(), 0.0) / sum;

  rms = 0;
  for (size_t i = 0; i < pe_vector.size(); i++)
    rms += pe_vector[i] * (pos_vector[i] - mean) * (pos_vector[i] - mean);

  rms = std::sqrt(rms) / sum;
}

void opdet::FlashUtilities::GetPosition(const std::vector<float>& pe_vector,
                                        const std::vector<float>& pos_vector,
                                        double& mean,
                                        double& rms)
{
  float fmean, frms;
  GetPosition(pe_vector, pos_vector, fmean, frms);
  mean = double(fmean);
  rms = double(frms);
}
