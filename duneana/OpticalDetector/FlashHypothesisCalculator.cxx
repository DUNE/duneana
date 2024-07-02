#include "FlashHypothesisCalculator.h"

#include "larana/OpticalDetector/FlashHypothesis.h"
#include "lardataalg/Utilities/MappedContainer.h"

#include "TVector3.h"

std::vector<double> opdet::FlashHypothesisCalculator::SegmentMidpoint(TVector3 const& pt1,
                                                                      TVector3 const& pt2,
                                                                      float XOffset)
{
  std::vector<double> xyz_segment(3);
  xyz_segment[0] = 0.5 * (pt2.x() + pt1.x()) + XOffset;
  xyz_segment[1] = 0.5 * (pt2.y() + pt1.y());
  xyz_segment[2] = 0.5 * (pt2.z() + pt1.z());
  return xyz_segment;
}

void opdet::FlashHypothesisCalculator::FillFlashHypothesis(const float& yield,
                                                           const float& dEdx,
                                                           const TVector3& pt1,
                                                           const TVector3& pt2,
                                                           const std::vector<float>& qe_vector,
                                                           phot::MappedCounts_t const& vis_vector,
                                                           FlashHypothesis& hyp)
{

  if (qe_vector.size() != hyp.GetVectorSize() || !vis_vector)
    throw std::runtime_error("ERROR in FlashHypothesisCalculator: vector sizes not equal!");

  const float total_yield = yield * dEdx * (pt2 - pt1).Mag();

  for (size_t i_chan = 0; i_chan < hyp.GetVectorSize(); i_chan++)
    hyp.SetHypothesisAndError(i_chan, total_yield * vis_vector[i_chan] * qe_vector[i_chan]);
}
