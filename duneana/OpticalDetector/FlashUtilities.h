#ifndef FLASHUTILITIES_H
#define FLASHUTILITIES_H

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

#include <vector>

namespace opdet {

  class FlashHypothesis;

  class FlashUtilities {

  public:
    FlashUtilities() {}

    float CompareByError(const FlashHypothesis&, const std::vector<float>&, std::vector<float>&);
    float CompareByFraction(const FlashHypothesis&, const std::vector<float>&, std::vector<float>&);
    float CompareByFraction(const std::vector<float>&,
                            const std::vector<float>&,
                            std::vector<float>&);
    void GetPosition(const std::vector<float>&, const std::vector<float>&, float&, float&);
    void GetPosition(const std::vector<float>&, const std::vector<float>&, double&, double&);

  private:
  };

}

#endif
