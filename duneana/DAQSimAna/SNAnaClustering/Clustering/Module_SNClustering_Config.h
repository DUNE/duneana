#ifndef SNCLUSTERING_CONFIG_H
#define SNCLUSTERING_CONFIG_H

std::vector<int>   cut_AdjChanTolerance = {1,2,2,2,2,2};
std::vector<int>   cut_HitsInWindow     = {2,3,3,4,5,6};
std::vector<int>   cut_MinChannels      = {2,2,2,2,2,2};
std::vector<int>   cut_MinChanWidth     = {0,0,0,0,0,0};
std::vector<float> cut_TimeWindowSize   = {20,20,20,20,20,20};
std::vector<float> cut_TotalADC         = {350,400,450,400,400,0};

unsigned int NConfigs = cut_AdjChanTolerance.size();
unsigned int NCuts    = 6;

double detectorScaling = 0.12;

#endif
