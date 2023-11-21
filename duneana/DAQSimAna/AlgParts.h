#ifndef ALGPARTS_H
#define ALGPARTS_H

#include <vector>

void do_frugal_update(short& median, int& runningDiff, const short sample, const int ncontig)
{
    if(sample>median) ++runningDiff;
    if(sample<median) --runningDiff;
    
    if(runningDiff > ncontig){
        ++median;
        runningDiff=0;
    }
    if(runningDiff < -1*ncontig){
        --median;
        runningDiff=0;
    }
}

std::vector<short> frugal_pedestal_sigkill(const std::vector<short>& raw_in,
                                           const int lookahead, const int threshold,
                                           const int ncontig)
{
    short median=raw_in[0];
    std::vector<short> ped(raw_in.size(), 0);
    int runningDiff=0;
    // Are we updating the median (because we're not within a hit)?
    bool updating=true;
    for(size_t i=0; i<raw_in.size()-lookahead; ++i){
        short s=raw_in[i]; // The current sample
        short sig_cand=raw_in[i+lookahead]; // The sample a little way ahead
        // Do we later go over the threshold?
        bool cand_above=(sig_cand>median+threshold);
        // Are we currently below the threshold?
        bool current_below=(s<median+threshold);
        // Is there an upcoming transition over threshold? If so, freeze the pedestal...
        if(updating && cand_above){
            updating=false;
        }
        // ...until we fall below the pedestal again
        if( (!updating) && (current_below)){
            updating=true;
        }
        // Do the frugal streaming if we're not in a hit
        if(updating){
            do_frugal_update(median, runningDiff, s, ncontig);
        }
        ped[i]=median;
    }
    // Get the last few samples, which we couldn't do in the main loop
    // because we'd go out-of-bounds
    for(size_t i=raw_in.size()-lookahead; i<raw_in.size(); ++i){
        ped[i]=median;
    }
    return ped;
}

std::vector<short> frugal_iqr(const std::vector<short>& raw_in,
                              const std::vector<short>& median,
                              const int ncontig)
{
    std::vector<short> iqr(raw_in.size(), 0);

    int runningDiffLo=0;
    int runningDiffHi=0;

    short quartileLo=median[0]-1;
    short quartileHi=median[0]+1;

    for(size_t i=0; i<raw_in.size(); ++i){
        short s=raw_in[i]; // The current sample
        short m=median[i];

        if(s<m)
            do_frugal_update(quartileLo, runningDiffLo, s, ncontig);
        if(s>m)
            do_frugal_update(quartileHi, runningDiffHi, s, ncontig);

        iqr[i]=quartileHi-quartileLo;
    }

    return iqr;
}

std::vector<short> frugal_pedestal(const std::vector<short>& raw_in,
                                   const int ncontig)
{
    short median=raw_in[0];
    std::vector<short> ped(raw_in.size(), 0);

    int runningDiff=0;
    for(size_t i=0; i<raw_in.size(); ++i){
        short s=raw_in[i]; // The current sample

        do_frugal_update(median, runningDiff, s, ncontig);

        ped[i]=median;
    }

    return ped;
}

std::vector<short> apply_fir_filter(const std::vector<short>& input,
                                    const size_t ntaps, const short* taps)
{
    std::vector<short> filtered(input.size(), 0);
    for(size_t i=0; i<input.size(); ++i){
        for(size_t j=0; j<ntaps; ++j){
            const size_t index=i>j ? i-j : 0;
            filtered[i]+=input[index]*taps[j];
        }
    }
    return filtered;
}



#endif
