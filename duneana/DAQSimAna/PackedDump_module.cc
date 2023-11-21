////////////////////////////////////////////////////////////////////////
// Class:       PackedDump
// Plugin Type: analyzer (art v2_10_03)
// File:        PackedDump_module.cc
//
// Generated at Thu Jun 14 08:21:34 2018 by Philip Rodrigues using cetskelgen
// from cetlib version v3_02_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RawData/RawDigit.h"
#include "dunepdlegacy/Services/ChannelMap/PdspChannelMapService.h"

#include <fstream>
#include <arpa/inet.h>

class PackedDump;


class PackedDump : public art::EDAnalyzer {
public:
    explicit PackedDump(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    PackedDump(PackedDump const &) = delete;
    PackedDump(PackedDump &&) = delete;
    PackedDump & operator = (PackedDump const &) = delete;
    PackedDump & operator = (PackedDump &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;

private:

    art::ServiceHandle<dune::PdspChannelMapService> m_channelMap;
    std::string m_outputFilename;
    std::ofstream m_outputFile;
};


PackedDump::PackedDump(fhicl::ParameterSet const & p)
    :
    EDAnalyzer(p),
    m_outputFilename(p.get<std::string>("OutputFile")),
    m_outputFile(m_outputFilename, std::ios::out | std::ios::binary)
{
    // A little test code to check we can do the bit packing properly

    // const size_t nsamples=16;
    // const uint16_t samples[nsamples]={
    //     0, 1,  2,  3,  4,  5,  6,  7,
    //     8, 9, 10, 11, 12, 13, 14, 15
    // };
    // for(size_t i=0; i<nsamples; i+=8){
    //     uint32_t word0=samples[i] | (samples[i+1] << 12) | (samples[i+2] << 24);
    //     uint32_t word1=(samples[i+2] >> 4) | (samples[i+3] << 4) | (samples[i+4] << 16) | (samples[i+5] << 28);
    //     uint32_t word2=(samples[i+5] >> 4) | (samples[i+6] << 8) | (samples[i+7] << 20);

    //     m_outputFile.write((char*)&word0, sizeof(uint32_t));
    //     m_outputFile.write((char*)&word1, sizeof(uint32_t));
    //     m_outputFile.write((char*)&word2, sizeof(uint32_t));

    //     uint32_t word0_net=htonl(word0);
    //     uint32_t word1_net=htonl(word1);
    //     uint32_t word2_net=htonl(word2);

    //     m_outputFile.write((char*)&word0_net, sizeof(uint32_t));
    //     m_outputFile.write((char*)&word1_net, sizeof(uint32_t));
    //     m_outputFile.write((char*)&word2_net, sizeof(uint32_t));
    // }
}

void PackedDump::analyze(art::Event const & e)
{
    // Create a map from offline channel number to RawDigit
    std::vector<const raw::RawDigit*> inputDigitsHandle;
    std::map<unsigned int, const raw::RawDigit*> channelToDigit;
    e.getView("daq", inputDigitsHandle);
    for(size_t c=0; c<inputDigitsHandle.size(); ++c){
        const raw::RawDigit* dig=inputDigitsHandle[c];
        channelToDigit[dig->Channel()] = dig;
    }
    
    // Set a bunch of these to 1 to just do a subset to start with
    const size_t nCrate=1;
    const size_t nSlot=1;
    const size_t nFiber=1;
    const size_t nFEMBChannel=128;
    const size_t nTDC=4492;

    // We loop down to the level of a board and then, for a given
    // board, emit frame. Each frame starts with 0xdeadbeef, followed
    // by the tdc, followed by each 12-bit FMB channel value, in order
    // of FEMB channel
    for(unsigned int crate=0; crate<nCrate; ++crate){
        for(unsigned int slot=0; slot<nSlot; ++slot){
            for(unsigned int fiber=0; fiber<nFiber; ++fiber){
                for(uint32_t tdc=0; tdc<nTDC; ++tdc){
                    uint32_t header=0xdeadbeef;
                    m_outputFile.write((char*)&header, sizeof(uint32_t));
                    m_outputFile.write((char*)&tdc, sizeof(uint32_t));
                    for(unsigned int fembChannel=0; fembChannel<nFEMBChannel-1; fembChannel+=8){

                        unsigned short samples[8];
                        for(int i=0; i<8; ++i){
                            // args are: unsigned int crate, unsigned int slot, unsigned int fiber, unsigned int fembchannel
                            unsigned int offlineChan=m_channelMap->GetOfflineNumberFromDetectorElements(crate, slot, fiber, fembChannel+i, dune::PdspChannelMapService::kRCE);
                            // std::cout << "Fiber " << fiber << " FEMB channel " << (fembChannel+8) << " maps to offline " << offlineChan << std::endl;
                            samples[i]=channelToDigit[offlineChan]->ADC(tdc);
                        }
                        uint32_t word0=samples[0] | (samples[1] << 12) | (samples[2] << 24);
                        uint32_t word1=(samples[2] >> 8) | (samples[3] << 4) | (samples[4] << 16) | (samples[5] << 28);
                        uint32_t word2=(samples[5] >> 4) | (samples[6] << 8) | (samples[7] << 20);

                        m_outputFile.write((char*)&word0, sizeof(uint32_t));
                        m_outputFile.write((char*)&word1, sizeof(uint32_t));
                        m_outputFile.write((char*)&word2, sizeof(uint32_t));
                    }
                }
            }
        }
    }
    uint32_t eof=0xffffffff;
    m_outputFile.write((char*)&eof, sizeof(uint32_t));
}

DEFINE_ART_MODULE(PackedDump)

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
