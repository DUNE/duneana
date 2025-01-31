#ifndef G4PROCESSUTIL_H
#define G4PROCESSUTIL_H

enum g4_process_
    {
	kG4primary=0,
	kG4CoupledTransportation=1,
	kG4FastScintillation=2,
	kG4Decay=3,
	kG4anti_neutronInelastic=4,
	kG4neutronInelastic=5,
	kG4anti_protonInelastic=6,
	kG4protonInelastic=7,
	kG4hadInelastic=8,
	kG4pipInelastic=9,
	kG4pimInelastic=10,
	kG4xipInelastic=11,
	kG4ximInelastic=12,
	kG4kaonpInelastic=13,
	kG4kaonmInelastic=14,
	kG4sigmapInelastic=15,
	kG4sigmamInelastic=16,
	kG4kaon0LInelastic=17,
	kG4kaon0SInelastic=18,
	kG4lambdaInelastic=19,
	kG4anti_lambdaInelastic=20,
	kG4He3Inelastic=21,
	kG4ionInelastic=22,
	kG4xi0Inelastic=23,
	kG4alphaInelastic=24,
	kG4tInelastic=25,
	kG4dInelastic=26,
	kG4anti_neutronElastic=27,
	kG4neutronElastic=28,
	kG4anti_protonElastic=29,
	kG4protonElastic=30,
	kG4hadElastic=31,
	kG4pipElastic=32,
	kG4pimElastic=33,
	kG4kaonpElastic=34,
	kG4kaonmElastic=35,
	kG4conv=36,
	kG4phot=37,
	kG4annihil=38,
	kG4nCapture=39,
	kG4nKiller=40,
	kG4muMinusCaptureAtRest=41,
	kG4muIoni=42,
	kG4eBrem=43,
	kG4CoulombScat=44,
	kG4hBertiniCaptureAtRest=45,
	kG4hFritiofCaptureAtRest=46,
	kG4photonNuclear=47,
	kG4muonNuclear=48,
	kG4electronNuclear=49,
	kG4positronNuclear=50,
	kG4compt=51,
	kG4eIoni=52,
	kG4muBrems=53,
	kG4hIoni=54,
	kG4muPairProd=55,
	kG4hPairProd=56,
	kG4LArVoxelReadoutScoringProcess=57,
	kG4ionIoni=58,
	kG4hBrems=59,
	kG4Transportation=60,
	kG4msc=61,
	kG4StepLimiter=62,
	kG4UNKNOWN=63
    };// g4_process_


g4_process_ GetG4ProcessID(const std::string &process_name) {
#define MATCH_PROCESS(name) if (process_name == #name) {return kG4 ## name;}
#define MATCH_PROCESS_NAMED(strname, id) if (process_name == #strname) {return kG4 ## id;}
    MATCH_PROCESS(primary)
	MATCH_PROCESS(CoupledTransportation)
	MATCH_PROCESS(FastScintillation)
	MATCH_PROCESS(Decay)
	MATCH_PROCESS(anti_neutronInelastic)
	MATCH_PROCESS(neutronInelastic)
	MATCH_PROCESS(anti_protonInelastic)
	MATCH_PROCESS(protonInelastic)
	MATCH_PROCESS(hadInelastic)
	MATCH_PROCESS_NAMED(kaon+Inelastic, kaonpInelastic)
	MATCH_PROCESS_NAMED(kaon-Inelastic, kaonmInelastic)
	MATCH_PROCESS_NAMED(kaon+Inelastic, kaonpInelastic)
	MATCH_PROCESS_NAMED(kaon-Inelastic, kaonmInelastic)
	MATCH_PROCESS_NAMED(sigma+Inelastic, sigmapInelastic)
	MATCH_PROCESS_NAMED(sigma-Inelastic, sigmamInelastic)
	MATCH_PROCESS_NAMED(pi+Inelastic, pipInelastic)
	MATCH_PROCESS_NAMED(pi-Inelastic, pimInelastic)
	MATCH_PROCESS_NAMED(xi+Inelastic, xipInelastic)
	MATCH_PROCESS_NAMED(xi-Inelastic, ximInelastic)
	MATCH_PROCESS(kaon0LInelastic)
	MATCH_PROCESS(kaon0SInelastic)
	MATCH_PROCESS(lambdaInelastic)
	MATCH_PROCESS_NAMED(anti-lambdaInelastic, anti_lambdaInelastic)
	MATCH_PROCESS(He3Inelastic)
	MATCH_PROCESS(ionInelastic)
	MATCH_PROCESS(xi0Inelastic)
	MATCH_PROCESS(alphaInelastic)
	MATCH_PROCESS(tInelastic)
	MATCH_PROCESS(dInelastic)
	MATCH_PROCESS(anti_neutronElastic)
	MATCH_PROCESS(neutronElastic)
	MATCH_PROCESS(anti_protonElastic)
	MATCH_PROCESS(protonElastic)
	MATCH_PROCESS(hadElastic)
	MATCH_PROCESS_NAMED(kaon+Elastic, kaonpElastic)
	MATCH_PROCESS_NAMED(kaon-Elastic, kaonmElastic)
	MATCH_PROCESS_NAMED(pi+Elastic, pipElastic)
	MATCH_PROCESS_NAMED(pi-Elastic, pimElastic)
	MATCH_PROCESS(conv)
	MATCH_PROCESS(phot)
	MATCH_PROCESS(annihil)
	MATCH_PROCESS(nCapture)
	MATCH_PROCESS(nKiller)
	MATCH_PROCESS(muMinusCaptureAtRest)
	MATCH_PROCESS(muIoni)
	MATCH_PROCESS(eBrem)
	MATCH_PROCESS(CoulombScat)
	MATCH_PROCESS(hBertiniCaptureAtRest)
	MATCH_PROCESS(hFritiofCaptureAtRest)
	MATCH_PROCESS(photonNuclear)
	MATCH_PROCESS(muonNuclear)
	MATCH_PROCESS(electronNuclear)
	MATCH_PROCESS(positronNuclear)
	MATCH_PROCESS(compt)
	MATCH_PROCESS(eIoni)
	MATCH_PROCESS(muBrems)
	MATCH_PROCESS(hIoni)
	MATCH_PROCESS(ionIoni)
	MATCH_PROCESS(hBrems)
	MATCH_PROCESS(muPairProd)
	MATCH_PROCESS(hPairProd)
	MATCH_PROCESS(LArVoxelReadoutScoringProcess)
	MATCH_PROCESS(Transportation)
	MATCH_PROCESS(msc)
	MATCH_PROCESS(StepLimiter)
	std::cerr << "Error: Process name with no match (" << process_name << ")\n";
    assert(false);
    return kG4UNKNOWN; // unreachable in debug mode
#undef MATCH_PROCESS
#undef MATCH_PROCESS_NAMED

}//GetG4ProcessID

#endif
