# ./baseline -massless -chs -hard /home/ipht/gsoyez/software/extra/2014PileupWorkshop/events/lhc14-pythia8-4C-Zprime500-nobdecay-noUE-nev1e5.pu14.gz -pileup /home/ipht/gsoyez/software/extra/2014PileupWorkshop/events/lhc14-pythia8-4C-minbias-nev1e7.pu14.gz -out ref-results/Zp500-npu105.res -nev 100000 -maxprintout 3 -npu 105 
#
# event mixer: Event mixer using hard events from /home/ipht/gsoyez/software/extra/2014PileupWorkshop/events/lhc14-pythia8-4C-Zprime500-nobdecay-noUE-nev1e5.pu14.gz and 105 pileup events from /home/ipht/gsoyez/software/extra/2014PileupWorkshop/events/lhc14-pythia8-4C-minbias-nev1e7.pu14.gz with CHS (rescaling charged PU by a factor 1e-60) and massless particles
# default fcut = 0
# rapidity rescaling for rho = 1
# jet_def: Longitudinally invariant anti-kt algorithm with R = 1 and E scheme recombination
# subjet_def: Longitudinally invariant kt algorithm with R = 0.3 and E scheme recombination
# area_def: Active area (explicit ghosts) with ghosts of area 0.00997331 (had requested 0.01), placed up to y = 4, scattered wrt to perfect grid by (rel) 1, mean_ghost_pt = 1e-100, rel pt_scatter =  0.1, n repetitions of ghost distributions =  1
# sel_particles: |rap| <= 4
# sel_jets: (|rap| <= 2.5 * pt >= 150)
# area subtractor: Subtractor that uses the following background estimator to determine rho: GridMedianBackgroundEstimator, with rectangular grid with rapidity extent -4 < rap < 4, tile size drap x dphi = 0.533333 x 0.571199; including mass safety tests; using known vertex selection: is_charged and leading vertex selection: vertex number == 0
# trimmer: Filter with subjet_def = Longitudinally invariant kt algorithm with R = 0.3 and E scheme recombination, selection pt >= 0* pt_ref, subtractor: Subtractor that uses the following background estimator to determine rho: GridMedianBackgroundEstimator, with rectangular grid with rapidity extent -4 < rap < 4, tile size drap x dphi = 0.533333 x 0.571199; including mass safety tests; using known vertex selection: is_charged and leading vertex selection: vertex number == 0
# subjet_zeroing: Filter with subjet_def = Longitudinally invariant kt algorithm with R = 0.3 and E scheme recombination, selection (pt >= 0* pt_ref && has charge from leading vertex), subtractor: Subtractor that uses the following background estimator to determine rho: GridMedianBackgroundEstimator, with rectangular grid with rapidity extent -4 < rap < 4, tile size drap x dphi = 0.533333 x 0.571199; including mass safety tests; using known vertex selection: is_charged and leading vertex selection: vertex number == 0
# subjet_protected: Filter with subjet_def = Longitudinally invariant kt algorithm with R = 0.3 and E scheme recombination, selection (pt >= 0* pt_ref && (has charge from leading vertex || pt >= 20)), subtractor: Subtractor that uses the following background estimator to determine rho: GridMedianBackgroundEstimator, with rectangular grid with rapidity extent -4 < rap < 4, tile size drap x dphi = 0.533333 x 0.571199; including mass safety tests; using known vertex selection: is_charged and leading vertex selection: vertex number == 0
# npc subtractor: Neutral-proportional-to-Charged subtraction using known fraction: 1.57732e-60, known vertex selection: is_charged, leading vertex selection: vertex number == 0 and safety checks for pt and mass.
# npc trimmer: Filter with subjet_def = Longitudinally invariant kt algorithm with R = 0.3 and E scheme recombination, selection pt >= 0* pt_ref, subtractor: Neutral-proportional-to-Charged subtraction using known fraction: 1.57732e-60, known vertex selection: is_charged, leading vertex selection: vertex number == 0 and safety checks for pt and mass.
# npc subjet_zeroing: Filter with subjet_def = Longitudinally invariant kt algorithm with R = 0.3 and E scheme recombination, selection (pt >= 0* pt_ref && has charge from leading vertex), subtractor: Neutral-proportional-to-Charged subtraction using known fraction: 1.57732e-60, known vertex selection: is_charged, leading vertex selection: vertex number == 0 and safety checks for pt and mass.
# subjet_protected: Filter with subjet_def = Longitudinally invariant kt algorithm with R = 0.3 and E scheme recombination, selection (pt >= 0* pt_ref && (has charge from leading vertex || pt >= 20)), subtractor: Neutral-proportional-to-Charged subtraction using known fraction: 1.57732e-60, known vertex selection: is_charged, leading vertex selection: vertex number == 0 and safety checks for pt and mass.
# linear_cleanser: JetCleanser [Linear mode, input = neutral and charged separate]
 Trimming: fcut = 0
 g0_mean = 0.612


# iev = 95238 at 2015-01-20 13:34:35 (CET)
# <npu> = 105
# fraction of events with (njet >= 1) = 0.756799 +- 0.00139018
dijet.mass.area-sub:          <diff> =  -1.07 +-  0.143; std.dev =   33.5 +-  0.385; correl.coeff =  0.967; nentries = 54495
dijet.mass.area-subj-zeroed:  <diff> =   -3.4 +-  0.185; std.dev =   43.3 +-  0.859; correl.coeff =  0.947; nentries = 54495
dijet.mass.area-subj-zprot:   <diff> =  0.768 +-  0.107; std.dev =   24.9 +-  0.479; correl.coeff =  0.982; nentries = 54495
dijet.mass.area-trimmed:      <diff> =     35 +-   0.13; std.dev =   30.3 +-  0.419; correl.coeff =  0.973; nentries = 54495
dijet.mass.full:              <diff> =    297 +-  0.302; std.dev =   70.4 +-  0.452; correl.coeff =  0.902; nentries = 54495
dijet.mass.linclean:          <diff> =   5.89 +-   0.19; std.dev =   44.5 +-   0.84; correl.coeff =  0.945; nentries = 54495
dijet.mass.npc-sub:           <diff> = -0.171 +-  0.141; std.dev =     33 +-  0.382; correl.coeff =  0.968; nentries = 54495
dijet.mass.npc-subj-zeroed:   <diff> = -0.903 +-  0.185; std.dev =   43.2 +-  0.848; correl.coeff =  0.947; nentries = 54495
dijet.mass.npc-subj-zprot:    <diff> =   3.09 +-  0.106; std.dev =   24.8 +-   0.48; correl.coeff =  0.982; nentries = 54495
dijet.mass.npc-trimmed:       <diff> =   43.9 +-  0.126; std.dev =   29.5 +-  0.427; correl.coeff =  0.975; nentries = 54495
jet.mass.area-sub:            <diff> =   1.54 +- 0.0534; std.dev =   19.2 +-  0.103; correl.coeff =  0.772; nentries = 129865
jet.mass.area-subj-zeroed:    <diff> =   1.23 +- 0.0368; std.dev =   13.3 +-  0.143; correl.coeff =  0.884; nentries = 129865
jet.mass.area-subj-zprot:     <diff> =   1.73 +- 0.0364; std.dev =   13.1 +-  0.144; correl.coeff =  0.887; nentries = 129865
jet.mass.area-trimmed:        <diff> =   22.6 +- 0.0474; std.dev =   17.1 +-  0.112; correl.coeff =   0.79; nentries = 129865
jet.mass.full:                <diff> =    110 +- 0.0656; std.dev =   23.6 +- 0.0995; correl.coeff =  0.631; nentries = 129865
jet.mass.linclean:            <diff> =   6.91 +- 0.0362; std.dev =     13 +-  0.146; correl.coeff =  0.888; nentries = 129865
jet.mass.npc-sub:             <diff> =   2.18 +- 0.0546; std.dev =   19.7 +-  0.104; correl.coeff =  0.767; nentries = 129865
jet.mass.npc-subj-zeroed:     <diff> =   2.77 +- 0.0372; std.dev =   13.4 +-  0.142; correl.coeff =  0.882; nentries = 129865
jet.mass.npc-subj-zprot:      <diff> =   3.19 +- 0.0365; std.dev =   13.1 +-  0.144; correl.coeff =  0.887; nentries = 129865
jet.mass.npc-trimmed:         <diff> =     27 +- 0.0468; std.dev =   16.9 +-  0.115; correl.coeff =  0.793; nentries = 129865
jet.pt.area-sub:              <diff> = -0.468 +- 0.0422; std.dev =   15.2 +-  0.117; correl.coeff =  0.963; nentries = 129865
jet.pt.area-subj-zeroed:      <diff> =  -1.03 +- 0.0482; std.dev =   17.4 +-  0.259; correl.coeff =  0.952; nentries = 129865
jet.pt.area-subj-zprot:       <diff> =  0.184 +- 0.0306; std.dev =     11 +-  0.156; correl.coeff =   0.98; nentries = 129865
jet.pt.area-trimmed:          <diff> =   11.2 +- 0.0364; std.dev =   13.1 +-  0.136; correl.coeff =  0.972; nentries = 129865
jet.pt.full:                  <diff> =    102 +-  0.059; std.dev =   21.3 +-  0.104; correl.coeff =  0.931; nentries = 129865
jet.pt.linclean:              <diff> =   2.15 +- 0.0503; std.dev =   18.1 +-  0.264; correl.coeff =  0.948; nentries = 129865
jet.pt.npc-sub:               <diff> = -0.313 +- 0.0427; std.dev =   15.4 +-  0.115; correl.coeff =  0.962; nentries = 129865
jet.pt.npc-subj-zeroed:       <diff> = -0.162 +- 0.0486; std.dev =   17.5 +-  0.258; correl.coeff =  0.951; nentries = 129865
jet.pt.npc-subj-zprot:        <diff> =  0.999 +-  0.031; std.dev =   11.2 +-  0.155; correl.coeff =  0.979; nentries = 129865
jet.pt.npc-trimmed:           <diff> =   14.2 +- 0.0356; std.dev =   12.8 +-  0.139; correl.coeff =  0.973; nentries = 129865
