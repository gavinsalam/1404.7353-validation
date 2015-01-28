///////////////////////////////////////////////////////////////////////////////////////////////////
//
// This routine is meant to provide an example implementation of some
// of analysis elements described in arXiv:1404.7353
//
// Run it with
//  ./baseline -chs -massless -hard ??? -pileup ???
//
//  ./baseline -massless -chs -hard events/lhc14-pythia8-4C-Zprime500-noUE-nev1e5.pu14.gz -pileup events/lhc14-pythia8-4C-minbias-nev1e6.pu14.gz -nev 300 -npu 100
//
//
// Possible improvements to framework?
// - make it possible to reuse PU events if hard event rejected?
// - make it possible to set massless & chs from C++ interface?
// 
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "EventMixer.hh"
#include "CmdLine.hh"
#include "PU14.hh"
#include "AverageAndError.hh"
#include "CorrelationCoefficient.hh"
#include "ProfileHist.hh"
#include "MatchableSet.hh"

#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/contrib/JetCleanser.hh"
#include "SafeSubtractor.hh" // only for NpC.

#include <iomanip>      // std::setprecision
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>

#include "boost/foreach.hpp"
#define foreach BOOST_FOREACH

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;


/// class for maintaining various comparisons between a reference and
/// "alternate" value (e.g. between the original pt in the hard event
/// and the pt after pileup addition and subtraction).
///
/// Currently maintains the average offset and dispersion, as well as
/// the correlation coefficient.
class Comparison {
public:
  Comparison() {}
  // add an entry, where x is the reference value and y is the quantity being compared to x
  void add_entry(double x, double y) {
    correl.add_entry(x,y);
    averr.add_entry(y-x);
  }
  AverageAndError averr;
  CorrelationCoefficient correl;
};

ostream & operator<<(ostream & ostr, const Comparison & comp) {
  ostr << setprecision(3);
  int w = 6;
  ostr << "<diff> = "  << setw(w) << comp.averr.average()
       <<      " +- "  << setw(w) << comp.averr.error() << "; ";
  ostr << "std.dev = " << setw(w) << comp.averr.sd()
       <<       " +- " << setw(w) << comp.averr.error_on_sd() << "; ";
  ostr << "correl.coeff = " << setw(w) << comp.correl.r() << "; ";
  ostr << "nentries = " << comp.averr.n_entries();
  return ostr;
}

// pre-declaration of a few helpers (implementation at the end of this file)
void do_output(int iev, int & iev_write_interval);
void do_output(int iev);
void do_comparison(const string & name, const vector<PseudoJet> & hard, const vector<PseudoJet> & full_or_corrected);
void print_jets(const string & event_tag, const vector<PseudoJet> & jets, const unsigned int njets_to_print);
Selector SelectorHasLVCharge();

// our main global types
AverageAndError npu, pass_fraction;           //< statistics of npu and matched jets
map<string,Comparison> comp;                  //< stores most of the comparison results
typedef pair<string,Comparison> comp_it_type; //< a helper to give comparisons a name
ostream * ostr = &cout;                       //< stream to output results
ostringstream header;                         //< header containing information about the run
string filename;                              //< output filename
bool do_chs;                                  //< whether or not we work with CHS events

CmdLine * cmdline_ptr;

////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char ** argv) {
  // get fastjet banner out of the way
  ClusterSequence::print_banner();

  // read options from the command line and initialise the event mixer ---------------

  CmdLine cmdline(argc,argv);
  cmdline_ptr = &cmdline;

  filename = cmdline.value<string>("-out", "");
  header << "# " << cmdline.command_line() << "\n#" << endl;

  // create mixer that will construct events by mixing hard and pileup
  // events read from files given from command line using 
  // -hard hard_events_file(.gz) -pileup pileup_events_file(.gz)
  EventMixer mixer(&cmdline);  
  header << "# event mixer: " << mixer.description() << endl;
  do_chs = (mixer.chs_rescaling_factor() < 1);

  // protection to make sure runs are performed only with the massless option
  if (!cmdline.present("-massless") && !cmdline.present("-force-masses")) {
    cerr << "***** ERROR: this program is only set up to run with massless particles (you must include the \"-massless\" command-line option)" << endl;
    exit(-1);
  }
  
  // inputs read from command line
  int nev = cmdline.value<double>("-nev",1.0);  // first argument: command line option; second argument: default value
  int maxprintout = cmdline.value<int>("-maxprintout",1);  

  // particle and jet selection & R
  double particle_rapmax = cmdline.value<double>("-particle-rapmax", 4.0);
  double jet_rapmax      = cmdline.value<double>("-jet-rapmax", 2.5);
  double jet_ptmin       = cmdline.value<double>("-jet-ptmin",  150);
  double R               = cmdline.value<double>("-R",    1.0);
  double Rsub            = cmdline.value<double>("-Rsub", 0.3);
  double fcut            = cmdline.value<double>("-fcut", 0.0);
  header << "# default fcut = " << fcut << endl;
  bool   untrimmed_ref = cmdline.present("-untrimmed-ref");
  if (fcut > 0) {
    if (untrimmed_ref) header << "# using UNTRIMMED reference for groomed jet comparisons" << endl;
    else               header << "# using TRIMMED reference for groomed jet comparisons" << endl;
  }
  
  // rapidity rescaling of rho
  bool rescale = ! cmdline.present("-norescale");
  header << "# rapidity rescaling for rho = " << rescale << endl;
  
  // make sure there are no unused command-line arguments
  cmdline.assert_all_options_used();

  // initialise the objects needed for the event analysis ----------------------------
  
  // some definitions
  JetDefinition jet_def(antikt_algorithm,R);             // the jet definition
  JetDefinition subjet_def(kt_algorithm,Rsub);           // the subjet definition
  AreaDefinition area_def(active_area_explicit_ghosts, GhostedAreaSpec(particle_rapmax));  // the area definition
  header << "# jet_def: "    << jet_def.description() << endl;
  header << "# subjet_def: " << subjet_def.description() << endl;
  header << "# area_def: "   << area_def.description() << endl;           

  // take particles only within |y|<4
  Selector sel_particles = SelectorAbsRapMax(particle_rapmax);
  header << "# sel_particles: " << sel_particles.description() << endl;

  // select the jets
  Selector sel_jets = SelectorAbsRapMax(jet_rapmax)*SelectorPtMin(jet_ptmin);
  header << "# sel_jets: " << sel_jets.description() << endl;

  // set up (safe) area--median subtraction
  // --------------------------------------
  // define background estimator (grid type, does not need to recluster)
  GridMedianBackgroundEstimator * gmbge
              = new GridMedianBackgroundEstimator(particle_rapmax,0.55);
  // define (and then apply) function for rapidity rescaling of rho.
  // NB These parameters have actually been determined at 13 TeV, but they vary slowly
  // and should therefore also be aproproate for 14 TeV
  FunctionOfPseudoJet<double> * rescaling =
    new BackgroundRescalingYPolynomial(1.1685397, 0, -0.0246807, 0, 5.94119e-05);
  if ( rescale) { gmbge -> set_rescaling_class(rescaling); }
  // define subtractor
  Subtractor areasub(gmbge);
  // set up safety (using native FJ3.1 options)
  areasub.set_safe_mass(true);
  // with specific info needed to get a sensible safe CHS where relevant
  if (do_chs) areasub.set_known_selectors(SelectorIsCharged(), SelectorVertexNumber(0));
  header << "# area subtractor: "  << areasub.description() << endl;

  // set up trimming
  // ---------------
  Filter trimmer_nosub(subjet_def, SelectorPtFractionMin(fcut));
  Filter trimmer(subjet_def, SelectorPtFractionMin(fcut));
  trimmer.set_subtractor(&areasub);
  header << "# trimmer: " << trimmer.description() << endl;
  
  // set up subjet zeroing (including trimming)
  // ------------------------------------------
  Filter subjet_zeroing(subjet_def, SelectorPtFractionMin(fcut) && SelectorHasLVCharge());
  subjet_zeroing.set_subtractor(&areasub);
  header << "# subjet_zeroing: " << subjet_zeroing.description() << endl;

  // set up subjet protected zeroing
  // -------------------------------
  // This is the simple version, whereby "protection" simply means that
  // a subjet that doesn't contain any tracks from the LV is still accepted
  // if its subtracted pt > 20GeV. A full implementation would probably
  // seek to make this threshold a function of rho & sigma and/or the charged pt
  // in the subjet each individual PU vertex (e.g. 10 times the maximum as found
  // across all the vertices).
  Filter subjet_protected(subjet_def,
                           SelectorPtFractionMin(fcut)
                           && (SelectorHasLVCharge() || SelectorPtMin(20.0)) );
  subjet_protected.set_subtractor(&areasub);
  header << "# subjet_protected: " << subjet_protected.description() << endl;

  
  // set up NpC subtraction including the various trimming options
  // ----------------------
  double ga0 = 0.612;
  double epsilon = mixer.chs_rescaling_factor();
  notcontrib::SafeNpCSubtractor npcsub(epsilon*ga0/(1-ga0+epsilon*ga0),
                                       SelectorIsCharged(), SelectorVertexNumber(0));
  header << "# npc subtractor: " << npcsub.description() << endl;

  Filter npc_trimmer(subjet_def, SelectorPtFractionMin(fcut));
  npc_trimmer.set_subtractor(&npcsub);
  header << "# npc trimmer: " << npc_trimmer.description() << endl;
  
  // set up npc with subjet zeroing
  Filter npc_subjet_zeroing(subjet_def, SelectorPtFractionMin(fcut) && SelectorHasLVCharge());
  npc_subjet_zeroing.set_subtractor(&npcsub);
  header << "# npc subjet_zeroing: " << npc_subjet_zeroing.description() << endl;

  // set up npc with protected subjet zeroing
  Filter npc_subjet_protected(subjet_def,
                               SelectorPtFractionMin(fcut)
                               && (SelectorHasLVCharge() || SelectorPtMin(20.0)) );
  npc_subjet_protected.set_subtractor(&npcsub);
  header << "# subjet_protected: " << npc_subjet_protected.description() << endl;
  
  // set up Linear Cleansing
  // -----------------------
  JetCleanser linear_cleanser(subjet_def,
                              JetCleanser::linear_cleansing, JetCleanser::input_nc_separate);
  linear_cleanser.SetLinearParameters(ga0);
  linear_cleanser.SetTrimming(fcut);
  header << "# linear_cleanser: " << linear_cleanser.description() << endl;

  // loop over events ----------------------------------------------------------------
  int iev = 0;
  int iev_write_interval = 10;
  while ( mixer.next_event() && iev < nev ) {
     // increment event number    
     iev++;
     
     // extract particles from event
     vector<PseudoJet> full_event = sel_particles(mixer.particles());
     // assign a unique user index to each particle to facilitate matching
     for (unsigned i = 0; i < full_event.size(); i++) {
       full_event[i].set_user_index(i);
     }

     if ( iev <= maxprintout ) { cerr << "\nEvent " << iev << endl; }
     if ( iev <= maxprintout ) { cerr << "nPU = " << mixer.npu() << endl; }
     npu.add_entry(mixer.npu());
     
     // extract hard event
     vector<PseudoJet> hard_event, pileup_event;
     SelectorIsHard().sift(full_event, hard_event, pileup_event); 

     // cluster hard event only (Note: area is not needed here)
     ClusterSequence cs_hard(hard_event,jet_def);
     //ClusterSequenceArea cs_hard(hard_event,jet_def,area_def);
     // select the jets in the hard event
     vector<PseudoJet> hard_jets = sel_jets(sorted_by_pt(cs_hard.inclusive_jets()));
     if (iev <= maxprintout) print_jets("Hard event", hard_jets, hard_jets.size());
     
     // if no hard jets pass the selection, continue on to the next event
     // (writing results if it's the right time)
     if (hard_jets.size() == 0) {
       pass_fraction += 0;
       do_output(iev, iev_write_interval);
       continue;
     }
     pass_fraction += 1;

     // now get (potentially) trimmed reference jets
     vector<PseudoJet> trimmed_hard_jets;
     if (fcut == 0 || untrimmed_ref) {
       trimmed_hard_jets = hard_jets;
     } else {
       trimmed_hard_jets = trimmer_nosub(hard_jets);
     }

     // cluster full event (hard + pileup)
     ClusterSequenceArea cs_full(full_event,jet_def,area_def);
          
     // get all jets in full event
     vector<PseudoJet> full_jets = sorted_by_pt(cs_full.inclusive_jets());
     if (iev <= maxprintout) print_jets("Full event", full_jets, 4U);

     //----------------------------------------------------------------------
     // then match them
     // ---------------
     MatchableSet matching_full(full_jets);
     vector<PseudoJet> matched_full_jets;
     vector<double>    matched_pts, matched_fractions;
     foreach (const PseudoJet &hard_jet, hard_jets){
       matching_full.match_probe(hard_jet, 0);
         
       matched_full_jets.push_back(*(matching_full.matched_jet_ptr()));
       matched_pts.push_back(matching_full.match().match_pt());
       matched_fractions.push_back(matching_full.match().match_pt()/
                                   SelectorIdentity().scalar_pt_sum(hard_jet.constituents()));
     }
     if ( iev <= maxprintout ) { // print it explicitly to include matched fractions
       cerr << "Matched full event" << endl;
       for (unsigned int i=0; i < 4U && i < matched_full_jets.size(); i++) {
         cerr << "  jet " << i << ": "  << matched_full_jets[i]
              << " (matched frac = " << matched_fractions[i] << ")"
              << endl;
       }
     }
     do_comparison("full", hard_jets, matched_full_jets);

     //----------------------------------------------------------------------
     // subtract the matched full jets
     // ------------------------------
     // 1. feed particles to background estimator
     gmbge->set_particles(full_event);
     // 2. subtract
     vector<PseudoJet> subtracted_jets = areasub(matched_full_jets); 
     if (iev <= maxprintout) print_jets("Subtracted event", subtracted_jets, subtracted_jets.size());
     do_comparison("area-sub", hard_jets, subtracted_jets);

     //----------------------------------------------------------------------
     // now the trimmed jets
     // --------------------
     // Note we've used the matched full jets, which means that fcut
     // multiplied the unsubtracted [CHS] pt; this was the procedure
     // adopted in 1404.7353.  We could have passed subtracted jets
     // and then fcut would multiply the subtracted pt
     //
     // We handle all three options: plain trimming, trimming+zeroing
     // and trimming+protected zeroing
     do_comparison("area-trimmed",     trimmed_hard_jets, trimmer(matched_full_jets));
     do_comparison("area-subj-zeroed", trimmed_hard_jets, subjet_zeroing(matched_full_jets));
     do_comparison("area-subj-zprot",  trimmed_hard_jets, subjet_protected(matched_full_jets));
     
     //----------------------------------------------------------------------
     // Now some NpC options
     // --------------------
     do_comparison("npc-sub",         hard_jets,         npcsub(matched_full_jets)); 
     do_comparison("npc-trimmed",     trimmed_hard_jets, npc_trimmer(matched_full_jets)); 
     do_comparison("npc-subj-zeroed", trimmed_hard_jets, npc_subjet_zeroing(matched_full_jets)); 
     do_comparison("npc-subj-zprot",  trimmed_hard_jets, npc_subjet_protected(matched_full_jets)); 

     //----------------------------------------------------------------------
     // now get the linear cleansed jets (only if we have CHS)
     // --------------------
     if (do_chs) {
       vector<PseudoJet> lin_cleansed_jets;
       foreach (const PseudoJet & jet, matched_full_jets) {
         // cleansing doesn't need ghosts, so get rid of them
         vector<PseudoJet> constituents = (!SelectorIsPureGhost())(jet.constituents());
         // sift the constituents into the neutrals and different kinds of tracks
         vector<PseudoJet> neutrals_all, tracks_all, tracks_lv, tracks_pu;
         SelectorIsCharged().sift(constituents, tracks_all, neutrals_all);
         SelectorVertexNumber(0).sift(tracks_all, tracks_lv, tracks_pu);
         // PU tracks had been scaled down (CHS), so scale them back up
         foreach (PseudoJet & track, tracks_pu) {track /= mixer.chs_rescaling_factor();}
         // and run linear cleansing to get a cleansed jet
         lin_cleansed_jets.push_back(linear_cleanser(neutrals_all, tracks_lv, tracks_pu));
       }
       if ( iev <= maxprintout ) print_jets("Cleansed event", lin_cleansed_jets, lin_cleansed_jets.size());
       do_comparison("linclean", trimmed_hard_jets, lin_cleansed_jets);
     }

     // write intermediate results
     do_output(iev, iev_write_interval);

  }  // end loop over events
  // output something only if we ran successfully...
  if (iev > 0) do_output(iev);
  
  //----------------------------------------------------------------------
  // do the output

  
  // free allocated memory
  delete rescaling;
  delete gmbge;
}

//----------------------------------------------------------------------
void do_output(int iev, int & iev_write_interval) {
  if (iev%iev_write_interval == 0) {
    cerr << "iev = " << iev << endl;
    do_output(iev);
    if (iev/iev_write_interval >= 25) iev_write_interval *= 10;
  }
}

void do_output(int iev) {
  ostream * ostr;
  if (filename == "") {
    ostr = &cout;
  } else {
    ostr = new ofstream(filename.c_str());
  }
  *ostr << header.str() << endl;
  *ostr << "# iev = " << iev << " at " << cmdline_ptr->time_stamp() << endl;
  *ostr << "# <npu> = " << npu.average() << endl;
  *ostr << "# fraction of events with (njet >= 1) = "
        << pass_fraction.average() << " +- " << pass_fraction.error() << endl;
  // loop over all our observables
  foreach (const comp_it_type & comp_it, comp) {
    *ostr << setw(30) << left << comp_it.first+": " << right << comp_it.second << endl;
  }
  if (ostr != &cout) delete ostr;
}

//----------------------------------------------------------------------
void do_comparison(const string & name,
                   const vector<PseudoJet> & hard,
                   const vector<PseudoJet> & full_or_corrected) {
  
  assert(hard.size() == full_or_corrected.size());
  for (unsigned i = 0; i < hard.size(); i++) {
    comp["jet.pt."+name ].add_entry(hard[i].pt(), full_or_corrected[i].pt());
    comp["jet.mass."+name ].add_entry(hard[i].m(), full_or_corrected[i].m());
  }
  if (hard.size() >= 2) {
    double hard_dijet_mass = (hard[0] + hard[1]).m();
    double full_dijet_mass = (full_or_corrected[0] + full_or_corrected[1]).m();
    comp["dijet.mass."+name].add_entry(hard_dijet_mass, full_dijet_mass);
  }
}

//----------------------------------------------------------------------
void print_jets(const string & event_tag,
                const vector<PseudoJet> & jets,
                const unsigned int njets_to_print){
  unsigned int nmax = njets_to_print;
  if (jets.size() < njets_to_print) nmax = jets.size();

  cerr << event_tag << endl; 
  for (unsigned int i=0; i < nmax; i++) {
    cerr << "  jet " << i << ": " << jets[i] << endl;
  }
}

//----------------------------------------------------------------------
// for selecting (sub)jets that have at least one charged track from the
// LV
class SelectorWorkerHasLVCharge : public SelectorWorker {
public:
  virtual bool pass(const PseudoJet & jet) const {
    int nLVCharge = (SelectorIsCharged() && SelectorVertexNumber(0)).count(jet.constituents());
    return nLVCharge > 0;
  }
  virtual string description() const {return "has charge from leading vertex";}
};
Selector SelectorHasLVCharge() {return new SelectorWorkerHasLVCharge();}
