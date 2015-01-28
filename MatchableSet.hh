#ifndef __MATCHABLESET_HH__
#define __MATCHABLESET_HH__
#include "fastjet/PseudoJet.hh"

//
/// Class to help find the matched pt between the probe jet and the
/// reference jet both in absolute terms, and as a fraction of the
/// reference jet scalar pt sum;
///
/// matching is based on coincidence of user indices, and only particles
/// with user_index() > 0 are considered
class Match {
public:
  Match() : _match_fraction(0), _match_pt(0) {}

  /// just establishes match between reference and probe jets
  Match(const fastjet::PseudoJet & ref_jet_in, const fastjet::PseudoJet & probe_jet_in);

  /// establishes the match given two vectors of constituents, assumed to have
  /// already been sorted according to constituent index order
  Match(const std::vector<fastjet::PseudoJet> & cns_ref,
	const std::vector<fastjet::PseudoJet> & cns_probe) {
    match(cns_ref, cns_probe);
  }

  /// establishes the match given two vectors of constituents, assumed to have
  /// already been sorted according to constituent index order
  void match(const std::vector<fastjet::PseudoJet> & cns_ref,
	     const std::vector<fastjet::PseudoJet> & cns_probe);

  void back_reaction(const std::vector<fastjet::PseudoJet> & cns_ref, 
		     const std::vector<fastjet::PseudoJet> & cns_prb,
		     double & br_gain, 
		     double & br_loss
		     ) const;

  /// initialise
  double match_fraction() const {return _match_fraction;}
  double match_pt() const {return _match_pt;}

private:
  double _match_fraction, _match_pt;
};


//----------------------------------------------------------------------
/// class which is initialised with a set of reference jets and which
/// is designed to make repeated matching checks with one or more probe
/// jets a little more efficient...
class MatchableSet {
public:
  MatchableSet(const std::vector<fastjet::PseudoJet> & ref_jets);
  MatchableSet(const std::vector<fastjet::SharedPtr<fastjet::PseudoJet> > & ref_jets);
  
  /// establish the best match among the reference jets for the
  /// supplied probe, ignoring cases where less than 50% of the
  /// reference jet's pt is contained in the probe.
  void match_probe(const fastjet::PseudoJet & probe_jet, double min_fraction = 0.5);

  /// return the class containing the matching information (match
  /// fraction and matched pt)
  const Match & match() const {return _match;}
  /// return a pointer to the matched (reference) jet
  const fastjet::PseudoJet * matched_jet_ptr() const {return _matched_jet;}

  /// work out the back reaction gain & loss for the current match
  /// (restricting the count to the reference jet's CS)
  void calculate_bach_reaction() const;
  
  /// returns the back reaction gain
  double back_reaction_gain() const {
    if (!_back_reaction_calculated) calculate_bach_reaction();
    return _br_gain;
  }

  /// returns the back reaction loss
  double back_reaction_loss() const {
    if (!_back_reaction_calculated) calculate_bach_reaction();
    return _br_loss;
  }

  /// returns the net back reaction 
  double back_reaction() const {
    if (!_back_reaction_calculated) calculate_bach_reaction();
    return _br_gain-_br_loss;
  }

private:
  Match _match;
  fastjet::PseudoJet * _matched_jet;
  std::vector<std::vector<fastjet::PseudoJet> > _cns_ref;
  std::vector<fastjet::PseudoJet>  _cns_prb;
  std::vector<fastjet::PseudoJet>  _ref_jets;
  int _matched_iref;
  mutable bool _back_reaction_calculated;
  mutable double _br_gain, _br_loss;
  
};


std::vector<fastjet::SharedPtr<fastjet::PseudoJet> > 
sorted_by_pt(const std::vector<fastjet::SharedPtr<fastjet::PseudoJet> > & jets);

#endif // __MATCHABLESET_HH__

