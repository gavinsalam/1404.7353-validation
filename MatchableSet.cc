#include "MatchableSet.hh"
#include "fastjet/ClusterSequence.hh"

using namespace fastjet;
using namespace std;


//--- stolen from fastjet code: should be made accessible directly from there? ---
/// given a vector of values with a one-to-one correspondence with the
/// vector of objects, sort objects into an order such that the
/// associated values would be in increasing order
template<class T> vector<T>  objects_sorted_by_values_cp(
                       const vector<T> & objects, 
		       const vector<double> & values) {

  assert(objects.size() == values.size());

  // get a vector of indices
  vector<int> indices(values.size());
  for (size_t i = 0; i < indices.size(); i++) {indices[i] = i;}
  
  // sort the indices
  sort_indices(indices, values);
  
  // copy the objects 
  vector<T> objects_sorted(objects.size());
  
  // place the objects in the correct order
  for (size_t i = 0; i < indices.size(); i++) {
    objects_sorted[i] = objects[indices[i]];
  }

  return objects_sorted;
}

vector<SharedPtr<PseudoJet> > sorted_by_pt(const vector<SharedPtr<PseudoJet> > & jets) {
  vector<double> minus_kt2(jets.size());
  for (size_t i = 0; i < jets.size(); i++) {minus_kt2[i] = -jets[i]->kt2();}
  return objects_sorted_by_values_cp(jets, minus_kt2);
}

//----------------------------------------------------------------------
/// return a vector of jets sorted into increasing user index
vector<PseudoJet> sorted_by_user_index(const vector<PseudoJet> & jets) {
  vector<double> user_index(jets.size());
  for (size_t i = 0; i < jets.size(); i++) {user_index[i] = jets[i].user_index();}

  return objects_sorted_by_values_cp(jets, user_index);
}


Match::Match(const PseudoJet & ref_jet, const PseudoJet & probe_jet) {
  vector<PseudoJet> cns_ref = sorted_by_user_index(ref_jet.constituents());
  vector<PseudoJet> cns_prb = sorted_by_user_index(probe_jet.constituents());
  
  match(cns_ref, cns_prb);
}

/// establish the match given the constituents (sorted in index order)
void Match::match(const vector<PseudoJet> & cns_ref, 
		  const vector<PseudoJet> & cns_prb) {
  _match_pt = 0;
  vector<PseudoJet>::const_iterator ref = cns_ref.begin();
  vector<PseudoJet>::const_iterator prb = cns_prb.begin();

  while (ref != cns_ref.end() && prb != cns_prb.end()) {
    if      (ref->user_index() < 0) {ref++;}
    else if (prb->user_index() < 0) {prb++;}
    else if (ref->user_index() < prb->user_index()) {ref++;}
    else if (ref->user_index() > prb->user_index()) {prb++;}
    else {
      // we have a match
      _match_pt += ref->perp();
      ref++;
      prb++;
    }
  }

  // get the reference pt
  double ref_pt = 0;
  for (ref = cns_ref.begin(); ref != cns_ref.end(); ref++) {
    if (ref->user_index() >= 0) ref_pt += ref->perp();
  }

  // set the fraction unless ref_pt = 0;
  _match_fraction = ref_pt > 0 ? _match_pt/ref_pt : 0;
}

/// Establish the back reaction for this jet.
///
/// The set of particles over which back-reaction is calculation is
/// determined as follows:
/// 
/// - take the first reference constituent 
/// - find its cluster sequence
/// - take as the range the user index of the first input particle to
///   that of the last input particle.
/// - check it corresponds to the number of particles!
///   (is this necessary?)
/// 
/// It is assumed that all reference constituents belong to
/// the same cluster sequence
void Match::back_reaction(const vector<PseudoJet> & cns_ref, 
			  const vector<PseudoJet> & cns_prb,
			  double & br_gain, 
			  double & br_loss
			  ) const {
  // get the reference cluster sequence
  assert(cns_ref.size() >0);
  assert(cns_ref[0].has_associated_cluster_sequence());
  const ClusterSequence * ref_cs = cns_ref[0].associated_cluster_sequence();

  // figure out a valid user-index range
  int ref_index_lo = ref_cs->jets()[0].user_index();
  unsigned i;
  for (i = 0; i < ref_cs->n_particles(); i++) {
    if (ref_cs->jets()[i].user_index() < 0) break;
  }
  assert(i>0);
  int ref_index_hi = ref_cs->jets()[i-1].user_index();
  assert(ref_index_hi - ref_index_lo == int(i)-1);

  br_gain = 0.0;
  br_loss = 0.0;

  vector<PseudoJet>::const_iterator ref = cns_ref.begin();
  vector<PseudoJet>::const_iterator prb = cns_prb.begin();

  while (ref != cns_ref.end() && prb != cns_prb.end()) {
    if      (ref->user_index() < 0) {ref++;}
    else if (prb->user_index() < 0) {prb++;}
    else if (ref->user_index() < prb->user_index()) {
      br_loss += ref->perp();
      ref++;
    }
    else if (ref->user_index() > prb->user_index()) {
      if (ref_index_lo <= prb->user_index() && prb->user_index() <= ref_index_hi) 
	br_gain += prb->perp();
      prb++;
    }
    else {
      // we have a match
      ref++;
      prb++;
    }
  }
}



//======================================================================
MatchableSet::MatchableSet(const std::vector<fastjet::PseudoJet> & ref_jets)
{
  _ref_jets.resize(ref_jets.size());
  _cns_ref.resize(ref_jets.size());
  // specific way of writing to allow for future change over to pointers?
  for (unsigned iref = 0; iref < _ref_jets.size(); iref++) {
    _ref_jets[iref] = ref_jets[iref];
    _cns_ref[iref] = sorted_by_user_index(ref_jets[iref].constituents());
  }
}

//======================================================================
MatchableSet::MatchableSet(const vector<SharedPtr<PseudoJet> >& ref_jets)
{
  _ref_jets.resize(ref_jets.size());
  _cns_ref.resize(ref_jets.size());
  // specific way of writing to allow for future change over to pointers?
  for (unsigned iref = 0; iref < _ref_jets.size(); iref++) {
    _ref_jets[iref] = *(ref_jets[iref]);
    _cns_ref[iref] = sorted_by_user_index(ref_jets[iref]->constituents());
  }
}


//----------------------------------------------------------------------
void MatchableSet::match_probe(const fastjet::PseudoJet & prb_jet, double min_fraction) {

  _back_reaction_calculated = false;
  _cns_prb = sorted_by_user_index(prb_jet.constituents());

  Match best_match;
  _matched_jet = 0;
  _matched_iref = -1;
  for (unsigned iref = 0; iref < _ref_jets.size(); iref++) {
    Match match(_cns_ref[iref], _cns_prb);
    if (match.match_fraction() > min_fraction &&
	match.match_pt() > best_match.match_pt()) {
      best_match = match;
      _matched_iref = iref;
      //_matched_jet = & _ref_jets[iref];
    }
  }
  _match = best_match;

  if (_matched_iref >= 0) _matched_jet = & _ref_jets[_matched_iref];
  else                    _matched_jet = 0;
}

//----------------------------------------------------------------------
void MatchableSet::calculate_bach_reaction() const {
  _match.back_reaction(_cns_ref[_matched_iref], _cns_prb, _br_gain, _br_loss);
  _back_reaction_calculated = true;
}
