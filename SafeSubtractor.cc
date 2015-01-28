// $Id: SafeSubtractor.cc 620 2014-04-28 22:36:41Z gsoyez $
//
// Copyright (c) 2014-, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#include "SafeSubtractor.hh"
#include "fastjet/config.h"

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace notcontrib {

//----------------------------------------------------------------------
// SubtractorBase implementation
//----------------------------------------------------------------------

// the result of the subtraction
//   - jet   the jet to be subtracted
PseudoJet SafeSubtractorBase::result(const PseudoJet &jet) const{
  // \todo Notes:
  //
  //  - the whole procedure is overkilling for "regular" areasub: no
  //    need for constituents, no need for all the sifting,...
  
  
  // pre-requirement: the jet should have constituents
  if (! jet.has_constituents()){
    throw Error("SubtractionBase: subtraction can only be applied for jets with constituents");
    return PseudoJet();
  }
  
  // separate the jet constituents in 3 groups:
  //   unknown vertex
  //   known vertex, leading vertex
  //   known vertex, non-leading vertex (PU)
  vector<PseudoJet> constits_unknown, constits_known;
  _selector_known_vertex.sift(jet.constituents(), 
                              constits_known,
                              constits_unknown);
  vector<PseudoJet> constits_known_lv, constits_known_pu;
  _selector_leading_vertex.sift(constits_known,
                                constits_known_lv,
                                constits_known_pu);

  // compute the momentum of these 3 sets of constiutuents
  //
  // For the parts related to the known vertices (LV or PU), we just
  // sum the 4-momenta. For the unknown part, we assign it the full
  // jet area.
  //
  //Selector sel_id = SelectorIdentity();
  //PseudoJet known_lv = (constits_known_lv.size() != 0)
  //  ? sel_id.sum(constits_known_lv) : 0.0*jet;
  //PseudoJet known_pu = (constits_known_pu.size() != 0)
  //  ? sel_id.sum(constits_known_pu) : 0.0*jet;
  //PseudoJet unknown = jet; // that keeps all info including area
  //unknown.reset_momentum((constits_unknown.size() != 0)
  //                       ? sel_id.sum(constits_unknown)  : 0.0*jet);
  //
  PseudoJet known_lv = _sum(constits_known_lv, jet);
  PseudoJet known_pu = _sum(constits_known_pu, jet);
  PseudoJet unknown = jet; // that keeps all info including area
  unknown.reset_momentum(_sum(constits_unknown, jet));

  // compute the 4-momentum that needs to be subtracted given these 3
  // components of the jet
  PseudoJet to_subtract = _amount_to_subtract(known_lv, known_pu, unknown);

  // now apply the subtraction itself (optionally performing safety
  // tests)
  PseudoJet subtracted_jet = jet;

  // for the moment only check that subtraction does not give a -ve
  // pt.
  //
  // Note that the lines below preserve the structural information
  // associated with the original jet
  if (to_subtract.pt2() <= subtracted_jet.pt2()){
    subtracted_jet -= to_subtract;
  } else {
    // we have subtracted more than what we know for sure comes from
    // the leading vertex, so we just return the part we know for sure
    // comes from the leading vertex
    subtracted_jet.reset_momentum(known_lv);

    // skip the safety tests below
    return subtracted_jet;
  }
  
  // now perform safety tests compared to the part we know comes from
  // the leading vertex
  //
  // Note that at this stage 'subtracted_jet' has already its momentum
  // subtracted (and has +ve pt)
  //
  // first check that the pt is at least the pt of the known particles
  if (subtracted_jet.pt2() < known_lv.pt2()){
    // we have subtracted more than what we know for sure comes
    // from the leading vertex [GS: not totally happy with this comment]
    // Note that in this case we can skip the mass test below!
    subtracted_jet.reset_momentum(known_lv);
  } else { // the pt safety test passed...
    // ... now perform the mass test
    if ((_safe_mass) && (subtracted_jet.m2() < known_lv.m2())){
      // in this case, we keep pt and phi as obtained from the
      // subtraction above and take rap and m from 'known_lv'
      subtracted_jet.reset_momentum(PtYPhiM(subtracted_jet.pt(),
					    known_lv.rap(),
					    subtracted_jet.phi(),
					    known_lv.m()));
    }
  }

  return subtracted_jet;
}

// given a vector of PseudoJet, this returns the 4-vector sum of the
// jets passing the selection creteria. If the vector is empty, return
// a 0-momentum PSeudoJet with structural info inherited from 'jet'
PseudoJet SafeSubtractorBase::_sum(const vector<PseudoJet> &particles,
				   const PseudoJet &jet) const{

  // note that in FJ3.1 we could use 
  //  return (particles.size() != 0) ? SelectoIdentity().sum(particles) : 0.0*jet;

  PseudoJet this_sum = 0.0*jet;
  for (unsigned i = 0; i < particles.size(); i++)
    this_sum += particles[i];
  
  return this_sum;
}


//----------------------------------------------------------------------
// AreaSubtractor implementation
//----------------------------------------------------------------------

// class description
string SafeAreaSubtractor::description() const{
  ostringstream oss;
  oss << "Area-median subtraction using background estimation: " << _bge->description();
  if (_bge_m)
    oss << ", background estimation for particle masses: " << _bge_m->description();
  oss << ", known vertex selection: " << _selector_known_vertex.description()
      << ", leading vertex selection: " << _selector_leading_vertex.description();
  oss << " and safety checks for pt";
  if (_safe_mass) oss << " and mass";
  oss << ".";

  return oss.str();
}


// compute the amount to be subtracted from the various bits in the jet:
//   - known_lv  momentum known as coming from the leading vertex
//   - known_pu  momentum known as coming from the pu vertices
//   - unknown   momentum of unknown vertex origin
PseudoJet SafeAreaSubtractor::_amount_to_subtract(const PseudoJet &known_lv,
						  const PseudoJet &known_pu,
						  const PseudoJet &unknown) const{
  // make sure the background estimator is present
  if (!_bge){
    throw Error("AreaSubtractor requires a non-zero background estimator");
    return PseudoJet();
  }

  // the amount to subtract receives 3 contributions:
  //  1. the known pu
  PseudoJet to_subtract = known_pu;

  //  2. the "pt" background from the unknown (obtained from the area)
  PseudoJet area = unknown.area_4vector();
  to_subtract += _bge->rho(unknown) * area;

  // 3. an optional contribution from the unknown particles masses
  if (_bge_m)
    to_subtract += _bge_m->rho(unknown) * PseudoJet(0.0, 0.0, area.pz(), area.E());

  return to_subtract;
}


//----------------------------------------------------------------------
// NpCSubtractor implementation
//----------------------------------------------------------------------

// class description
string SafeNpCSubtractor::description() const{
  ostringstream oss;
  oss << "Neutral-proportional-to-Charged subtraction using known fraction: ";

  // the fraction of known pt
  if (_gamma_fct)
    oss << _gamma_fct->description();
  else
    oss << _gamma;

  oss << ", known vertex selection: " << _selector_known_vertex.description()
      << ", leading vertex selection: " << _selector_leading_vertex.description();
  oss << " and safety checks for pt";
  if (_safe_mass) oss << " and mass";
  oss << ".";

  return oss.str();
}


// compute the amount to be subtracted from the various bits in the jet:
//   - known_lv  momentum known as coming from the leading vertex
//   - known_pu  momentum known as coming from the pu vertices
//   - unknown   momentum of unknown vertex origin
PseudoJet SafeNpCSubtractor::_amount_to_subtract(const PseudoJet &known_lv,
						 const PseudoJet &known_pu,
						 const PseudoJet &unknown) const{
  // the amount to subtract receives 2 contributions:
  //  1. the known pu
  //  2. the unknown PU obtained by scaling the known PU
  //
  // This is  known_pu + (1-g)/g known_pu = 1/g known_pu
  double g = (_gamma_fct) ? (*_gamma_fct)(known_pu) : _gamma;
  return 1.0/g * known_pu;
}


} // namespace notcontrib

FASTJET_END_NAMESPACE
