// $Id: SafeSubtractor.hh 630 2014-04-30 21:40:31Z gsoyez $
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

#ifndef __SAFE_SUBTRACTOR_HH__
#define __SAFE_SUBTRACTOR_HH__

// open questions
//
// - AreaSubtractor::use_common_bge_for_rho_and_rhom()?
// - examples: do we keep the current one and add a "massive" example
//             do we add an example_chs.cc?

#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/Transformer.hh"
#include "fastjet/tools/BackgroundEstimatorBase.hh"

FASTJET_BEGIN_NAMESPACE

namespace notcontrib {

//----------------------------------------------------------------------
// The main tools introduced in this file are the 
//    SafeAreaSubtractor  and    SafeNpCSubtractor
// classes. Their declaration and description is in this file, right
// after the specification of their common base class SafeSubtractorBase
//----------------------------------------------------------------------


//----------------------------------------------------------------------
/// \class SafeSubtractorBase
/// Base class for a series of subtraction methods
///
/// This class will split the jet constituents in 2 groups: the ones
/// with known vertex origin and the ones with unknown vertex
/// origin. The ones with known vertex origin are then further split
/// into leading-vertex and PU vertices particles.
///
/// Safety constraints are applied on the result of subtraction: the
/// subtracted pt and, optionally, the subtracted mass mass needs to
/// be at least the pt and mass we know to come from the leading
/// vertex. A few specific cases need special care:
///   - when the subtracted pt goes below the pt known to come from
///     the leading vertex, we set the pt, rapidity and phi to the
///     known contribution from the leading vertex.
///   - when the subtracted jet mass goes below the mass known to come
///     from the leading vertex, we set the pt and azimuthal angle to
///     the ones of the subtracted jet, the mass and rapidity to
///     corresponding values from the leading vertex contribution
///   - If rapidity and phi are ill-defined, we use the one from the
///     original, unsubtracted jet.
///
/// Derived classes should then overload '_amount_to_subtract' which
/// computes the amount to be subtracted from the full jet based on
/// the three momenta above (note that area information is attached to
/// the 'unknown' component).
///
/// When 'safe_mass' is set to false, the safety tests on the
/// subtracted jet mass described above are not performed (i.e. the
/// mass is allowed to go below the mass known to be associated with
/// the leading vertex)
///
/// By default, the Selectors are such that all the particles are of
/// unknown vertex origin, i.e. things behave as for a "full"
/// event. If one wants to apply this to a CHS event, one would
/// typically use a Selector that keeps charged particles as the
/// known-vertex selector (we know the vertex of origin for the
/// charged tracks but not for the neutrals) and a Selector that keeps
/// particles from the leading vertex for selecting, among the charged
/// tracks, the ones originating from the leading vertex.
class SafeSubtractorBase : public Transformer{
public:
  /// default ctor:
  ///   \param selector_known_vertex    selects particles of known vertex
  ///                                   (leading or pileup)
  ///   \param selector_leading_vertex  among particles of known
  ///                                   vertex, select the ones coming
  ///                                   from the leading vertex
  ///   \param safe_mass   impose that the mass of the resulting jet is
  ///                      not smaller than the known contribution from
  ///                      the leading vertex
  SafeSubtractorBase(Selector selector_known_vertex   = !SelectorIdentity(),
		     Selector selector_leading_vertex = SelectorIdentity(),
		     bool safe_mass=true)
    : _selector_known_vertex(selector_known_vertex),
      _selector_leading_vertex(selector_leading_vertex),
      _safe_mass(safe_mass){}

  /// the result of the subtraction
  ///   \param jet   the jet to be subtracted
  virtual PseudoJet result(const PseudoJet &jet) const;

  /// a reminder from the base class: one has to provide a description!
  virtual std::string description() const = 0;

protected: 
  /// compute the amount to be subtracted from the various bits in the jet:
  ///   \param known_lv  momentum known as coming from the leading vertex
  ///   \param known_pu  momentum known as coming from the pu vertices
  ///   \param unknown   momentum of unknown vertex origin
  ///
  /// each derived class will have to overload this
  virtual PseudoJet _amount_to_subtract(const PseudoJet & known_lv, 
                                        const PseudoJet & known_pu, 
                                        const PseudoJet & unknown) const = 0;

  Selector _selector_known_vertex;   ///< constituents of known vertex
  Selector _selector_leading_vertex; ///< constituents from the leading vertex
  bool _safe_mass;                   ///< pt and mass safety

private:
  /// given a selector and a vector of PseudoJet, this returns the
  /// 4-vector sum of the jets passing the selection creteria. If the
  /// vector is empty, return a 0-momentum PSeudoJet with structural
  /// info inherited from 'jet'
  PseudoJet _sum(const std::vector<PseudoJet> &particles,
		 const PseudoJet &jet) const;
};

//----------------------------------------------------------------------
/// \class SafeAreaSubtractor
///
/// A Transformer that takes an input jet contaminated by pileup,
/// returns a subtracted jet. 
///
/// Compared to FastJet's native Subtractor, its functionality differs
/// in two ways:
///
/// 1) it also includes the extra correction for massive particles
///    (rho_m in arXiv:1211.2811)
///
/// 2) in addition to imposing that the subtracted jet pt remains
///    positive, it also optionally ensures that the subtracted jet
///    has a non-negative mass.
///
/// 3) it also allows for the treatment of CHS-type of events, where
///    subtraction should only be applied to the part of the event
///    with no vertexing information. In this case, the pt and
///    optionally the mass of the subtracted jet are prevented from
///    going smaller than the part that is known for sure to come from
///    the leading vertex. See the SafeSubtractorBase class for more
///    details.
///
class SafeAreaSubtractor : public SafeSubtractorBase{
public:
  /// default ctor
  /// 
  /// Parameters specific to area subtraction :
  ///   \param bge    a (pointer to) the BackgroundEstimator used to
  ///                 estimate rho
  ///   \param bge_m  a (pointer to) the BackgroundEstimator used to
  ///                 estimate rho_m [the contriobution from the
  ///                 particle masses]. It should support density
  ///                 classes.
  ///
  /// Parameters inherited from the base class:
  ///   \param selector_known_vertex    selects particles of known vertex
  ///                                   (leading or pileup)
  ///   \param selector_leading_vertex  among particles of known
  ///                                   vertex, select the ones coming
  ///                                   from the leading vertex
  ///   \param safe_mass   impose that the mass of the resulting jet is
  ///                      not smaller than the known contribution from
  ///                      the leading vertex
  /// See the SafeSubtractorBase declaration above for more details.
  SafeAreaSubtractor(fastjet::BackgroundEstimatorBase * bge,
		     fastjet::BackgroundEstimatorBase * bge_m = 0,
		     Selector selector_known_vertex   = !SelectorIdentity(),
		     Selector selector_leading_vertex = SelectorIdentity(),
		     bool safe_mass=true)
    : SafeSubtractorBase(selector_known_vertex, selector_leading_vertex, safe_mass),
      _bge(bge), _bge_m(bge_m){}

  /// class description
  virtual std::string description() const;

protected:
  /// compute the amount to be subtracted from the various bits in the jet:
  ///   \param known_lv  momentum known as coming from the leading vertex
  ///   \param known_pu  momentum known as coming from the pu vertices
  ///   \param unknown   momentum of unknown vertex origin
  virtual PseudoJet _amount_to_subtract(const PseudoJet &known_lv,
                                        const PseudoJet &known_pu,
                                        const PseudoJet &unknown) const;

  mutable BackgroundEstimatorBase *_bge;   ///< estimation of the pileup pt
  mutable BackgroundEstimatorBase *_bge_m; ///< estimation of the contribution
                                           ///  from PU particle masses
};

//----------------------------------------------------------------------
/// \class NpCSubtractor
/// a version of SubtractorBase that estimates the amount to subtract
/// from the component of unknown vertex oigin from the contribution
/// that is known to come from PU.
///
/// The proportionality constant is obtained from the average fraction
/// of the PU that has a known vertex of origin. It can be specified
/// either as a constant of as a double-valued function of PseudoJet.
class SafeNpCSubtractor : public SafeSubtractorBase{
public:
  /// ctor with fraction of known PU specified as a constant
  /// 
  /// Parameters specific to area subtraction:
  ///
  ///   \param gamma  the average fraction of PU which is from a
  ///                 known vertex
  ///
  /// Note that if the known-vertex particles account for a fraction,
  /// say, gamma0 of the total, the estimated amount to subtract will
  /// be taken as 1/gamma0 (times the known pu). Now, if known PU
  /// particles have been scaled down by a factor epsilon (e.g. if one
  /// is working with the CHS event), then that fraction of known PU
  /// will be
  ///
  ///    epsilon gamma0 / (1 - gamma0 + epsilon gamma0)
  ///
  /// [GPS: I'm finding this a bit confusing. How about we have a
  /// parameter total_to_known_PU_ratio, and then we give two
  /// examples: if charged fraction is gamma, then without CHS we use
  ///
  ///   total_to_known_PU_ratio = 1/gamma
  ///
  /// Instead, in a CHS context
  ///
  ///   total_to_known_PU_ratio = (1-gamma)/(epsilon gamma) - 1
  ///
  /// ]
  ///
  /// Parameters inherited from the base class:
  ///   \param selector_known_vertex    selects particles of known vertex
  ///                                   (leading or pileup)
  ///   \param selector_leading_vertex  among particles of known
  ///                                   vertex, select the ones coming
  ///                                   from the leading vertex
  ///   \param safe_mass   impose that the mass of the resulting jet is
  ///                      not smaller than the known contribution from
  ///                      the leading vertex
  /// See the SafeSubtractorBase declaration above for more details.
  SafeNpCSubtractor(double gamma,
		    Selector selector_known_vertex   = !SelectorIdentity(),
		    Selector selector_leading_vertex = SelectorIdentity(),
		    bool safe_mass=true)
    : SafeSubtractorBase(selector_known_vertex, selector_leading_vertex, safe_mass),
      _gamma(gamma), _gamma_fct(0){}

  /// ctor with fraction of known PU specified as a function of the
  /// PseudoJet
  ///
  /// See above ctor for a detailed description.
  SafeNpCSubtractor(FunctionOfPseudoJet<double> * gamma_fct,
		    Selector selector_known_vertex   = !SelectorIdentity(),
		    Selector selector_leading_vertex = SelectorIdentity(),
		    bool safe_mass=true)
    : SafeSubtractorBase(selector_known_vertex, selector_leading_vertex, safe_mass),
      _gamma(1.0), _gamma_fct(gamma_fct){}

  /// class description
  virtual std::string description() const;

protected:
  /// compute the amount to be subtracted from the various bits in the jet:
  ///   \param known_lv  momentum known as coming from the leading vertex
  ///   \param known_pu  momentum known as coming from the pu vertices
  ///   \param unknown   momentum of unknown vertex origin
  virtual PseudoJet _amount_to_subtract(const PseudoJet &known_lv,
                                        const PseudoJet &known_pu,
                                        const PseudoJet &unknown) const;
  const double _gamma;   ///< the traction of PU from known vertices
  const FunctionOfPseudoJet<double> *_gamma_fct;  ///< the fraction of PU from known vertices 
                                                  ///< (now expressed as a function of a PJ)
};


} // namespace notcontrib

FASTJET_END_NAMESPACE

#endif  // __SAFE_SUBTRACTOR_HH__
