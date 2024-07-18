// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "Common/Core/RecoDecay.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/CaloClusters.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Qvectors.h"

#ifndef PWGEM_DILEPTON_DATAMODEL_DILEPTONTABLES_H_
#define PWGEM_DILEPTON_DATAMODEL_DILEPTONTABLES_H_

namespace o2::aod
{

namespace emevent
{
DECLARE_SOA_COLUMN(CollisionId, collisionId, int);
DECLARE_SOA_COLUMN(NeeULS, neeuls, int);
DECLARE_SOA_COLUMN(NeeLSpp, neelspp, int);
DECLARE_SOA_COLUMN(NeeLSmm, neelsmm, int);
DECLARE_SOA_COLUMN(Bz, bz, float);           //! kG
DECLARE_SOA_COLUMN(Q2xFT0M, q2xft0m, float); //! Qx for 2nd harmonics in FT0M
DECLARE_SOA_COLUMN(Q2yFT0M, q2yft0m, float); //! Qy for 2nd harmonics in FT0M
DECLARE_SOA_COLUMN(Q2xFT0A, q2xft0a, float); //! Qx for 2nd harmonics in FT0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q2yFT0A, q2yft0a, float); //! Qy for 2nd harmonics in FT0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q2xFT0C, q2xft0c, float); //! Qx for 2nd harmonics in FT0C (i.e. negative eta)
DECLARE_SOA_COLUMN(Q2yFT0C, q2yft0c, float); //! Qy for 2nd harmonics in FT0C (i.e. negative eta)
DECLARE_SOA_COLUMN(Q2xTPCpos, q2xbpos, float); //! Qx for 2nd harmonics in Barrel positive eta region
DECLARE_SOA_COLUMN(Q2yTPCpos, q2ybpos, float); //! Qy for 2nd harmonics in Barrel positive eta region
DECLARE_SOA_COLUMN(Q2xTPCneg, q2xbneg, float); //! Qx for 2nd harmonics in Barrel negative eta region
DECLARE_SOA_COLUMN(Q2yTPCneg, q2ybneg, float); //! Qy for 2nd harmonics in Barrel negative eta region
DECLARE_SOA_COLUMN(Q2xTPCall, q2xbtot, float); //! Qx for 2nd harmonics in Barrel full eta region
DECLARE_SOA_COLUMN(Q2yTPCall, q2ybtot, float); //! Qy for 2nd harmonics in Barrel full eta region
DECLARE_SOA_COLUMN(Q3xFT0M, q3xft0m, float); //! Qx for 3rd harmonics in FT0M
DECLARE_SOA_COLUMN(Q3yFT0M, q3yft0m, float); //! Qy for 3rd harmonics in FT0M
DECLARE_SOA_COLUMN(Q3xFT0A, q3xft0a, float); //! Qx for 3rd harmonics in FT0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q3yFT0A, q3yft0a, float); //! Qy for 3rd harmonics in FT0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q3xFT0C, q3xft0c, float); //! Qx for 3rd harmonics in FT0C (i.e. negative eta)
DECLARE_SOA_COLUMN(Q3yFT0C, q3yft0c, float); //! Qy for 3rd harmonics in FT0C (i.e. negative eta)
DECLARE_SOA_COLUMN(Q3xTPCpos, q3xbpos, float); //! Qx for 3rd harmonics in Barrel positive eta region
DECLARE_SOA_COLUMN(Q3yTPCpos, q3ybpos, float); //! Qy for 3rd harmonics in Barrel positive eta region
DECLARE_SOA_COLUMN(Q3xTPCneg, q3xbneg, float); //! Qx for 3rd harmonics in Barrel negative eta region
DECLARE_SOA_COLUMN(Q3yTPCneg, q3ybneg, float); //! Qy for 3rd harmonics in Barrel negative eta region
DECLARE_SOA_COLUMN(Q3xTPCall, q3xbtot, float); //! Qx for 3rd harmonics in Barrel full eta region
DECLARE_SOA_COLUMN(Q3yTPCall, q3ybtot, float); //! Qy for 3rd harmonics in Barrel full eta region
DECLARE_SOA_COLUMN(Q4xFT0M, q4xft0m, float); //! Qx for 4th harmonics in FT0M
DECLARE_SOA_COLUMN(Q4yFT0M, q4yft0m, float); //! Qy for 4th harmonics in FT0M
DECLARE_SOA_COLUMN(Q4xFT0A, q4xft0a, float); //! Qx for 4th harmonics in FT0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q4yFT0A, q4yft0a, float); //! Qy for 4th harmonics in FT0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q4xFT0C, q4xft0c, float); //! Qx for 4th harmonics in FT0C (i.e. negative eta)
DECLARE_SOA_COLUMN(Q4yFT0C, q4yft0c, float); //! Qy for 4th harmonics in FT0C (i.e. negative eta)
DECLARE_SOA_COLUMN(Q4xTPCpos, q4xbpos, float); //! Qx for 4th harmonics in Barrel positive eta region
DECLARE_SOA_COLUMN(Q4yTPCpos, q4ybpos, float); //! Qy for 4th harmonics in Barrel positive eta region
DECLARE_SOA_COLUMN(Q4xTPCneg, q4xbneg, float); //! Qx for 4th harmonics in Barrel negative eta region
DECLARE_SOA_COLUMN(Q4yTPCneg, q4ybneg, float); //! Qy for 4th harmonics in Barrel negative eta region
DECLARE_SOA_COLUMN(Q4xTPCall, q4xbtot, float); //! Qx for 4th harmonics in Barrel full eta region
DECLARE_SOA_COLUMN(Q4yTPCall, q4ybtot, float); //! Qy for 4th harmonics in Barrel full eta region

DECLARE_SOA_DYNAMIC_COLUMN(Sel8, sel8, [](uint64_t selection_bit) -> bool { return (selection_bit & BIT(o2::aod::evsel::kIsTriggerTVX)) && (selection_bit & BIT(o2::aod::evsel::kNoTimeFrameBorder)) && (selection_bit & BIT(o2::aod::evsel::kNoITSROFrameBorder)); });
DECLARE_SOA_DYNAMIC_COLUMN(EP2FT0M, ep2ft0m, [](float q2x, float q2y) -> float { return std::atan2(q2y, q2x) / 2.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP2FT0A, ep2ft0a, [](float q2x, float q2y) -> float { return std::atan2(q2y, q2x) / 2.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP2FT0C, ep2ft0c, [](float q2x, float q2y) -> float { return std::atan2(q2y, q2x) / 2.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP2TPCpos, ep2bpos, [](float q2x, float q2y) -> float { return std::atan2(q2y, q2x) / 2.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP2TPCneg, ep2bneg, [](float q2x, float q2y) -> float { return std::atan2(q2y, q2x) / 2.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP2TPCall, ep2btot, [](float q2x, float q2y) -> float { return std::atan2(q2y, q2x) / 2.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP3FT0M, ep3ft0m, [](float q3x, float q3y) -> float { return std::atan2(q3y, q3x) / 3.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP3FT0A, ep3ft0a, [](float q3x, float q3y) -> float { return std::atan2(q3y, q3x) / 3.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP3FT0C, ep3ft0c, [](float q3x, float q3y) -> float { return std::atan2(q3y, q3x) / 3.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP3TPCpos, ep3bpos, [](float q3x, float q3y) -> float { return std::atan2(q3y, q3x) / 3.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP3TPCneg, ep3bneg, [](float q3x, float q3y) -> float { return std::atan2(q3y, q3x) / 3.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP3TPCall, ep3btot, [](float q3x, float q3y) -> float { return std::atan2(q3y, q3x) / 3.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP4FT0M, ep4ft0m, [](float q4x, float q4y) -> float { return std::atan2(q4y, q4x) / 4.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP4FT0A, ep4ft0a, [](float q4x, float q4y) -> float { return std::atan2(q4y, q4x) / 4.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP4FT0C, ep4ft0c, [](float q4x, float q4y) -> float { return std::atan2(q4y, q4x) / 4.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP4TPCpos, ep4bpos, [](float q4x, float q4y) -> float { return std::atan2(q4y, q4x) / 4.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP4TPCneg, ep4bneg, [](float q4x, float q4y) -> float { return std::atan2(q4y, q4x) / 4.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP4TPCall, ep4btot, [](float q4x, float q4y) -> float { return std::atan2(q4y, q4x) / 4.0; });
} // namespace emevent

DECLARE_SOA_TABLE(EMEvents, "AOD", "EMEVENT", //!   Main event information table
                  o2::soa::Index<>, emevent::CollisionId, bc::RunNumber, bc::GlobalBC, evsel::Alias, evsel::Selection, timestamp::Timestamp,
                  collision::PosX, collision::PosY, collision::PosZ,
                  collision::NumContrib, evsel::NumTracksInTimeRange, emevent::Sel8<evsel::Selection>);
using EMEvent = EMEvents::iterator;

DECLARE_SOA_TABLE(EMEventsCov, "AOD", "EMEVENTCOV", //! joinable to EMEvents
                  collision::CovXX, collision::CovXY, collision::CovXZ, collision::CovYY, collision::CovYZ, collision::CovZZ, collision::Chi2);
using EMEventCov = EMEventsCov::iterator;

DECLARE_SOA_TABLE(EMEventsBz, "AOD", "EMEVENTBZ", emevent::Bz); // joinable to EMEvents
using EMEventBz = EMEventsBz::iterator;

DECLARE_SOA_TABLE(EMEventsMult, "AOD", "EMEVENTMULT", //!   event multiplicity table, joinable to EMEvents
                  mult::MultFT0A, mult::MultFT0C,
                  mult::MultTPC, mult::MultNTracksPV, mult::MultNTracksPVeta1, mult::MultNTracksPVetaHalf,
                  mult::IsInelGt0<mult::MultNTracksPVeta1>, mult::IsInelGt1<mult::MultNTracksPVeta1>, mult::MultFT0M<mult::MultFT0A, mult::MultFT0C>);
using EMEventMult = EMEventsMult::iterator;

DECLARE_SOA_TABLE(EMEventsCent, "AOD", "EMEVENTCENT", //!   event centrality table, joinable to EMEvents
                  cent::CentFT0M, cent::CentFT0A, cent::CentFT0C, cent::CentNTPV);
using EMEventCent = EMEventsCent::iterator;

DECLARE_SOA_TABLE(EMEventsQvec, "AOD", "EMEVENTQVEC", //!   event q vector table, joinable to EMEvents
                  emevent::Q2xFT0M, emevent::Q2yFT0M, emevent::Q2xFT0A, emevent::Q2yFT0A, emevent::Q2xFT0C, emevent::Q2yFT0C,
                  emevent::Q2xTPCpos, emevent::Q2yTPCpos, emevent::Q2xTPCneg, emevent::Q2yTPCneg, emevent::Q2xTPCall, emevent::Q2yTPCall,
                  emevent::Q3xFT0M, emevent::Q3yFT0M, emevent::Q3xFT0A, emevent::Q3yFT0A, emevent::Q3xFT0C, emevent::Q3yFT0C,
                  emevent::Q3xTPCpos, emevent::Q3yTPCpos, emevent::Q3xTPCneg, emevent::Q3yTPCneg, emevent::Q3xTPCall, emevent::Q3yTPCall,

                  // Dynamic columns
                  emevent::EP2FT0M<emevent::Q2xFT0M, emevent::Q2yFT0M>,
                  emevent::EP2FT0A<emevent::Q2xFT0A, emevent::Q2yFT0A>,
                  emevent::EP2FT0C<emevent::Q2xFT0C, emevent::Q2yFT0C>,
                  emevent::EP2TPCpos<emevent::Q2xTPCpos, emevent::Q2yTPCpos>,
                  emevent::EP2TPCneg<emevent::Q2xTPCneg, emevent::Q2yTPCneg>,
                  emevent::EP2TPCall<emevent::Q2xTPCall, emevent::Q2yTPCall>,
                  emevent::EP3FT0M<emevent::Q3xFT0M, emevent::Q3yFT0M>,
                  emevent::EP3FT0A<emevent::Q3xFT0A, emevent::Q3yFT0A>,
                  emevent::EP3FT0C<emevent::Q3xFT0C, emevent::Q3yFT0C>,
                  emevent::EP3TPCpos<emevent::Q3xTPCpos, emevent::Q3yTPCpos>,
                  emevent::EP3TPCneg<emevent::Q3xTPCneg, emevent::Q3yTPCneg>,
                  emevent::EP3TPCall<emevent::Q3xTPCall, emevent::Q3yTPCall>);
using EMEventQvec = EMEventsQvec::iterator;

DECLARE_SOA_TABLE(EMEventsNee, "AOD", "EMEVENTNEE", emevent::NeeULS, emevent::NeeLSpp, emevent::NeeLSmm); // joinable to EMEvents
using EMEventNee = EMEventsNee::iterator;

namespace emmcevent
{
DECLARE_SOA_COLUMN(McCollisionId, mcCollisionId, int);
} // namespace emmcevent

DECLARE_SOA_TABLE(EMMCEvents, "AOD", "EMMCEVENT", //!   MC event information table
                  o2::soa::Index<>, emmcevent::McCollisionId, mccollision::GeneratorsID,
                  mccollision::PosX, mccollision::PosY, mccollision::PosZ,
                  mccollision::T, mccollision::ImpactParameter,

                  // dynamic column
                  mccollision::GetGeneratorId<mccollision::GeneratorsID>,
                  mccollision::GetSubGeneratorId<mccollision::GeneratorsID>,
                  mccollision::GetSourceId<mccollision::GeneratorsID>);

using EMMCEvent = EMMCEvents::iterator;

namespace emmceventlabel
{
DECLARE_SOA_INDEX_COLUMN(EMMCEvent, emmcevent); //! MC collision
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);   //! Bit mask to indicate collision mismatches (bit ON means mismatch). Bit 15: indicates negative label
} // namespace emmceventlabel

DECLARE_SOA_TABLE(EMMCEventLabels, "AOD", "EMMCEVENTLABEL", //! Table joined to the EMEvents table containing the MC index
                  emmceventlabel::EMMCEventId, emmceventlabel::McMask);
using EMMCEventLabel = EMMCEventLabels::iterator;

namespace emmcparticle
{
DECLARE_SOA_INDEX_COLUMN(EMMCEvent, emmcevent);
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(Mothers, mothers);     //! Mother tracks (possible empty) array. Iterate over mcParticle.mothers_as<aod::McParticles>())
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(Daughters, daughters); //! Daughter tracks (possibly empty) array. Check for non-zero with mcParticle.has_daughters(). Iterate over mcParticle.daughters_as<aod::McParticles>())
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float px, float py) -> float { return RecoDecay::sqrtSumOfSquares(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](float px, float py, float pz) -> float { return RecoDecay::eta(std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](float px, float py) -> float { return RecoDecay::phi(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float px, float py, float pz) -> float { return RecoDecay::sqrtSumOfSquares(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(Y, y, //! Particle rapidity
                           [](float pz, float e) -> float {
                             if ((e - pz) > static_cast<float>(1e-7)) {
                               return 0.5f * std::log((e + pz) / (e - pz));
                             } else {
                               return -999.0f;
                             }
                           });
} // namespace emmcparticle

// This table contains all MC truth tracks (both v0 and calos)
DECLARE_SOA_TABLE_FULL(EMMCParticles, "EMMCParticles", "AOD", "EMMCPARTICLE", //!  MC track information (on disk)
                       o2::soa::Index<>, emmcparticle::EMMCEventId,
                       mcparticle::PdgCode, mcparticle::Flags,
                       emmcparticle::MothersIds, emmcparticle::DaughtersIds,
                       mcparticle::Px, mcparticle::Py, mcparticle::Pz, mcparticle::E,
                       mcparticle::Vx, mcparticle::Vy, mcparticle::Vz,

                       // dynamic column
                       emmcparticle::Pt<mcparticle::Px, mcparticle::Py>,
                       emmcparticle::Eta<mcparticle::Px, mcparticle::Py, mcparticle::Pz>,
                       emmcparticle::Phi<mcparticle::Px, mcparticle::Py>,

                       emmcparticle::P<mcparticle::Px, mcparticle::Py, mcparticle::Pz>,
                       emmcparticle::Y<mcparticle::Pz, mcparticle::E>,
                       mcparticle::ProducedByGenerator<mcparticle::Flags>,
                       mcparticle::FromBackgroundEvent<mcparticle::Flags>,
                       mcparticle::IsPhysicalPrimary<mcparticle::Flags>);

using EMMCParticle = EMMCParticles::iterator;

namespace emprimaryelectronmclabel
{
DECLARE_SOA_INDEX_COLUMN(EMMCParticle, emmcparticle); //!
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);
} // namespace emprimaryelectronmclabel

// NOTE: MC labels. This table has one entry for each reconstructed track (joinable with EMPrimaryElectrons table)
DECLARE_SOA_TABLE(EMPrimaryElectronMCLabels, "AOD", "EMPRMELMCLABEL", //!
                  emprimaryelectronmclabel::EMMCParticleId, emprimaryelectronmclabel::McMask);
using EMPrimaryElectronMCLabel = EMPrimaryElectronMCLabels::iterator;

namespace emprimarymuonmclabel
{
DECLARE_SOA_INDEX_COLUMN(EMMCParticle, emmcparticle); //!
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);
} // namespace emprimarymuonmclabel

// NOTE: MC labels. This table has one entry for each reconstructed track (joinable with EMPrimaryMuons table)
DECLARE_SOA_TABLE(EMPrimaryMuonMCLabels, "AOD", "EMPRMMUMCLABEL", //!
                  emprimarymuonmclabel::EMMCParticleId, emprimarymuonmclabel::McMask);
using EMPrimaryMuonMCLabel = EMPrimaryMuonMCLabels::iterator;

namespace emprimaryelectron
{
DECLARE_SOA_INDEX_COLUMN(EMEvent, emevent);        //!
DECLARE_SOA_COLUMN(CollisionId, collisionId, int); //!
DECLARE_SOA_COLUMN(TrackId, trackId, int);         //!
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(AmbiguousElectrons, ambiguousElectrons);
DECLARE_SOA_COLUMN(IsAssociatedToMPC, isAssociatedToMPC, bool); //! is associated to most probable collision
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                         //!
DECLARE_SOA_COLUMN(PrefilterBit, pfb, uint8_t);                 //!
DECLARE_SOA_DYNAMIC_COLUMN(Signed1Pt, signed1Pt, [](float pt, int8_t sign) -> float { return sign * 1. / pt; });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float pt, float eta) -> float { return pt * std::cosh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, [](float pt, float phi) -> float { return pt * std::cos(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, [](float pt, float phi) -> float { return pt * std::sin(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, [](float pt, float eta) -> float { return pt * std::sinh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(MeanClusterSizeITS, meanClusterSizeITS, [](uint32_t itsClusterSizes) -> float {
  int total_cluster_size = 0, nl = 0;
  for (unsigned int layer = 0; layer < 7; layer++) {
    int cluster_size_per_layer = (itsClusterSizes >> (layer * 4)) & 0xf;
    if (cluster_size_per_layer > 0) {
      nl++;
    }
    total_cluster_size += cluster_size_per_layer;
  }
  if (nl > 0) {
    return static_cast<float>(total_cluster_size) / static_cast<float>(nl);
  } else {
    return 0;
  }
});
DECLARE_SOA_DYNAMIC_COLUMN(MeanClusterSizeITSib, meanClusterSizeITSib, [](uint32_t itsClusterSizes) -> float {
  int total_cluster_size = 0, nl = 0;
  for (unsigned int layer = 0; layer < 3; layer++) {
    int cluster_size_per_layer = (itsClusterSizes >> (layer * 4)) & 0xf;
    if (cluster_size_per_layer > 0) {
      nl++;
    }
    total_cluster_size += cluster_size_per_layer;
  }
  if (nl > 0) {
    return static_cast<float>(total_cluster_size) / static_cast<float>(nl);
  } else {
    return 0;
  }
});
DECLARE_SOA_DYNAMIC_COLUMN(MeanClusterSizeITSob, meanClusterSizeITSob, [](uint32_t itsClusterSizes) -> float {
  int total_cluster_size = 0, nl = 0;
  for (unsigned int layer = 3; layer < 7; layer++) {
    int cluster_size_per_layer = (itsClusterSizes >> (layer * 4)) & 0xf;
    if (cluster_size_per_layer > 0) {
      nl++;
    }
    total_cluster_size += cluster_size_per_layer;
  }
  if (nl > 0) {
    return static_cast<float>(total_cluster_size) / static_cast<float>(nl);
  } else {
    return 0;
  }
});
} // namespace emprimaryelectron
DECLARE_SOA_TABLE(EMPrimaryElectrons, "AOD", "EMPRIMARYEL", //!
                  o2::soa::Index<>, emprimaryelectron::CollisionId,
                  emprimaryelectron::TrackId, emprimaryelectron::Sign,
                  track::Pt, track::Eta, track::Phi, track::DcaXY, track::DcaZ,
                  track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows,
                  track::TPCChi2NCl, track::TPCInnerParam,
                  track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaMu, pidtpc::TPCNSigmaPi, pidtpc::TPCNSigmaKa, pidtpc::TPCNSigmaPr,
                  pidtofbeta::Beta, pidtof::TOFNSigmaEl, pidtof::TOFNSigmaMu, pidtof::TOFNSigmaPi, pidtof::TOFNSigmaKa, pidtof::TOFNSigmaPr,
                  track::ITSClusterSizes, track::ITSChi2NCl, track::DetectorMap,
                  track::X, track::Alpha, track::Y, track::Z, track::Snp, track::Tgl, emprimaryelectron::IsAssociatedToMPC,

                  // dynamic column
                  track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::v001::ITSClusterMap<track::ITSClusterSizes>, track::v001::ITSNCls<track::ITSClusterSizes>, track::v001::ITSNClsInnerBarrel<track::ITSClusterSizes>,
                  track::HasITS<track::DetectorMap>, track::HasTPC<track::DetectorMap>,
                  track::HasTRD<track::DetectorMap>, track::HasTOF<track::DetectorMap>,
                  emprimaryelectron::Signed1Pt<track::Pt, emprimaryelectron::Sign>,
                  emprimaryelectron::P<track::Pt, track::Eta>,
                  emprimaryelectron::Px<track::Pt, track::Phi>,
                  emprimaryelectron::Py<track::Pt, track::Phi>,
                  emprimaryelectron::Pz<track::Pt, track::Eta>,
                  emprimaryelectron::MeanClusterSizeITS<track::ITSClusterSizes>,
                  emprimaryelectron::MeanClusterSizeITSib<track::ITSClusterSizes>,
                  emprimaryelectron::MeanClusterSizeITSob<track::ITSClusterSizes>);
// iterators
using EMPrimaryElectron = EMPrimaryElectrons::iterator;

DECLARE_SOA_TABLE(EMPrimaryElectronsCov, "AOD", "EMPRIMARYELCOV", //!
                  aod::track::CYY,
                  aod::track::CZY,
                  aod::track::CZZ,
                  aod::track::CSnpY,
                  aod::track::CSnpZ,
                  aod::track::CSnpSnp,
                  aod::track::CTglY,
                  aod::track::CTglZ,
                  aod::track::CTglSnp,
                  aod::track::CTglTgl,
                  aod::track::C1PtY,
                  aod::track::C1PtZ,
                  aod::track::C1PtSnp,
                  aod::track::C1PtTgl,
                  aod::track::C1Pt21Pt2);
// iterators
using EMPrimaryElectronCov = EMPrimaryElectronsCov::iterator;

DECLARE_SOA_TABLE(EMPrimaryElectronEMEventIds, "AOD", "PRMELMEVENTID", emprimaryelectron::EMEventId); // To be joined with EMPrimaryElectrons table at analysis level.
// iterators
using EMPrimaryElectronEMEventId = EMPrimaryElectronEMEventIds::iterator;

DECLARE_SOA_TABLE(EMPrimaryElectronsPrefilterBit, "AOD", "PRMELPFB", emprimaryelectron::PrefilterBit); // To be joined with EMPrimaryElectrons table at analysis level.
// iterators
using EMPrimaryElectronPrefilterBit = EMPrimaryElectronsPrefilterBit::iterator;

DECLARE_SOA_TABLE(EMAmbiguousElectronSelfIds, "AOD", "EMAMBELSELFID", emprimaryelectron::AmbiguousElectronsIds); // To be joined with EMPrimaryElectrons table at analysis level.
// iterators
using EMAmbiguousElectronSelfId = EMAmbiguousElectronSelfIds::iterator;

namespace emprimarymuon
{
DECLARE_SOA_INDEX_COLUMN(EMEvent, emevent);                                                     //!
DECLARE_SOA_COLUMN(CollisionId, collisionId, int);                                              //!
DECLARE_SOA_COLUMN(FwdTrackId, fwdtrackId, int);                                                //!
DECLARE_SOA_SELF_INDEX_COLUMN_FULL(MCHTrack, matchMCHTrack, int, "EMPRIMARYMUs_MatchMCHTrack"); //! Index of matched MCH track for GlobalMuonTracks and GlobalForwardTracks
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(AmbiguousMuons, ambiguousMuons);
DECLARE_SOA_COLUMN(IsAssociatedToMPC, isAssociatedToMPC, bool); //! is associated to most probable collision
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                         //!
DECLARE_SOA_COLUMN(Chi2MFTsa, chi2MFTsa, float);                //! chi2 of MFT standalone track
DECLARE_SOA_DYNAMIC_COLUMN(Signed1Pt, signed1Pt, [](float pt, int8_t sign) -> float { return sign * 1. / pt; });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float pt, float eta) -> float { return pt * std::cosh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, [](float pt, float phi) -> float { return pt * std::cos(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, [](float pt, float phi) -> float { return pt * std::sin(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, [](float pt, float eta) -> float { return pt * std::sinh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(DcaXY, dcaXY, [](float dcaX, float dcaY) -> float { return std::sqrt(dcaX * dcaX + dcaY * dcaY); });
DECLARE_SOA_DYNAMIC_COLUMN(NClustersMFT, nClustersMFT, //! Number of MFT clusters
                           [](uint64_t mftClusterSizesAndTrackFlags) -> uint8_t {
                             uint8_t nClusters = 0;
                             for (int layer = 0; layer < 10; layer++) {
                               if ((mftClusterSizesAndTrackFlags >> (layer * 6)) & 0x3F) {
                                 nClusters++;
                               }
                             }
                             return nClusters;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(MFTClusterMap, mftClusterMap, //! MFT cluster map, one bit per a layer, starting from the innermost
                           [](uint64_t mftClusterSizesAndTrackFlags) -> uint16_t {
                             uint16_t clmap = 0;
                             for (unsigned int layer = 0; layer < 10; layer++) {
                               if ((mftClusterSizesAndTrackFlags >> (layer * 6)) & 0x3f) {
                                 clmap |= (1 << layer);
                               }
                             }
                             return clmap;
                           });
} // namespace emprimarymuon
DECLARE_SOA_TABLE(EMPrimaryMuons, "AOD", "EMPRIMARYMU", //!
                  o2::soa::Index<>, emprimarymuon::CollisionId,
                  emprimarymuon::FwdTrackId, fwdtrack::TrackType,
                  fwdtrack::Pt, fwdtrack::Eta, fwdtrack::Phi, emprimarymuon::Sign,
                  fwdtrack::FwdDcaX, fwdtrack::FwdDcaY,
                  fwdtrack::X, fwdtrack::Y, fwdtrack::Z, fwdtrack::Tgl,

                  fwdtrack::NClusters, fwdtrack::PDca, fwdtrack::RAtAbsorberEnd,
                  fwdtrack::Chi2, fwdtrack::Chi2MatchMCHMID, fwdtrack::Chi2MatchMCHMFT,
                  // fwdtrack::MatchScoreMCHMFT, fwdtrack::MFTTrackId, fwdtrack::MCHTrackId,
                  emprimarymuon::MCHTrackId,
                  fwdtrack::MCHBitMap, fwdtrack::MIDBitMap, fwdtrack::MIDBoards,
                  fwdtrack::MFTClusterSizesAndTrackFlags, emprimarymuon::Chi2MFTsa, emprimarymuon::IsAssociatedToMPC,

                  // dynamic column
                  emprimarymuon::Signed1Pt<fwdtrack::Pt, emprimarymuon::Sign>,
                  emprimarymuon::NClustersMFT<fwdtrack::MFTClusterSizesAndTrackFlags>,
                  emprimarymuon::MFTClusterMap<fwdtrack::MFTClusterSizesAndTrackFlags>,
                  emprimarymuon::P<fwdtrack::Pt, fwdtrack::Eta>,
                  emprimarymuon::Px<fwdtrack::Pt, fwdtrack::Phi>,
                  emprimarymuon::Py<fwdtrack::Pt, fwdtrack::Phi>,
                  emprimarymuon::Pz<fwdtrack::Pt, fwdtrack::Eta>,
                  emprimarymuon::DcaXY<fwdtrack::FwdDcaX, fwdtrack::FwdDcaY>);
// iterators
using EMPrimaryMuon = EMPrimaryMuons::iterator;

DECLARE_SOA_TABLE(EMPrimaryMuonsCov, "AOD", "EMPRIMARYMUCOV", //!
                  aod::fwdtrack::CXX,
                  aod::fwdtrack::CXY,
                  aod::fwdtrack::CYY,
                  aod::fwdtrack::CPhiX,
                  aod::fwdtrack::CPhiY,
                  aod::fwdtrack::CPhiPhi,
                  aod::fwdtrack::CTglX,
                  aod::fwdtrack::CTglY,
                  aod::fwdtrack::CTglPhi,
                  aod::fwdtrack::CTglTgl,
                  aod::fwdtrack::C1PtX,
                  aod::fwdtrack::C1PtY,
                  aod::fwdtrack::C1PtPhi,
                  aod::fwdtrack::C1PtTgl,
                  aod::fwdtrack::C1Pt21Pt2);
// iterators
using EMPrimaryMuonCov = EMPrimaryMuonsCov::iterator;

DECLARE_SOA_TABLE(EMPrimaryMuonEMEventIds, "AOD", "PRMMUEMEVENTID", emprimarymuon::EMEventId); // To be joined with EMPrimaryMuons table at analysis level.
// iterators
using EMPrimaryMuonEMEventId = EMPrimaryMuonEMEventIds::iterator;

DECLARE_SOA_TABLE(EMAmbiguousMuonSelfIds, "AOD", "EMAMBMUSELFID", emprimarymuon::AmbiguousMuonsIds); // To be joined with EMPrimaryMuons table at analysis level.
// iterators
using EMAmbiguousMuonSelfId = EMAmbiguousMuonSelfIds::iterator;

} // namespace o2::aod

#endif // PWGEM_DILEPTON_DATAMODEL_DILEPTONTABLES_H_
