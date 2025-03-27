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

/// \file filterlambdalambda.cxx
/// \brief event selection to exclusively collect two reconstructed (anti-)lambdas
/// \author Junlee Kim, (junlee.kim@cern.ch)

#include <Framework/Configurable.h>
#include <Math/GenVector/Boost.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <TMath.h>
#include <fairlogger/Logger.h>

#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "EventFiltering/Zorro.h"
#include "EventFiltering/ZorroSummary.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "CommonConstants/MathConstants.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"

#include "../filterTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct filterlambdalambda {

  // Produce derived tables
  Produces<aod::LambdaLambdaFilters> tags;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};

  Configurable<bool> cfgUseGlobalTrack{"cfgUseGlobalTrack", true, "use Global track"};
  Configurable<float> cfgCutPt{"cfgCutPt", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 0.2f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 0.2f, "DCAz range for tracks"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 50, "Number of TPC cluster"};

  Configurable<float> cfgv0radiusMin{"cfgv0radiusMin", 1.2, "minimum decay radius"};
  Configurable<float> cfgDCAPosToPVMin{"cfgDCAPosToPVMin", 0.1, "minimum DCA to PV for positive track"};
  Configurable<float> cfgDCANegToPVMin{"cfgDCANegToPVMin", 0.1, "minimum DCA to PV for negative track"};
  Configurable<float> cfgv0CosPA{"cfgv0CosPA", 0.995, "minimum v0 cosine"};
  Configurable<float> cfgDCAV0Dau{"cfgDCAV0Dau", 1.0, "maximum DCA between daughters"};

  Configurable<float> cfgV0PtMin{"cfgV0PtMin", 0, "minimum pT for lambda"};
  Configurable<float> cfgV0EtaMin{"cfgV0EtaMin", -0.5, "maximum rapidity"};
  Configurable<float> cfgV0EtaMax{"cfgV0EtaMax", 0.5, "maximum rapidity"};
  Configurable<float> cfgV0LifeTime{"cfgV0LifeTime", 30., "maximum lambda lifetime"};

  Configurable<int> cfgDaughTPCnclsMin{"cfgDaughTPCnclsMin", 50, "minimum fired crossed rows"};
  Configurable<float> cfgDaughPIDCutsTPCPr{"cfgDaughPIDCutsTPCPr", 3, "proton nsigma for TPC"};
  Configurable<float> cfgDaughPIDCutsTPCPi{"cfgDaughPIDCutsTPCPi", 3, "pion nsigma for TPC"};
  Configurable<float> cfgDaughEtaMin{"cfgDaughEtaMin", -0.8, "minimum daughter eta"};
  Configurable<float> cfgDaughEtaMax{"cfgDaughEtaMax", 0.8, "maximum daughter eta"};
  Configurable<float> cfgDaughPrPt{"cfgDaughPrPt", 0.5, "minimum daughter proton pt"};
  Configurable<float> cfgDaughPiPt{"cfgDaughPiPt", 0.2, "minimum daughter pion pt"};

  Configurable<float> cfgHypMassWindow{"cfgHypMassWindow", 0.01, "Lambda mass selection window"};

  Configurable<float> cfgV0V0RapMax{"cfgV0V0RapMax", 0.5, "V0V0 rapidity selection"};
  Configurable<bool> cfgV0V0Sel{"cfgV0V0Sel", true, "application of V0V0 selections"};
  Configurable<float> cfgV0V0Radius{"cfgV0V0Radius", 10.0, "maximum radius of v0v0"};
  Configurable<float> cfgV0V0CPA{"cfgV0V0CPA", 0.96, "minimum CPA of v0v0"};
  Configurable<float> cfgV0V0Distance{"cfgV0V0Distance", 10, "minimum distance of v0v0"};
  Configurable<float> cfgV0V0DCA{"cfgV0V0DCA", 0.1, "maximum DCA of v0v0"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPt);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTPCFullPi, aod::pidTPCFullPr>>;

  OutputObj<TH1D> hProcessedEvents{TH1D("hProcessedEvents", ";; Number of events", 3, 0.0f, 3.0f)};

  double massLambda = o2::constants::physics::MassLambda;
  double massPr = o2::constants::physics::MassProton;
  double massPi = o2::constants::physics::MassPionCharged;

  void init(o2::framework::InitContext&)
  {
    hProcessedEvents->GetXaxis()->SetBinLabel(1, "All events");
    hProcessedEvents->GetXaxis()->SetBinLabel(2, "Events with F1");
    hProcessedEvents->GetXaxis()->SetBinLabel(3, aod::filtering::TriggerEventLambdaLambda::columnLabel());
  }

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (cfgUseGlobalTrack && !(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.tpcNClsFound() > cfgTPCcluster)) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool selectionPID(const T& candidate, int pid)
  {
    if (pid == 0) {
      if (std::abs(candidate.tpcNSigmaPi()) > cfgDaughPIDCutsTPCPi) {
        return false;
      }
    } else if (pid == 2) {
      if (std::abs(candidate.tpcNSigmaPr()) > cfgDaughPIDCutsTPCPr) {
        return false;
      }
    }
    return true;
  }

  template <typename V0>
  bool selectionV0(V0 const& candidate)
  {
    if (candidate.v0radius() < cfgv0radiusMin)
      return false;
    if (std::abs(candidate.dcapostopv()) < cfgDCAPosToPVMin)
      return false;
    if (std::abs(candidate.dcanegtopv()) < cfgDCANegToPVMin)
      return false;
    if (candidate.v0cosPA() < cfgv0CosPA)
      return false;
    if (std::abs(candidate.dcaV0daughters()) > cfgDCAV0Dau)
      return false;
    if (candidate.pt() < cfgV0PtMin)
      return false;
    if (candidate.yLambda() < cfgV0EtaMin)
      return false;
    if (candidate.yLambda() > cfgV0EtaMax)
      return false;

    return true;
  }

  template <typename T>
  bool selectionV0Daughter(T const& track, int pid) // pid 0: proton, pid 1: pion
  {
    if (track.tpcNClsFound() < cfgDaughTPCnclsMin)
      return false;
    if (pid == 0 && std::abs(track.tpcNSigmaPr()) > cfgDaughPIDCutsTPCPr)
      return false;
    if (pid == 1 && std::abs(track.tpcNSigmaPi()) > cfgDaughPIDCutsTPCPi)
      return false;
    if (track.eta() > cfgDaughEtaMax)
      return false;
    if (track.eta() < cfgDaughEtaMin)
      return false;
    if (pid == 0 && track.pt() < cfgDaughPrPt)
      return false;
    if (pid == 1 && track.pt() < cfgDaughPiPt)
      return false;

    return true;
  }

  template <typename V01, typename V02>
  float getDCAofV0V0(V01 const& v01, V02 const& v02)
  {
    ROOT::Math::XYZVector v01pos, v02pos, v01mom, v02mom;
    v01pos.SetXYZ(v01.x(), v01.y(), v01.z());
    v02pos.SetXYZ(v02.x(), v02.y(), v02.z());
    v01mom.SetXYZ(v01.px(), v01.py(), v01.pz());
    v02mom.SetXYZ(v02.px(), v02.py(), v02.pz());

    ROOT::Math::XYZVector posdiff = v02pos - v01pos;
    ROOT::Math::XYZVector cross = v01mom.Cross(v02mom);
    ROOT::Math::XYZVector dcaVec = (posdiff.Dot(cross) / cross.Mag2()) * cross;
    return std::sqrt(dcaVec.Mag2());
  }

  template <typename V01, typename V02>
  float getCPA(V01 const& v01, V02 const& v02)
  {
    ROOT::Math::XYZVector v01mom, v02mom;
    v01mom.SetXYZ(v01.px() / v01.p(), v01.py() / v01.p(), v01.pz() / v01.p());
    v02mom.SetXYZ(v02.px() / v02.p(), v02.py() / v02.p(), v02.pz() / v02.p());
    return v01mom.Dot(v02mom);
  }

  template <typename V01, typename V02>
  float getDistance(V01 const& v01, V02 const& v02)
  {
    ROOT::Math::XYZVector v01pos, v02pos;
    v01pos.SetXYZ(v01.x(), v01.y(), v01.z());
    v02pos.SetXYZ(v02.x(), v02.y(), v02.z());
    ROOT::Math::XYZVector posdiff = v02pos - v01pos;
    return std::sqrt(posdiff.Mag2());
  }

  template <typename V01, typename V02>
  float getRadius(V01 const& v01, V02 const& v02)
  {
    ROOT::Math::XYZVector v01pos, v02pos, v01mom, v02mom;
    v01pos.SetXYZ(v01.x(), v01.y(), v01.z());
    v02pos.SetXYZ(v02.x(), v02.y(), v02.z());
    v01mom.SetXYZ(v01.px() / v01.p(), v01.py() / v01.p(), v01.pz() / v01.p());
    v02mom.SetXYZ(v02.px() / v02.p(), v02.py() / v02.p(), v02.pz() / v02.p());
    ROOT::Math::XYZVector posdiff = v02pos - v01pos;

    float d = 1. - TMath::Power(v01mom.Dot(v02mom), 2);
    if (d < 1e-5)
      return 999;
    float t = posdiff.Dot(v01mom - v01mom.Dot(v02mom) * v02mom) / d;
    float s = -posdiff.Dot(v02mom - v01mom.Dot(v02mom) * v01mom) / d;
    ROOT::Math::XYZVector dca = v01pos + v02pos + t * v01mom + s * v02mom;
    dca /= 2.;
    return std::sqrt(dca.Mag2());
  }

  template <typename V01, typename V02>
  bool isSelectedV0V0(V01 const& v01, V02 const& v02)
  {
      return false;
    if (getDCAofV0V0(v01, v02) > cfgV0V0DCA)
      return false;
    if (getCPA(v01, v02) < cfgV0V0CPA)
      return false;
    if (getDistance(v01, v02) > cfgV0V0Distance)
      return false;
    if (getRadius(v01, v02) > cfgV0V0Radius)
      return false;

    return true;
  }

  ROOT::Math::PxPyPzMVector RecoV01, RecoV02, RecoV0V0;
  void processLLReducedTable(EventCandidates::iterator const& collision, TrackCandidates const& /*tracks*/, aod::V0Datas const& V0s, aod::BCsWithTimestamps const&)
  {
    bool keepEventDoubleLambda = false;
    hProcessedEvents->Fill(0.5);

    if (!collision.sel8())
      return;

    for (auto& v01 : V0s) {
      auto postrack_v01 = v01.template posTrack_as<TrackCandidates>();
      auto negtrack_v01 = v01.template negTrack_as<TrackCandidates>();

      int LambdaTag = 0;
      int aLambdaTag = 0;

      if (selectionV0Daughter(postrack_v01, 0) && selectionV0Daughter(negtrack_v01, 1)) {
        LambdaTag = 1;
      }
      if (selectionV0Daughter(negtrack_v01, 0) && selectionV0Daughter(postrack_v01, 1)) {
        aLambdaTag = 1;
      }

      if (LambdaTag == aLambdaTag)
        continue;

      if (!selectionV0(v01))
        continue;

      if (LambdaTag) {
        if (std::abs(massLambda - v01.mLambda()) > cfgHypMassWindow)
          continue;
        RecoV01 = ROOT::Math::PxPyPzMVector(v01.px(), v01.py(), v01.pz(), v01.mLambda());
      } else if (aLambdaTag) {
        if (std::abs(massLambda - v01.mAntiLambda()) > cfgHypMassWindow)
          continue;
        RecoV01 = ROOT::Math::PxPyPzMVector(v01.px(), v01.py(), v01.pz(), v01.mAntiLambda());
      }

      for (auto& v02 : V0s) {
        if (v01.v0Id() <= v02.v0Id())
          continue;
        auto postrack_v02 = v02.template posTrack_as<TrackCandidates>();
        auto negtrack_v02 = v02.template negTrack_as<TrackCandidates>();

        LambdaTag = 0;
        aLambdaTag = 0;

        if (selectionV0Daughter(postrack_v02, 0) && selectionV0Daughter(negtrack_v02, 1)) {
          LambdaTag = 1;
        }
        if (selectionV0Daughter(negtrack_v02, 0) && selectionV0Daughter(postrack_v02, 1)) {
          aLambdaTag = 1;
        }

        if (LambdaTag == aLambdaTag)
          continue;

        if (!selectionV0(v02))
          continue;

        if (postrack_v01.globalIndex() == postrack_v02.globalIndex() || postrack_v01.globalIndex() == negtrack_v02.globalIndex() || negtrack_v01.globalIndex() == postrack_v02.globalIndex() || negtrack_v01.globalIndex() == negtrack_v02.globalIndex())
          continue; // no shared decay products

        if (LambdaTag) {
          if (std::abs(massLambda - v02.mLambda()) > cfgHypMassWindow)
            continue;
          RecoV02 = ROOT::Math::PxPyPzMVector(v02.px(), v02.py(), v02.pz(), v02.mLambda());
        } else if (aLambdaTag) {
          if (std::abs(massLambda - v02.mAntiLambda()) > cfgHypMassWindow)
            continue;
          RecoV02 = ROOT::Math::PxPyPzMVector(v02.px(), v02.py(), v02.pz(), v02.mAntiLambda());
        }

        RecoV0V0 = RecoV01 + RecoV02;

        if (std::abs(RecoV0V0.Rapidity()) > cfgV0V0RapMax)
          continue;

        hProcessedEvents->Fill(1.5);

        if (isSelectedV0V0(v01, v02))
          keepEventDoubleLambda = true;
      }
    }
    if (keepEventDoubleLambda) {
      hProcessedEvents->Fill(2.5);
    }
  } // process
  PROCESS_SWITCH(filterlambdalambda, processLLReducedTable, "Process table creation for double ll", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<filterlambdalambda>(cfg, TaskName{"lf-lambdalambda-filter"})};
}

