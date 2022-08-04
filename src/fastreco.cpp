#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TGeoTrd2.h>
#include <TGeoTube.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TRandom3.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

#include "TG4Event.h"
#include "TG4HitSegment.h"

#include "TPRegexp.h"
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

#include "struct.h"
#include "utils.h"

using namespace fast_reco;

TRandom3 r(0);

int getStrawID(double x, double y, double z) {
  TGeoNode *nod = gGeoManager->GetCurrentNode();
  TString path = gGeoManager->GetPath();
  TString name = nod->GetName();
  // /volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/MagIntVol_volume_PV_0/sand_inner_volume_PV_0/STTtracker_PV_0/STT_00_TrkMod_PV_0/STT_00_TrkMod_hh_PV_1
  // /volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/MagIntVol_volume_PV_0/sand_inner_volume_PV_0/STTtracker_PV_0/STT_26_C3H6Mod_PV_0/STT_26_C3H6Mod_ST_PV_0/STT_26_C3H6Mod_ST_vv_PV_0
  // const char *const rSTplane_string =
  //     "STT_([0-9]+)_(Trk|C3H6|C)Mod(_ST|)_(hh|vv)_PV_([0-9]+)(/|)";
  // const char *const rSTmod_string =
  //     "STT_([0-9]+)_(Trk|C3H6|C)Mod_PV_([0-9]+)(/|)";

  const char *const strawCopy_inDL_hh = "horizontalST_(Xe|Ar)_PV_([0-9]+)(/|)";
  const char *const strawCopy_inDL_vv = "vv_ST_PV_([0-9]+)(/|)";
  const char *const strawCopy_singleLayer = "2straw_PV_([0-9]+)(/|)";

  // TPRegexp *rSTplane = new TPRegexp(rSTplane_string);
  // TPRegexp *rSTmod = new TPRegexp(rSTmod_string);
  TPRegexp *rstrawCopy_hh = new TPRegexp(strawCopy_inDL_hh);
  TPRegexp *rstrawCopy_vv = new TPRegexp(strawCopy_inDL_vv);
  TPRegexp *rstrawCopy_singleLayer = new TPRegexp(strawCopy_singleLayer);

  TObjArray *obja = rSTplane->MatchS(path);
  TObjArray *obja2 = rSTmod->MatchS(path);
  // std::cout << "path:" << path << std::endl;
  TObjArray *obj3 = rstrawCopy_hh->MatchS(path);
  TObjArray *obj4 = rstrawCopy_vv->MatchS(path);
  TObjArray *obj5 = rstrawCopy_singleLayer->MatchS(path);

  // for (Int_t i = 0; i < obja->GetLast() + 1; i++) {
  //   const TString subStr = ((TObjString *)obja->At(i))->GetString();
  //   std::cout << "\"" << subStr << "\" ";
  // }
  // std::cout << std::endl;
  int mod = (reinterpret_cast<TObjString *>(obja->At(1)))->GetString().Atoi();
  int icopymod =
      (reinterpret_cast<TObjString *>(obja2->At(3)))->GetString().Atoi();
  //// icopy =1 only for trackingmod hh
  int icopy = (reinterpret_cast<TObjString *>(obja->At(5)))
                  ->GetString()
                  .Atoi(); //// icopy =1 only for trackingmod hh

  bool hor =
      (reinterpret_cast<TObjString *>(obja->At(4)))->GetString().EqualTo("hh");
  int doubleLayer = (hor ? 0 : 1) + 2 * icopy;
  int strawcopy_in_doubleLayer =
      hor ? (reinterpret_cast<TObjString *>(obj3->At(2)))->GetString().Atoi()
          : (reinterpret_cast<TObjString *>(obj4->At(1)))->GetString().Atoi();
  int strawcopy_in_singleLayer =
      (reinterpret_cast<TObjString *>(obj5->At(1)))->GetString().Atoi();

  // std::cout<<" mod:"<<mod<<" modcopy:"<<icopymod<<"  dl:"<<doubleLayer<<"
  // icopy:"<<icopy<<" iindl:"<<strawcopy_in_doubleLayer<<"
  // id:"<<strawcopy_in_singleLayer<<std::endl;
  int modID = icopymod == 0 ? mod : (60 - mod);
  int planeID = doubleLayer * 2 + strawcopy_in_doubleLayer;
  int strawID = modID * 10000 + planeID * 1000 + strawcopy_in_singleLayer;
  // std::cout<<" modID:"<<modID<<"  planeID:"<<planeID<<"
  // strawcopy_in_singleLayer:"<<strawcopy_in_singleLayer<<"
  // ---->:"<<strawID<<std::endl;
  return strawID;
}

void showAll(TG4Event *event, TGeoManager *geo) {
  std::cout << "=============================================" << std::endl;
  for (std::vector<TG4Trajectory>::iterator t = event->Trajectories.begin();
       t != event->Trajectories.end(); ++t) {
    std::cout << "   Traj " << t->TrackId;
    std::cout << " " << t->ParentId;
    std::cout << " " << t->Name;
    std::cout << " " << t->Points.size();
    std::cout << " E:" << t->GetInitialMomentum().E();
    std::cout << " beginpro:" << t->Points.begin()->Process << " "
              << t->Points.begin()->Subprocess
              << " endpro:" << (t->Points.end() - 1)->Process << " "
              << (t->Points.end() - 1)->Subprocess;
    std::cout << std::endl;
  }
  for (auto d = event->SegmentDetectors.begin();
       d != event->SegmentDetectors.end(); ++d) {
    std::cout << "   det " << d->first;
    std::cout << " " << d->second.size();
    // int count = 10;
    // std::cout << " up to " << count << " segments";
    std::cout << std::endl;
    if (d->first != "Straw")
      continue;

    int i = 0;
    int preStid = -1;
    for (std::vector<TG4HitSegment>::iterator h = d->second.begin();
         h != d->second.end(); ++h) {

      TLorentzVector mid = (h->Start + h->Stop) * 0.5;
      TString name = geo->FindNode(mid.X(), mid.Y(), mid.Z())->GetName();
      int stid = getStrawID(0, 0, 0);
      if (stid == preStid)
        std::cout << " !!! repeatTube:";
      std::cout << "   " << stid << " " << i;
      i++;
      std::cout << " P: " << h->PrimaryId << " " << h->Contrib[0];
      std::cout << " E: " << h->EnergyDeposit;
      std::cout << " S: " << h->SecondaryDeposit;
      std::cout << " C: " << h->Contrib.size() << "->";
      for (unsigned long j = 0; j < h->Contrib.size(); j++) {
        std::cout << " " << h->Contrib[j];
      }
      //      std::cout<<" name:"<<h->GetVolName();
      //            std::cout << " L: " << h->TrackLength;

      std::cout << " " << name;
      std::cout << " start:" << h->Start.X() << " " << h->Start.Y() << " "
                << h->Start.Z() << " " << h->Start.T()
                << " endT:" << h->Stop.T();
      if ((h + 1) != d->second.end() && (h + 1)->Start.T() < h->Start.T())
        std::cout << "   !!!!!!! time reverted";
      std::cout << std::endl;
      preStid = stid;
    }
  }

  std::cout << "=============================================" << std::endl;
}

// value of parameter of segment (y1,z1,y2,z2)
// corresponding to minimal distance to point (y,z)
double getT(double y1, double y2, double y, double z1, double z2, double z) {
  double t = 0;
  if (y1 != y2 || z1 != z2) {
    t = -((y1 - y) * (y2 - y1) + (z1 - z) * (z2 - z1)) /
        ((y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));
  }
  if (t < 0)
    return 0;
  else if (t > 1)
    return 1;
  else
    return t;
}

// Group hits into tube
void CollectHits(TG4Event *ev, TGeoManager *geo, int NHits,
                 Int_t DetType[10000], Float_t xPos[10000], Float_t yPos[10000],
                 Float_t zPos[10000],
                 std::map<int, std::vector<hit>> &hits2Tube) {
  hits2Tube.clear();

  for (unsigned int j = 0; j < ev->SegmentDetectors["Straw"].size(); j++) {
    const TG4HitSegment &hseg = ev->SegmentDetectors["Straw"].at(j);

    double x = 0.5 * (hseg.Start.X() + hseg.Stop.X());
    double y = 0.5 * (hseg.Start.Y() + hseg.Stop.Y());
    double z = 0.5 * (hseg.Start.Z() + hseg.Stop.Z());

    std::string sttname = "NULL";
    int stid = -999; // should be implemented for FLUKA

    sttname = geo->FindNode(x, y, z)->GetName();

    stid = getStrawID(x, y, z);
    // std::cout << " stid:" << stid;
    if (stid == -999) {
      std::cout << " stid is 999, something is wrong !!!!!!!!!!!!!!! check "
                << std::endl;
      continue;
    }

    hit h;
    h.det = sttname;
    h.did = stid;
    h.x1 = hseg.Start.X();
    h.y1 = hseg.Start.Y();
    h.z1 = hseg.Start.Z();
    h.t1 = hseg.Start.T();
    h.x2 = hseg.Stop.X();
    h.y2 = hseg.Stop.Y();
    h.z2 = hseg.Stop.Z();
    h.t2 = hseg.Stop.T();
    h.de = hseg.EnergyDeposit;
    h.pid = hseg.PrimaryId;
    h.index = j;

    hits2Tube[stid].push_back(h);
  }
}

// for each tube simulate tdc and adc
// tdc is the time of closest point to wire + drift time
// adc is the sum of energy deposit within integration time window
void Hits2Digit(std::map<int, std::vector<hit>> &hits2Tube,
                std::map<int, dg_tube> &digits) {
  digits.clear();

  for (std::map<int, std::vector<hit>>::iterator it = hits2Tube.begin();
       it != hits2Tube.end(); ++it) {
    double min_time_tub = 1E9; // mm
    int did = it->first;

    // int mod, tub, type, pla, plloc;
    double dwire = 0.;

    // decodeSTID(did, pla, tub);
    // decodePlaneID(pla, mod, plloc, type);
    int modID = did / 10000;
    int planeID = did % 10000 / 1000;
    int modDPlaneID = modID * 10 + planeID / 2;

    TVector2 wire = tubePos[did];

    dg_tube d;
    d.det = it->second[0].det;
    d.did = did;
    d.de = 0;
    d.hor = (modDPlaneID % 2 == 0);
    d.t0 = fast_reco::t0[modDPlaneID];

    if (d.hor == true) {
      d.x = fast_reco::stt_center[0];
      d.y = wire.Y();
      d.z = wire.X();
      // dwire = d.x - 0.5 * fast_reco::stL[did];

    } else {
      d.x = wire.Y();
      d.y = fast_reco::stt_center[1];
      d.z = wire.X();
      // dwire = d.y - 0.5 * fast_reco::stL[did];
    }

    for (unsigned int i = 0; i < it->second.size(); i++) {
      double x1 = it->second[i].z1;
      double x2 = it->second[i].z2;
      double t1 = it->second[i].t1;
      double t2 = it->second[i].t2;

      double y1, y2;

      if (d.hor) {
        y1 = it->second[i].y1;
        y2 = it->second[i].y2;
        dwire = (it->second[i].x1 + it->second[i].x2) / 2. - stt_center[0] +
                stL[modDPlaneID] / 2.;
      } else {
        y1 = it->second[i].x1;
        y2 = it->second[i].x2;
        dwire = stt_center[1] - (it->second[i].y1 + it->second[i].y2) / 2. +
                stL[modDPlaneID] / 2.;
      }

      ///// no need to use getT() 
      // double l = getT(y1, y2, wire.Y(), x1, x2, wire.X());
      // double x = x1 + (x2 - x1) * l;
      // double y = y1 + (y2 - y1) * l;
      // double t = t1 + (t2 - t1) * l;
      double x=(x1+x2)/2.;
      double y=(y1+y2)/2.;
      double t=(t1+t2)/2.;
      //std::cout<<" x:"<<x<<" y:"<<y<<" t:"<<t<<"  xm:"<<(x1+x2)/2.<<" ym:"<<(y1+y2)/2.<<"  tm:"<<(t1+t2)/2.<<std::endl;
      TVector2 min_dist_point(x, y);
      double min_dist_hit = (min_dist_point - wire).Mod();
      double min_time_hit =
          t + (min_dist_hit - fast_reco::wire_radius) / fast_reco::v_drift +
          dwire / fast_reco::v_signal_inwire;

      if (min_time_hit < min_time_tub)
        min_time_tub = min_time_hit;

      if (t - d.t0 < fast_reco::stt_int_time)
        d.de += it->second[i].de;

      d.hindex.push_back(it->second[i].index);

      std::cout<<"digitize: did:"<<did<<" t:"<<t<<"  min_dist_hit:"<<min_dist_hit<<" driftT:"<<(min_dist_hit - fast_reco::wire_radius) / fast_reco::v_drift<<" wirePropT:"<<dwire / fast_reco::v_signal_inwire<<" dwire:"<<dwire<<" wirelhalf:"<<stL[modDPlaneID] / 2.<<" tdc(nogaus):"<<min_time_hit<<std::endl;
    }

    d.tdc = min_time_tub; // + r.Gaus(0, fast_reco::tm_stt_smearing);
    d.adc = d.de;

    digits[did] = d;
  }
}

void digitizeSTT(TG4Event *event, TGeoManager *geo,
                 std::map<int, track> &track_map) {

  Int_t NHits;
  Float_t xHits[100000];
  Float_t yHits[100000];
  Float_t zHits[100000];
  Int_t DetType[100000];
  std::map<int, std::vector<hit>> hits2Tube;
  std::map<int, dg_tube> digits;

  CollectHits(event, geo, NHits, DetType, xHits, yHits, zHits, hits2Tube);
  Hits2Digit(hits2Tube, digits);

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////// NOTE: do not record discontinous hits
  //////// NOTE: do not record when time start to reverse
  //////// NOTE: it's possible same track cross same straw twice in connected
  /// two segment hits,avoid insert repeated tube into "track.clX/clY"

  int preTrackID = -1;
  int preStid = -1;
  double preT = -999.;
  std::unordered_set<int> done_tracks;
  std::unordered_set<int> done_tubes_for1track;

  for (unsigned int i = 0; i < event->SegmentDetectors["Straw"].size(); i++) {

    const TG4HitSegment &hseg = event->SegmentDetectors["Straw"].at(i);
    int id = hseg.PrimaryId;

    if (id != preTrackID) {
      done_tracks.insert(preTrackID);
      done_tubes_for1track.clear();
      preTrackID = id;
      preT=-999.;
    }
    if (done_tracks.find(id) == done_tracks.end()) {

      double x = 0.5 * (hseg.Start.X() + hseg.Stop.X());
      double y = 0.5 * (hseg.Start.Y() + hseg.Stop.Y());
      double z = 0.5 * (hseg.Start.Z() + hseg.Stop.Z());
      double t = 0.5 * (hseg.Start.T() + hseg.Stop.T());

      if (t < preT ) {
        std::cout << " time reverse for this track, dump the rest of hits for "
                     "trackid:"
                  << id << std::endl;
        done_tracks.insert(id);
      }

      geo->FindNode(x, y, z);
      int stid = getStrawID(x, y, z);
      if (stid == preStid) {
        std::cout << " !!!!! repeat:";
        geo->FindNode(hseg.Start.X(), hseg.Start.Y(), hseg.Start.Z());
        std::cout << " start:" << getStrawID(x, y, z);
        geo->FindNode(hseg.Stop.X(), hseg.Stop.Y(), hseg.Stop.Z());
        std::cout << "  stop:" << getStrawID(x, y, z);
      }
      std::cout << " stid:" << stid << " -> " << gGeoManager->GetPath()
                << std::endl;
      if (done_tubes_for1track.find(stid) == done_tubes_for1track.end()) {
        if (digits[stid].hor) {
          track_map[id].clY.push_back(std::move(digits[stid]));

        } else {
          track_map[id].clX.push_back(std::move(digits[stid]));
        }
      }
      done_tubes_for1track.insert(stid);
      preStid = stid;
      preT = t;
    }
  }
  auto cmp = [](dg_tube &a, dg_tube &b) { return a.did < b.did; };
  for (std::map<int, track>::iterator it = track_map.begin();
       it != track_map.end(); ++it) {

        it->second.tid=it->first;
    std::cout << "  track id:" << it->first
              << "  clX size:" << it->second.clX.size()
              << "  clY size:" << it->second.clY.size() << std::endl;
    std::sort(it->second.clX.begin(), it->second.clX.end(), cmp);
    std::sort(it->second.clY.begin(), it->second.clY.end(), cmp);
  }
}

double
get_traverse_from_neighborStraw(const std::vector<dg_tube *> &tubesP,
                                const std::vector<int> &neighbor_LayerList) {
  for (auto il : neighbor_LayerList) {
    for (auto tubep : tubesP) {
      int il2 = tubep->did % 10000 / 1000;
      if (il2 == il) {
        if (tubep->hor) {
          return tubep->y;
        } else {
          return tubep->x;
        }
      }
    }
  }
  return 99999.;
}

void get_driftDis_byModule(std::map<int, std::vector<dg_tube *>> &mod2tubes1,
                           std::map<int, std::vector<dg_tube *>> &mod2tubes2) {

  std::vector<std::vector<int>> neighborLayers_needList{
      {2, 3}, {2, 3}, {1, 0, 4, 5}, {4, 5, 1, 0}, {3, 2}, {3, 2}};
  for (std::map<int, std::vector<dg_tube *>>::iterator it = mod2tubes1.begin();
       it != mod2tubes1.end(); ++it) {
    int modid = it->first;
    for (unsigned int ix = 0; ix < it->second.size(); ix++) {
      dg_tube *tube = it->second[ix];
      int planeID = tube->did % 10000 / 1000;
      int modDPlaneID = modid * 10 + planeID / 2;
      double wireTravel = stL[modDPlaneID] / 2.;
      if (mod2tubes2.find(modid) != mod2tubes2.end()) {
        double dis = get_traverse_from_neighborStraw(
            mod2tubes2[modid], neighborLayers_needList[planeID]);
        if (dis != 99999.)
          wireTravel = tube->hor ? (wireTravel + dis - stt_center[0])
                                 : (wireTravel + stt_center[1] - dis);
        double driftT =
            tube->tdc - tube->t0 - wireTravel / fast_reco::v_signal_inwire;
        std::cout << " tdc:" << tube->tdc << "  t0:" << tube->t0
                  << "   wireTravel:" << wireTravel <<" wireLenhalf:"<<stL[modDPlaneID] / 2.<<" wirepropT:"<<wireTravel / fast_reco::v_signal_inwire<<"   driftT:" << driftT;
        tube->dis2center_reco = driftT * fast_reco::v_drift + fast_reco::wire_radius;
        std::cout << "  tube id:" << tube->did
                  << "  driftDis:" << tube->dis2center_reco << std::endl;
        if(driftT<0) std::cout<<" !!!!!!!!!!!!! driftT <0: "<<driftT<<std::endl;
      }
    }
  }
}
void cal_tube_drift_dis_forEachTrack(track &trk) {
  ////////////// only consider tubes within each track
  /////// layer id: 0-5
  ///  0, 1 : 2->3
  ///  4, 5 : 3->2
  ///  2:  1->0->4->5
  ///  3:  4->5 -> 1->0

  std::map<int, std::vector<dg_tube *>> mod2tubesX, mod2tubesY;

  for (unsigned int i = 0; i < trk.clX.size(); i++) {
    int did = trk.clX[i].did;
    int modid = did / 10000;
    mod2tubesX[modid].push_back(&trk.clX[i]);
  }
  for (unsigned int i = 0; i < trk.clY.size(); i++) {
    int did = trk.clY[i].did;
    int modid = did / 10000;
    mod2tubesY[modid].push_back(&trk.clY[i]);
  }
  std::cout << " ----- work on xstraws ----------" << std::endl;
  get_driftDis_byModule(mod2tubesX, mod2tubesY);
  std::cout << " ----- work on ystraws ----------" << std::endl;
  get_driftDis_byModule(mod2tubesY, mod2tubesX);
}

void reconstructSTT(TG4Event *event, TGeoManager *geo,
                    std::map<int, track> &track_map) {

  for (auto &trk : track_map) {
    std::cout << " ------------------------------  trackid:" << trk.first
              << std::endl;
    cal_tube_drift_dis_forEachTrack(trk.second);
  }
}

void processOneEntry(TG4Event *event, TGeoManager *geo) {
  std::map<int, track> track_map;
  digitizeSTT(event, geo, track_map);
  reconstructSTT(event, geo, track_map);
}

void FastReco(const char *finname, const char *foutname) {

  TFile f(finname, "READ");
  TTree *t = (TTree *)f.Get("EDepSimEvents");

  TG4Event *ev = new TG4Event;
  t->SetBranchAddress("Event", &ev);

  TGeoManager *geo = 0;

  geo = (TGeoManager *)f.Get("EDepSimGeometry");

  fast_reco::init(geo);

  TFile fout(foutname, "RECREATE");

  int startEntry = 1;
  int endEntry = 2; // t->GetEntries();
  for (int i = startEntry; i < endEntry; i++) {
    t->GetEntry(i);
    std::cout << " ----------------------------------- i:" << i
              << "----------------------------------" << std::endl;
    initT0(ev);
    showAll(ev, geo);
    processOneEntry(ev, geo);
  }
}

void help() {}

int main(int argc, char *argv[]) {
  if (argc != 3)
    help();
  else
    FastReco(argv[1], argv[2]);
}
