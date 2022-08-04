#include "utils.h"
#include "struct.h"

#include "TG4Event.h"
#include "TG4HitSegment.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TGeoManager.h>
#include <TGeoTrd2.h>
#include <TGeoTube.h>
#include <TObjString.h>
#include <TRandom3.h>

//#include "transf.h"
#include <iostream>

namespace fast_reco {
bool flukatype = false;

double stt_center[3];

TPRegexp *rST;
TPRegexp *r2ST;
TPRegexp *rSTplane;
TPRegexp *rSTmod;

std::map<int, std::map<double, int>> stX;
std::map<int, double> stL;
std::map<int, std::map<int, TVector2>> stPos;
std::map<int, TVector2> tubePos;
std::map<int, double> t0;
} // namespace fast_reco

void fast_reco::decodeModDPlaneID(int modDplaneID, int &moduleid,
                                  int &DPlaneid) {

  moduleid = modDplaneID / 10;
  DPlaneid = modDplaneID % 10;
}

void fast_reco::getSTinfo(TGeoNode *nod, TGeoHMatrix mat, int modDPlaneid,
                          std::map<double, int> &stX,
                          std::map<int, double> &stL,
                          std::map<int, TVector2> &stPos) {
                          
  int modid, dPlaneid;
  decodeModDPlaneID(modDPlaneid, modid, dPlaneid);
  int ic = (dPlaneid == 1) ? 0 : 1;

  if (ic != 0 && ic != 1)
    std::cout << "Error: ic expected 0 or 1 -> " << ic << std::endl;
  //std::cout<<" modDPlaneid:"<<modDPlaneid<<"  ndaughter:"<<nod->GetNdaughters()<<std::endl;
  for (int i = 0; i < nod->GetNdaughters(); i++) {
    auto n2straw = nod->GetDaughter(i);
    auto obja = fast_reco::r2ST->MatchS(n2straw->GetName());

    int n2straw_id =
        (reinterpret_cast<TObjString *>(obja->At(5)))->GetString().Atoi();
    delete obja;

    TGeoMatrix *n2strawmat = n2straw->GetMatrix();
    TGeoHMatrix n2strawhmat = mat * (*n2strawmat);

    for (int j = 0; j < n2straw->GetNdaughters(); j++) {
      TGeoNode *dau = n2straw->GetDaughter(j);
      TGeoTube *tub = (TGeoTube *)dau->GetVolume()->GetShape();
      double lenght = 2 * tub->GetDz();
      TString name = dau->GetName();

      if (!isST(name))
        std::cout << "Error: expected ST but not -> " << name.Data()
                  << std::endl;

      TGeoMatrix *thismat = n2straw->GetDaughter(j)->GetMatrix();
      TGeoHMatrix mymat = n2strawhmat * (*thismat);

      auto obja2 = fast_reco::rST->MatchS(name);

      int tid =
          (reinterpret_cast<TObjString *>(obja2->At(obja2->GetEntries() - 2)))
              ->GetString()
              .Atoi();
      delete obja2;

      int id = n2straw_id * 2 + tid;

      TVector2 v;
      v.SetX(mymat.GetTranslation()[2]);
      v.SetY(mymat.GetTranslation()[ic]);
      //// X is actually z here, Y() is y for horizontal straw, x for vertical straw
      stX[v.Y()] = id;

      /////////////
      /////////// the following is wrong stL with id as index
      //stL[id] = lenght;
      stL[modDPlaneid] = lenght;
      stPos[id] = v;
      //std::cout<<"  modDPlaneid:"<<modDPlaneid<<"   id:"<<id<<"  v.X:"<<v.X()<<" v.y:"<<v.Y()<<std::endl;
      // std::cout<<" fill mstX, id:"<<id<<"  v.Y():"<<v.Y()<<"
      // v.X():"<<v.X()<<" length:"<<lenght<<std::endl;
    }
  }
}

int fast_reco::getModDPlaneID(TString path) {
  //      TGeoNode *nod = gGeoManager->GetCurrentNode();
  //   TString path = gGeoManager->GetPath();
  //   TString name = nod->GetName();
  // /volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/MagIntVol_volume_PV_0/sand_inner_volume_PV_0/STTtracker_PV_0/STT_00_TrkMod_PV_0/STT_00_TrkMod_hh_PV_1
  // /volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/MagIntVol_volume_PV_0/sand_inner_volume_PV_0/STTtracker_PV_0/STT_26_C3H6Mod_PV_0/STT_26_C3H6Mod_ST_PV_0/STT_26_C3H6Mod_ST_vv_PV_0
  //   const char *const rSTplane_string =
  //       "STT_([0-9]+)_(Trk|C3H6|C)Mod(_ST|)_(hh|vv)_PV_([0-9]+)(/|)";
  //   const char *const rSTmod_string =
  //       "STT_([0-9]+)_(Trk|C3H6|C)Mod_PV_([0-9]+)(/|)";

  // const char *const strawCopy_inDL_hh =
  // "horizontalST_(Xe|Ar)_PV_([0-9]+)(/|)"; const char *const strawCopy_inDL_vv
  // = "vv_ST_PV_([0-9]+)(/|)";
  // const char *const strawCopy_singleLayer = "2straw_PV_([0-9]+)(/|)";

  // TPRegexp *rSTplane = new TPRegexp(rSTplane_string);
  // TPRegexp *rSTmod = new TPRegexp(rSTmod_string);
  // TPRegexp *rstrawCopy_hh = new TPRegexp(strawCopy_inDL_hh);
  // TPRegexp *rstrawCopy_vv = new TPRegexp(strawCopy_inDL_vv);
  // TPRegexp *rstrawCopy_singleLayer = new TPRegexp(strawCopy_singleLayer);

  TObjArray *obja = fast_reco::rSTplane->MatchS(path);
  TObjArray *obja2 = fast_reco::rSTmod->MatchS(path);
  // std::cout << "path:" << path << std::endl;
  // TObjArray *obj3 = rstrawCopy_hh->MatchS(path);
  // TObjArray *obj4 = rstrawCopy_vv->MatchS(path);
  // TObjArray *obj5 = rstrawCopy_singleLayer->MatchS(path);

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
  // int strawcopy_in_doubleLayer =
  //     hor ? (reinterpret_cast<TObjString *>(obj3->At(2)))->GetString().Atoi()
  //         : (reinterpret_cast<TObjString
  //         *>(obj4->At(1)))->GetString().Atoi();

  // std::cout<<" mod:"<<mod<<" modcopy:"<<icopymod<<"  dl:"<<doubleLayer<<"
  // icopy:"<<icopy<<" iindl:"<<strawcopy_in_doubleLayer<<"
  // id:"<<strawcopy_in_singleLayer<<std::endl;
  int modID = icopymod == 0 ? mod : (60 - mod);
  // int planeID = doubleLayer * 2 + strawcopy_in_doubleLayer;
  // int strawID=modID * 10000 + planeID * 1000 + strawcopy_in_singleLayer;

  return modID * 10 + doubleLayer;
}

// get position of the center of the tube for each plane
void fast_reco::getSTPlaneinfo(TGeoHMatrix mat,
                               std::map<int, std::map<double, int>> &stX,
                               std::map<int, double> &stL,
                               std::map<int, std::map<int, TVector2>> &stPos) {
  TGeoNode *nod = gGeoManager->GetCurrentNode();
  TString path = gGeoManager->GetPath();
  TString name = nod->GetName();
  TGeoMatrix *thismat = nod->GetMatrix();
  TGeoHMatrix mymat = mat * (*thismat);

  int modDplaneID = 0;
  double x = 0;
  double z = 0;
  // std::cout<<" path:"<<path<<std::endl;
  if (isSTPlane(name)) {
    modDplaneID = getModDPlaneID(path);
    //std::cout<<path<<" -> isSTPlane, modDplaneID:"<<modDplaneID<<std::endl;
    std::map<double, int> mstX;
    std::map<int, TVector2> mstPos;

    getSTinfo(nod, mymat, modDplaneID, mstX, stL, mstPos);

    stX[modDplaneID] = mstX;
    stPos[modDplaneID] = mstPos;

  } else {

    for (int i = 0; i < nod->GetNdaughters(); i++) {
      // std::cout<<" cd down i:"<<i<<std::endl;
      gGeoManager->CdDown(i);
      getSTPlaneinfo(mymat, stX, stL, stPos);
      gGeoManager->CdUp();
    }
  }
}

void fast_reco::init(TGeoManager *geo) {
  std::cout << " --------------------------- geo init "
               "-----------------------------------"
            << std::endl;
  //TGeoTube *ec = (TGeoTube *)geo->FindVolumeFast("STTtracker_PV")->GetShape();
  geo->cd("/volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/MagIntVol_volume_PV_0/sand_inner_volume_PV_0");

  TGeoHMatrix mat = *(gGeoManager->GetCurrentMatrix()); //*gGeoIdentity;
  
  TGeoNode *nod = gGeoManager->GetCurrentNode();
  for (int i = 0; i < nod->GetNdaughters(); i++) {
    TString name=nod->GetDaughter(i)->GetName();
    if(name.Contains("STTtracker_PV")){
      //std::cout<<" cd down to i:"<<i<<std::endl;
      gGeoManager->CdDown(i);
    }
  }
  double local[3]={0.,0.,0.};
  geo->LocalToMaster(local, stt_center);
  std::cout<<"--------------------------- stt_center, x: "<<stt_center[0]<<" y:"<<stt_center[1]<<" z:"<<stt_center[2]<<std::endl;

  fast_reco::rST = new TPRegexp(rST_string);
  fast_reco::r2ST = new TPRegexp(r2ST_string);
  fast_reco::rSTplane = new TPRegexp(rSTplane_string);
  fast_reco::rSTmod = new TPRegexp(rSTmod_string);

  getSTPlaneinfo(mat, fast_reco::stX, fast_reco::stL, fast_reco::stPos);

  for (std::map<int, std::map<int, TVector2>>::iterator it =
           fast_reco::stPos.begin();
       it != fast_reco::stPos.end(); it++) {
    int modDPlaneid = it->first;
    for (std::map<int, TVector2>::iterator ite = it->second.begin();
         ite != it->second.end(); ite++) {
      int dp_strawid = ite->first;
      int i_in_dplane = dp_strawid % 2;
      int i_in_straw = dp_strawid / 2;
      int DPlane = modDPlaneid % 10;
      int modID = modDPlaneid / 10;
      int planeID = DPlane * 2 + i_in_dplane;
      int id = modID * 10000 + planeID * 1000 + i_in_straw;
      tubePos[id] = ite->second;
    }
  }
  std::cout<<" ----- tubePos size:"<<tubePos.size()<<std::endl;
  std::cout << " --------------------------- geo init "
               "end-----------------------------------"
            << std::endl;
}

void fast_reco::initT0(TG4Event* ev)
{
  TRandom3 r(0);
  fast_reco::t0.clear();

  double t0_beam = ev->Primaries[0].Position.T() -
                   ev->Primaries[0].Position.Z() / fast_reco::c; 
                   // +
                   //r.Gaus(0, fast_reco::bucket_rms);
  //std::cout<<"   t0 size:"<<t0.size()<<"  t0_beam:"<<t0_beam<<"  ev->Primaries[0].Position.T():"<<ev->Primaries[0].Position.T()<<" z:"<<ev->Primaries[0].Position.Z()<<std::endl;

  for (std::map<int, std::map<int, TVector2> >::iterator it = stPos.begin();
       it != stPos.end(); it++) {
    fast_reco::t0[it->first] =
        t0_beam + it->second.begin()->second.X() / fast_reco::c;
        //std::cout<<it->first<<"  t:"<<fast_reco::t0[it->first]<<"  z:"<<it->second.begin()->second.X()<<std::endl;
  }
}


bool fast_reco::isST(TString name) { return name.Contains(*fast_reco::rST); }

bool fast_reco::isSTPlane(TString name) {
  return name.Contains(*fast_reco::rSTplane);
}