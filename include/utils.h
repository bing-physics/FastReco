#ifndef UTILS_H
#define UTILS_H
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TPRegexp.h>
#include <TString.h>
#include <TVector2.h>

#include "struct.h"

#include <TG4Event.h>


namespace fast_reco
{
const bool debug = false;

extern bool flukatype;  // for FLUKA

extern std::map<int, std::map<double, int> > stX;
extern std::map<int, double> stL;
extern std::map<int, std::map<int, TVector2> > stPos;
extern std::map<int, TVector2> tubePos;
extern std::map<int, double> t0;

extern double stt_center[3];

const double k = 0.299792458;
const double B = 0.6;
const double GeV_to_MeV = 1000.;
const double c = k * 1E3;  // mm/ns

const double bucket_rms = 1.;         // ns
const double tm_stt_smearing = 1;   // ns
const double wire_radius = 0.02;      // mm
const double v_drift = 0.05;          // mm/ns
const double v_signal_inwire = 200.;  // mm/ns
const double stt_int_time = 400.;     // ns




void getSTinfo(TGeoNode* nod, TGeoHMatrix mat, int planemodid,
               std::map<double, int>& stX, std::map<int, double>& stL,
               std::map<int, TVector2>& stPos);
void getSTPlaneinfo(TGeoHMatrix mat, std::map<int, std::map<double, int> >& stX,
                    std::map<int, double>& stL,
                    std::map<int, std::map<int, TVector2> >& stPos);

void decodeModDPlaneID(int modDPlaneid, int& moduleid, int& DPlaneid);

void init(TGeoManager* geo);
void initT0(TG4Event* ev);

const char* const rST_string =
    "(horizontalST_(Ar|Xe)|STT_([0-9]+)_(Trk|C3H6|C)Mod(_ST|)_vv_ST)_PV_([0-9]+"
    ")(/|)";
// "_(C3H6|C|Tr)Mod_([0-9]+)_(ST_|)(hor|ver|hor2)_ST_stGas_(Xe|Ar)19_vol_PV_(["
// "0-9]+)";
const char* const r2ST_string =
    "STT_([0-9]+)_(Trk|C3H6|C)Mod_(ST_|)(hh|vv)_2straw_PV_([0-9]+)(/|)";

const char* const rSTplane_string =
    "STT_([0-9]+)_(Trk|C3H6|C)Mod(_ST|)_(hh|vv)_PV_([0-9]+)(/|)";
// "_(C3H6|C|Tr)Mod_([0-9]+)_(ST_|)(hor|ver|hor2)_vol_PV_0";
const char* const rSTmod_string =
    "STT_([0-9]+)_(Trk|C3H6|C)Mod_PV_([0-9]+)(/|)";

extern TPRegexp* rST;
extern TPRegexp* rSTmod;
extern TPRegexp* r2ST;
extern TPRegexp* rSTplane;


bool isST(TString name);
bool isSTPlane(TString name);
int getModDPlaneID(TString path);

}




#endif
