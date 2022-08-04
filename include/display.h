#ifndef DISPLAY_H
#define DISPLAY_H

#include <TApplication.h>
#include <TArrow.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TFile.h>
#include <TGeoEltu.h>
#include <TGeoManager.h>
#include <TGeoTrd2.h>
#include <TGeoTube.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TLine.h>
#include <TMarker.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

#include <iomanip>
#include <iostream>
#include <map>
#include <stdlib.h>

#include "struct.h"

#include "TG4Event.h"
#include "TG4HitSegment.h"


namespace display
{
const bool debug = false;

static const int nMod = 24;
static const int nLay = 5;
static const int nCel = 12;
static const int nLay_ec = 5;
static const int nCel_ec = 90;

static const int nTotCells = nMod * nLay * nCel;
static const int nCellModule = nLay * nCel;

static const double dt = 500;


TCanvas* cev = 0;


int palette = 87;


}

void init(){

    gStyle->SetPalette(palette);

}

void show(){
     if (!initialized) {
    std::cout << "not initialized" << std::endl;
    return;
  }

  if (cev == 0) {
    cev = new TCanvas("cev", TString::Format("Event: %d", index).Data(), 1200,
                      600);
    cev->Divide(2, 1);
  } else {
    cev->SetTitle(TString::Format("Event: %d", index).Data());
  }
}

#endif