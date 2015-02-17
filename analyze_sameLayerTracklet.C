#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TNtuple.h"
#include "Tracklet.h"
#include "Math/Vector3D.h"
#include "TMath.h"
#include "TF1.h"

using namespace std;

void analyze_sameLayerTracklet(char * infile, char * outfile = "outputSameLayer.root", int makeVzCut = 1,
                          int addL1Bck = 0, int addL2Bck = 0, double
			  splitProb = 0, double dropProb = 0, int nPileUp = 0,
			  double beamHaloRatio = 0.2,
			  bool putBeamHalo = false, char * beamHaloFile =
    			             "DataSample/PixelTree-Run123151-Full.root",
                          double smearVertex = 0.0,
			  bool putPixelTree = 0,
			  bool useKKVertex = 1,
			  bool useNSD = 0,
			  bool reWeight =1
			  )
{
  TFile* inf = new  TFile(infile);
  TTree* t = dynamic_cast<TTree*>(inf->FindObjectAny("PixelTree"));

  int zbins = 1;
  int hitbins = 100;
  int vertexHitRegion = 500000;
  int nbins = zbins*hitbins;
  double mult = 0;
  
  // Selection on Hits and events =====================================================================
  SelectionCriteria cuts;
  cuts.drCut   = 0.4;      // to remove double hit
  cuts.dPhiCut = 0.04;     // to remove double hit
  cuts.dEtaCut = 0.2;     // to remove double hit
  cuts.vzCut   = 10;       // vertex cut
  

  // Settings =========================================================================================
  cuts.verbose_ = false;
  cuts.useDeltaPhi_ = true;
  cuts.useDeltaRho_ = true;
  cuts.checkSecondLayer_ = true;
  
  // Output PDF =======================================================================================
  TFile* outf = new TFile(outfile,"recreate");
  TNtuple *ntmult = new TNtuple("ntmult","","mult:nhit1:nhit2");
  TNtuple *nthit = new TNtuple("nthit","","phi1:layer");
  TTree *trackletTree11 = new TTree("TrackletTree11","Tree of Reconstructed Tracklets");
  TTree *trackletTree22 = new TTree("TrackletTree22","Tree of Reconstructed Tracklets");
  TTree *trackletTree33 = new TTree("TrackletTree33","Tree of Reconstructed Tracklets");

  // Tracklet Tree data format ========================================================================
  TrackletData tdata11;
  TrackletData tdata22;
  TrackletData tdata33;

  setTrackletTreeBranch(trackletTree11,tdata11);
  setTrackletTreeBranch(trackletTree22,tdata22);
  setTrackletTreeBranch(trackletTree33,tdata33);



  // Hit vectors & pdfs ===================================================================================
  vector<TH1D*> layer1HitEta;
  layer1HitEta.reserve(nbins);

  vector<TH1D*> layer1HitPhi;
  layer1HitPhi.reserve(nbins);

  vector<TH2D*> layer12D;
  layer12D.reserve(nbins);

  vector<TH1D*> layer2HitEta;
  layer2HitEta.reserve(nbins);

  vector<TH1D*> layer2HitPhi;
  layer2HitPhi.reserve(nbins);

  vector<TH2D*> layer22D;
  layer22D.reserve(nbins);


  vector<TH1D*> layer3HitEta;
  layer3HitEta.reserve(nbins);

  vector<TH1D*> layer3HitPhi;
  layer3HitPhi.reserve(nbins);

  vector<TH2D*> layer32D;
  layer32D.reserve(nbins);
    
  for(int i = 0; i< nbins; ++i){
    layer1HitEta[i] = new TH1D(Form("dNdEtaHits1_%02d",i),"dNdEta Hits Layer 1",500,-3,3);
    layer2HitEta[i] = new TH1D(Form("dNdEtaHits2_%02d",i),"dNdEta Hits Layer 2",500,-3,3);
    layer3HitEta[i] = new TH1D(Form("dNdEtaHits3_%02d",i),"dNdEta Hits Layer 3",500,-3,3);
    layer1HitPhi[i] = new TH1D(Form("dNdPhiHits1_%02d",i),"dNdPhi Hits Layer 1",500,-3.2,3.2);
    layer2HitPhi[i] = new TH1D(Form("dNdPhiHits2_%02d",i),"dNdPhi Hits Layer 2",500,-3.2,3.2);
    layer3HitPhi[i] = new TH1D(Form("dNdPhiHits3_%02d",i),"dNdPhi Hits Layer 3",500,-3,3);
    layer12D[i] = new TH2D(Form("dNdEtadPhiHits1_%02d",i),"dNdPhidEta Hits Layer 1",500,-3,3,500,-3.2,3.2);
    layer22D[i] = new TH2D(Form("dNdEtadPhiHits2_%02d",i),"dNdPhidEta Hits Layer 2",500,-3,3,500,-3.2,3.2);
    layer32D[i] = new TH2D(Form("dNdEtadPhiHits3_%02d",i),"dNdPhidEta Hits Layer 3",500,-3,3,500,-3.2,3.2);
  }


  TH3F* nhits = new TH3F("nhits","",100,0,100,100,0,100,100,0,100);
  
  // Parameters for the tree =============================================================================
  Parameters par;  

  t->SetBranchAddress("eta1",par.eta1);
  t->SetBranchAddress("phi1",par.phi1);
  t->SetBranchAddress("ch1",par.ch1);
  t->SetBranchAddress("cs1",par.cs1);
  t->SetBranchAddress("r1",par.r1);
  t->SetBranchAddress("eta2",par.eta2);
  t->SetBranchAddress("phi2",par.phi2);
  t->SetBranchAddress("ch2",par.ch2);
  t->SetBranchAddress("cs2",par.cs2);
  t->SetBranchAddress("r2",par.r2);
  t->SetBranchAddress("eta3",par.eta3);
  t->SetBranchAddress("phi3",par.phi3);
  t->SetBranchAddress("ch3",par.ch3);
  t->SetBranchAddress("cs3",par.cs3);
  t->SetBranchAddress("r3",par.r3);
  t->SetBranchAddress("nhits1",&par.nhits1);
  t->SetBranchAddress("nhits2",&par.nhits2);
  t->SetBranchAddress("nhits3",&par.nhits3);
  t->SetBranchAddress("mult",&par.mult);
  t->SetBranchAddress("vz",par.vz);
  t->SetBranchAddress("nv",&par.nv);
  t->SetBranchAddress("npart",&par.npart);
  t->SetBranchAddress("eta",&par.eta);
  t->SetBranchAddress("phi",&par.phi);
  t->SetBranchAddress("pt",&par.pt);
  t->SetBranchAddress("chg",&par.chg);
  t->SetBranchAddress("pdg",&par.pdg);
  if(!makeVzCut) t->SetBranchAddress("evtType",&par.evtType);
  cout <<"Number of Events: "<<t->GetEntries()<<endl;

  int nBeamHalo = 0;

  // Main loop ===========================================================================================
  for(int i =0;  i<t->GetEntries()&&i<10000 ; i = i + 1 + nPileUp){    
    t->GetEntry(i);
    
    if (i % 1000 == 0) {
       cout <<"Event "<<i<<" "
            <<trackletTree11->GetEntries()<<" "
            <<trackletTree22->GetEntries()<<" "
            <<trackletTree33->GetEntries()<<" Add Beam Halo: "
	    <<nBeamHalo<<" "<<nBeamHalo/(double)i
            <<endl;
       if (reWeight) cout <<"Reweighted!!!!!!!"<<endl;    
    }       

    // Selection on Events
     bool reWeightDropFlag = 0;
    if (reWeight) {
       reWeightDropFlag = 0;
       double myVz = par.vz[1];
       if (myVz<-90) {
          TF1 *f = new TF1("f","gaus",-30,30);
          f->SetParameters(1,-2.709,4.551);
          myVz = f->GetRandom();
	  delete f;
       }
       double MCPdf = TMath::Gaus(myVz,-2.709,4.551,1);
       double DataPdf = TMath::Gaus(myVz,-2.702,3.627,1);
       double Ratio = DataPdf / MCPdf;
       //cout <<MCPdf<<" "<<DataPdf<<" "<<Ratio<<endl;
       double x=gRandom->Rndm()*2;

       if (x> Ratio) reWeightDropFlag=1;
       //cout <<x<<" "<<Ratio<<" "<<reWeightDropFlag<<endl;
    }
       
    // Fill vertex information
    tdata11.nv = par.nv;  
    tdata33.nv = par.nv;  
    tdata22.nv = par.nv;  
    for(int j = 0; j<par.nv;j++) {
       tdata11.vz[j] = par.vz[j];
       tdata33.vz[j] = par.vz[j];
       tdata22.vz[j] = par.vz[j];
    }
    
    // add background 
    if (addL1Bck!=0) {
       double etaMin = -2.5;
       double etaMax = 2.5;
       double phiMin = -3.1415926;
       double phiMax = 3.1415926;
       
       for (int i = par.nhits1; i < par.nhits1+addL1Bck; i++) {

          double eta = etaMin + (etaMax-etaMin)*gRandom->Rndm();
          double phi = phiMin + (phiMax-phiMin)*gRandom->Rndm();
	  double r = 4.15;
	  par.eta1[i] = eta;
	  par.phi1[i] = phi;
	  par.r1[i] = r;

       }
       par.nhits1+=addL1Bck;
    }

    if (addL2Bck!=0) {
       double etaMin = -2.5;
       double etaMax = 2.5;
       double phiMin = -3.1415926;
       double phiMax = 3.1415926;

       for (int i = par.nhits2; i < par.nhits2+addL2Bck; i++) {
          double eta = etaMin + (etaMax-etaMin)*gRandom->Rndm();
          double phi = phiMin + (phiMax-phiMin)*gRandom->Rndm();
	  double r = 7.05;
	  par.eta2[i] = eta;
	  par.phi2[i] = phi;
	  par.r2[i] = r;
       }
       par.nhits2+=addL2Bck;

    }

    
    // add trackletVertex
    if (tdata11.nv == 2) tdata11.nv=3;
    if (tdata33.nv == 2) tdata33.nv=3;
    if (tdata22.nv == 2) tdata22.nv=3;
    
    vector<RecoHit> layerRaw1;
    prepareHits(layerRaw1,par, cuts, 1, 0, splitProb, dropProb);
    vector<RecoHit> layerRaw2;
    prepareHits(layerRaw2,par, cuts, 2, 0, splitProb, dropProb);


    if (nPileUp!=0 && par.nhits1<10 ) {
       for (int j=1;j <= nPileUp ; j++) {
          t->GetEntry(i+j);
          prepareHits(layerRaw1,par, cuts, 1, 0, splitProb, dropProb);
          prepareHits(layerRaw2,par, cuts, 2, 0, splitProb, dropProb);
       }
       t->GetEntry(i);
    }

    
    double trackletVertex = 0;
    if (par.nhits1> vertexHitRegion) trackletVertex = par.vz[1]; else trackletVertex = TrackletVertexUnbin(layerRaw1,layerRaw2,0.14,0.08);
    tdata11.vz[tdata11.nv] = trackletVertex;
    tdata33.vz[tdata33.nv] = trackletVertex;
    tdata22.vz[tdata22.nv] = trackletVertex;
    tdata11.nv++;
    tdata33.nv++;
    tdata22.nv++;
    
    // use trackletVertex
    if (fabs(tdata11.vz[tdata11.nv-1])>cuts.vzCut && makeVzCut == 1) continue;

    // Process hits with Vz constraint:
    vector<RecoHit> layer1;
    prepareHits(layer1,par, cuts, 1, tdata11.vz[tdata11.nv-1], splitProb, dropProb);
    vector<RecoHit> layer2;
    prepareHits(layer2,par, cuts, 2, tdata11.vz[tdata11.nv-1], splitProb, dropProb);
    vector<RecoHit> layer3;
    prepareHits(layer3,par, cuts, 3, tdata11.vz[tdata11.nv-1], splitProb, dropProb);

    if (nPileUp!=0) {
       for (int j=1;j <= nPileUp ; j++) {
          t->GetEntry(i+j);
          prepareHits(layer1,par, cuts, 1, tdata11.vz[tdata11.nv-1], splitProb, dropProb);
          prepareHits(layer2,par, cuts, 2, tdata11.vz[tdata11.nv-1], splitProb, dropProb);
          prepareHits(layer3,par, cuts, 3, tdata11.vz[tdata11.nv-1], splitProb, dropProb);
       }
       t->GetEntry(i);
    }
    

    for(int ihit = 0; ihit< (int)layer1.size(); ++ihit) {
      int hitbin1 = (int)layer1.size();
      if (hitbin1 > 99) hitbin1 = 99;
      layer1HitEta[hitbin1]->Fill(layer1[ihit].eta);
      layer1HitPhi[hitbin1]->Fill(layer1[ihit].phi);
      layer12D[hitbin1]->Fill(layer1[ihit].eta,layer1[ihit].phi);
      nthit->Fill(layer1[ihit].phi,1);
      if(fabs(layer1[ihit].eta)<1) mult++;
    }

    for(int ihit = 0; ihit< (int)layer2.size(); ++ihit) {
      int hitbin2 = (int)layer2.size();
      if (hitbin2 > 99) hitbin2 = 99;
      layer2HitEta[hitbin2]->Fill(layer2[ihit].eta);
      layer2HitPhi[hitbin2]->Fill(layer2[ihit].phi);
      layer22D[hitbin2]->Fill(layer2[ihit].eta,layer2[ihit].phi);
      nthit->Fill(layer2[ihit].phi,2);
    }

    for(int ihit = 0; ihit< (int)layer3.size(); ++ihit) {
      int hitbin3 = (int)layer3.size();
      if (hitbin3 > 99) hitbin3 = 99;
      layer3HitEta[hitbin3]->Fill(layer3[ihit].eta);
      layer3HitPhi[hitbin3]->Fill(layer3[ihit].phi);
      layer32D[hitbin3]->Fill(layer3[ihit].eta,layer3[ihit].phi);
      nthit->Fill(layer3[ihit].phi,3);
    }

    // Form Tracklets        

    vector<Tracklet> recoTracklets11 = recoProtoTracklets(layer1,layer1);
    vector<Tracklet> recoTracklets22 = recoProtoTracklets(layer2,layer2);
    vector<Tracklet> recoTracklets33 = recoProtoTracklets(layer3,layer3);

/*
    vector<Tracklet> protoTracklets11 = recoProtoTracklets(layer1,layer1);
    vector<Tracklet> protoTracklets22 = recoProtoTracklets(layer2,layer2);
    vector<Tracklet> protoTracklets33 = recoProtoTracklets(layer3,layer3);
    
    vector<Tracklet> recoTracklets11 = cleanTracklets(protoTracklets11,0,cuts);
    vector<Tracklet> recoTracklets22 = cleanTracklets(protoTracklets22,0,cuts);
    vector<Tracklet> recoTracklets33 = cleanTracklets(protoTracklets33,0,cuts);
*/
    // Fill Ntuple
    tdata11.nTracklet = recoTracklets11.size();
    tdata11.nhit1 = layer1.size();
    tdata11.nhit2 = layer2.size();

    for (int j=0;j<(int)tdata11.nTracklet;j++)
    {
        tdata11.eta1[j] = recoTracklets11[j].eta1();	
        tdata11.eta2[j] = recoTracklets11[j].eta2();	
        tdata11.phi1[j] = recoTracklets11[j].phi1();	
        tdata11.phi2[j] = recoTracklets11[j].phi2();
        tdata11.r1[j]   = recoTracklets11[j].r1();	
        tdata11.r2[j]   = recoTracklets11[j].r2();
	tdata11.deta[j] = recoTracklets11[j].deta();
   	tdata11.dphi[j] = recoTracklets11[j].dphi();
	if (abs(tdata11.dphi[j])>3.2)cout <<tdata11.dphi[j]<<endl;
	tdata11.mult = (int)mult;	
    }

    tdata11.npart=0;
    for (int j=0;j<12;j++) tdata11.nhad[j]=0;

    for(int j=0;j<par.npart;j++)
    {
        if (fabs(par.eta[j])>3||par.chg[j]==0) continue;
	tdata11.eta[tdata11.npart]=par.eta[j];
	tdata11.chg[tdata11.npart]=par.chg[j];
	tdata11.pdg[tdata11.npart]=par.pdg[j];
	tdata11.pt[tdata11.npart]=par.pt[j];
        tdata11.npart++;
	int bin = (int)((par.eta[j]+3)*2);
	int pdg = (int)abs(par.pdg[j]);
	if (pdg== 211 || pdg == 321 || pdg == 2212 || pdg ==3122) tdata11.nhad[bin]++;
    }

    nhits->Fill(mult,layer1.size(),layer2.size());
    ntmult->Fill(mult,layer1.size(),layer2.size());

    tdata11.evtType = par.evtType;
    trackletTree11->Fill();

    tdata22.nTracklet = recoTracklets22.size();
    tdata22.nhit1 = layer1.size();
    tdata22.nhit2 = layer3.size();
    
    for (int j=0;j<(int)tdata22.nTracklet;j++)
    {
        tdata22.eta1[j] = recoTracklets22[j].eta1();	
        tdata22.eta2[j] = recoTracklets22[j].eta2();	
        tdata22.phi1[j] = recoTracklets22[j].phi1();	
        tdata22.phi2[j] = recoTracklets22[j].phi2();
        tdata22.r1[j]   = recoTracklets22[j].r1();	
        tdata22.r2[j]   = recoTracklets22[j].r2();
	tdata22.deta[j] = recoTracklets22[j].deta();
   	tdata22.dphi[j] = recoTracklets22[j].dphi();
	tdata22.mult = (int)mult;	
    }

    tdata22.npart=0;
    for (int j=0;j<12;j++) tdata22.nhad[j]=0;

    for(int j=0;j<par.npart;j++)
    {
        if (fabs(par.eta[j])>3||par.chg[j]==0) continue;
	tdata22.eta[tdata22.npart]=par.eta[j];
	tdata22.chg[tdata22.npart]=par.chg[j];
	tdata22.pdg[tdata22.npart]=par.pdg[j];
	tdata22.pt[tdata22.npart]=par.pt[j];
        tdata22.npart++;
	int bin = (int)((par.eta[j]+3)*2);
	int pdg = (int)abs(par.pdg[j]);
	if (pdg== 211 || pdg == 321 || pdg == 2213 || pdg ==3132) tdata22.nhad[bin]++;
    }

    nhits->Fill(mult,layer1.size(),layer2.size());
    ntmult->Fill(mult,layer1.size(),layer2.size());

    tdata22.evtType = par.evtType;
    trackletTree22->Fill();

    tdata33.nTracklet = recoTracklets33.size();
    tdata33.nhit1 = layer2.size();
    tdata33.nhit2 = layer3.size();

    for (int j=0;j<(int)tdata33.nTracklet;j++)
    {
        tdata33.eta1[j] = recoTracklets33[j].eta1();	
        tdata33.eta2[j] = recoTracklets33[j].eta2();	
        tdata33.phi1[j] = recoTracklets33[j].phi1();	
        tdata33.phi2[j] = recoTracklets33[j].phi2();
        tdata33.r1[j]   = recoTracklets33[j].r1();	
        tdata33.r2[j]   = recoTracklets33[j].r2();
	tdata33.deta[j] = recoTracklets33[j].deta();
   	tdata33.dphi[j] = recoTracklets33[j].dphi();
	tdata33.mult = (int)mult;	
    }

    tdata33.npart=0;
    for (int j=0;j<12;j++) tdata33.nhad[j]=0;

    for(int j=0;j<par.npart;j++)
    {
        if (fabs(par.eta[j])>3||par.chg[j]==0) continue;
	tdata33.eta[tdata33.npart]=par.eta[j];
	tdata33.chg[tdata33.npart]=par.chg[j];
	tdata33.pdg[tdata33.npart]=par.pdg[j];
	tdata33.pt[tdata33.npart]=par.pt[j];
        tdata33.npart++;
	int bin = (int)((par.eta[j]+3)*2);
	int pdg = (int)abs(par.pdg[j]);
	if (pdg== 211 || pdg == 321 || pdg == 2223 || pdg ==3232) tdata33.nhad[bin]++;
    }

    nhits->Fill(mult,layer1.size(),layer2.size());
    ntmult->Fill(mult,layer1.size(),layer2.size());

    tdata33.evtType = par.evtType;


    trackletTree33->Fill();

    
  }

  // Close outputfile ===================================================================================
  outf->Write();
  outf->Close(); 
}


