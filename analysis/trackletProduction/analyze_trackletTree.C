/*
#include <iostream>
#include <vector>
#include "TTree.h"
#include "TH1D.h"
#include "TH3F.h"
#include "Math/Vector3D.h"
#include "TMath.h"
*/

#include <TH2.h>
#include <TH3.h>
#include <TNtuple.h>
#include <TTimeStamp.h>
#include <TF1.h>
#include <TFile.h>

#include "pdfs.h"
#include "Tracklet.h"

using namespace std;

void analyze_trackletTree(char * infile, char * outfile = "output.root", long startEntry = 0, long endEntry = 1000000000,
                          int addL1Bck = 0, int addL2Bck = 0, int addL3Bck = 0,
   		          bool reWeight = 0,         // reweight to Run 123596 vtx distribution
                          bool useRandomVertex= 0,
			  bool cutOnClusterSize = 0,
			  bool mimicPixelCounting = 0, 
			  int makeVzCut = 0,
			  double splitProb = 0, double dropProb = 0, 
			  double nPileUp = 0,
			  double beamHaloRatio = 0.0,
			  bool putBeamHalo = false, char * beamHaloFile =
    			             "DataSample/PixelTree-Run123151-Full.root",
                          double smearVertex = 0,
			  bool putPixelTree = 0,
			  bool useKKVertex = 0,
			  bool useNSD = 0,
			  bool checkDoubleEvent = 0,
			  bool reweightMultiplicity = 0
			  )
{
  // Set Random Seed
  TTimeStamp myTime;
  gRandom->SetSeed(myTime.GetNanoSec());
  cout <<"Randomize "<<gRandom->Rndm()<<endl;
  
  // Input file =======================================================================================
  TFile* inf = new  TFile(infile);
  TTree* t = dynamic_cast<TTree*>(inf->FindObjectAny("PixelTree"));
  TFile* beamHaloInf;
  TTree* beamHaloTree;
  if (putBeamHalo) {
      cout <<"Add Beam Halo Background!!!!"<<endl;
      cout <<"Ratio = "<<beamHaloRatio<<endl;
      cout <<"File = "<<beamHaloFile<<endl;
      beamHaloInf = new TFile(beamHaloFile);
      beamHaloTree = dynamic_cast<TTree*>(beamHaloInf->FindObjectAny("PixelTree"));
  }
  
  // Output file =======================================================================================
  TFile* outf = new TFile(outfile,"recreate");
  TNtuple *ntmult = new TNtuple("ntmult","","mult:nhit1:nhit2");
  TNtuple *nthit = new TNtuple("nthit","","phi1:layer");
  TTree *trackletTree12 = new TTree("TrackletTree12","Tree of Reconstructed Tracklets");
  TTree *trackletTree13 = new TTree("TrackletTree13","Tree of Reconstructed Tracklets");
  TTree *trackletTree23 = new TTree("TrackletTree23","Tree of Reconstructed Tracklets");
  TTree *outTree;
  if (putPixelTree) {
     outTree = t->CloneTree();    
     cout <<"Put in Pixel Tree"<<endl;
  }
  int zbins = 1;
  int hitbins = 100;
  int vertexHitRegion = 500000;
  int nbins = zbins*hitbins;
  double mult = 0;
  bool isMC = 0;
  double vzShift = 0;
  
  TH1F *hRatio;
  if (reweightMultiplicity) {
     hRatio = getMultRatio();
  }
  
  // Selection on Hits and events =====================================================================
  SelectionCriteria cuts;
  cuts.drCut   = 0.4;      // to remove double hit
  cuts.dPhiCut = 0.04;     // to remove double hit
  cuts.dEtaCut = 0.2;     // to remove double hit
  cuts.vzCut   = 10;       // vertex cut
  

  // Settings =========================================================================================
  cuts.verbose_ = false;
  cuts.useDeltaPhi_ = false;
  cuts.useDeltaRho_ = false;
  cuts.checkSecondLayer_ = true;
  
  // Tracklet Tree data format ========================================================================
  TrackletData tdata12;
  TrackletData tdata13;
  TrackletData tdata23;

  setTrackletTreeBranch(trackletTree12,tdata12);
  setTrackletTreeBranch(trackletTree13,tdata13);
  setTrackletTreeBranch(trackletTree23,tdata23);

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

  TH3D* hLayer1Hit = new TH3D("hLayer1Hit","",75,0,15,60,-3,3,64,-3.2,3.2);
  TH3D* hLayer2Hit = new TH3D("hLayer2Hit","",75,0,15,60,-3,3,64,-3.2,3.2);
  TH3D* hLayer3Hit = new TH3D("hLayer3Hit","",75,0,15,60,-3,3,64,-3.2,3.2);

  // Prepare hit spectra for random hit 
  cout <<"Projecting...1"<<endl;
  if (addL1Bck!=0) t->Project("hLayer1Hit","phi1:eta1:r1");
  cout <<"Projecting...2"<<endl;
  if (addL2Bck!=0) t->Project("hLayer2Hit","phi2:eta2:r2");
  cout <<"Projecting...3"<<endl;
  if (addL3Bck!=0) t->Project("hLayer3Hit","phi3:eta3:r3");
  cout <<"Projecting...done"<<endl;

  if (t->FindBranch("npart")!=0) {
     
     isMC=true;
     cout <<"This is a Monte Carlo study."<<endl;
     vzShift = -0.4847;
     cout <<"vzShift = "<<vzShift<<endl;
  } else {
     cout <<"This is a data analysis."<<endl;
     smearVertex = 0;
  }
                                 
  // Event record
  vector <int> events[500];
  
  // Parameters for the tree =============================================================================
  Parameters par;  
  
  //  Parameters beamHaloPar;
  if (putBeamHalo) getPixelTreeBranch(beamHaloTree,par);  
  getPixelTreeBranch(t,par);
  
  if(!makeVzCut) t->SetBranchAddress("evtType",&par.evtType);
  cout <<"Number of Events: "<<t->GetEntries()<<endl;

  int nBeamHalo = 0;
  int doPileUp = 0;
  if (nPileUp!=0) {
     doPileUp=1;
     cout <<"Do pileup! With probability of "<<nPileUp<<endl;
  }   
  // Main loop ==========================================================================================
  for(int i = startEntry;  i<t->GetEntries()&&i<endEntry ; i = i + 1 + doPileUp){    
    t->GetEntry(i);
    if (i % 1000 == 0) {
       cout <<"Run "<<par.nRun<<" Event "<<i<<" "
            <<trackletTree12->GetEntries()<<" "
            <<trackletTree23->GetEntries()<<" "
            <<trackletTree13->GetEntries()<<" Add Beam Halo: "
	    <<nBeamHalo<<" "<<nBeamHalo/(double)i
            <<endl;
       if (reWeight) cout <<"Reweighted!!!!!!!"<<endl;    
    }       

    bool flagDoubleEvent = 0;

    if (checkDoubleEvent)
    {
       for (int j=0;j<events[par.nLumi].size();j++)
       {
          if (par.nEv==events[par.nLumi][j]) {
	     flagDoubleEvent = 1;
	     continue;
	  }
       }
       if (!flagDoubleEvent) events[par.nLumi].push_back(par.nEv);
    }

    if (flagDoubleEvent) continue;
//    if (par.nRun!=124023||(par.nRun==124033&&(par.nLumi<41||par.nLumi>96))) continue;
//    if (par.nRun!=132440||(par.nRun==132440&&(par.nLumi<141||par.nLumi>200))) continue;
    bool reWeightDropFlag = 0;
    if (reWeight) {
       reWeightDropFlag = 0;
       double myVz = par.vz[1];
       if (myVz<-90) {
          TF1 *f = new TF1("f","gaus",-30,30);
          f->SetParameters(1,-0.6536,4.438);
          myVz = f->GetRandom();
	  delete f;
       }
       
       // for early data 900 GeV
       /*
       double MCPdf = TMath::Gaus(myVz,-2.709,4.551,1);
       double DataPdf = TMath::Gaus(myVz,-2.702,3.627,1);
       */
       
       // for early data 7000 GeV Run 132440
       double MCPdf = TMath::Gaus(myVz,-0.6536,4.438,1);
       double DataPdf = TMath::Gaus(myVz,0.3533-vzShift,2.161,1);

//       double DataPdf = TMath::Gaus(myVz,-0.4623,2.731,1);
       double Ratio = DataPdf / MCPdf;
       //cout <<MCPdf<<" "<<DataPdf<<" "<<Ratio<<endl;
       double x=gRandom->Rndm()*2.5;

       if (x> Ratio) reWeightDropFlag=1;
       //cout <<x<<" "<<Ratio<<" "<<reWeightDropFlag<<endl;
    }
    
    if (reWeightDropFlag) continue;
    /*
    // Filter by evt selection cut //
    if (par.l1TBit[40]==0&&par.l1TBit[41]==0) continue;
    if (par.l1TBit[0]==0) continue;
    if (par.nLumi<69||par.nLumi>144) continue;
    */

    // Filter NSD events ==========================================================
    if ((par.evtType==92||par.evtType==93)&&useNSD) continue;
    
    // Filter HF coincidence
    if ((par.nHFn==0||par.nHFp==0||par.vz[1]<-99)&&reweightMultiplicity) continue;    
    
    // Beam Halo ==================================================================
    bool beamHaloFlag = false;
    
    if ( gRandom->Rndm() < beamHaloRatio && putBeamHalo ) {
       nBeamHalo++;
       beamHaloFlag = true;
       bool selectFlag = false;
       while (!selectFlag) {
          int nEntry = beamHaloTree->GetEntries();
          beamHaloTree->GetEntry(nEntry*gRandom->Rndm());
	  if (par.hltBit[67]==1) selectFlag=true;
       }
    }
    
    // Selection on Events
    
    // Fill reco vertex information
    tdata12.nv = par.nv+1;  
    tdata23.nv = par.nv+1;  
    tdata13.nv = par.nv+1;  
    for(int j = 1; j<par.nv;j++) {
       tdata12.vz[j+1] = par.vz[j];
       tdata23.vz[j+1] = par.vz[j];
       tdata13.vz[j+1] = par.vz[j];
       tdata12.vx[j+1] = par.vx[j];
       tdata23.vx[j+1] = par.vx[j];
       tdata13.vx[j+1] = par.vx[j];
       tdata12.vy[j+1] = par.vy[j];
       tdata23.vy[j+1] = par.vy[j];
       tdata13.vy[j+1] = par.vy[j];
    }
    // Fill MC vertex
    tdata12.vx[0] = par.vx[0];
    tdata23.vx[0] = par.vx[0];
    tdata13.vx[0] = par.vx[0];
    tdata12.vy[0] = par.vy[0];
    tdata23.vy[0] = par.vy[0];
    tdata13.vy[0] = par.vy[0];
    tdata12.vz[0] = par.vz[0];
    tdata23.vz[0] = par.vz[0];
    tdata13.vz[0] = par.vz[0];
    
    // add background 
    
    int bckHits=0;
    if (addL1Bck!=0||addL2Bck!=0||addL3Bck!=0){
       TF1 *fBck = new TF1("fBck","-0.00478376+0.000435517*x",0,1000);
       double val = fBck->Eval(par.nhits1)*2;
       for (int i = 0; i<par.nhits1;i++) {
          if (gRandom->Rndm()<val) {
	     bckHits++;
	  }
       }
       delete fBck;
    }
           
    if (addL1Bck!=0) {
       
       for (int i = par.nhits1; i < par.nhits1+bckHits; i++) {
          double eta, phi, r;
	  hLayer1Hit->GetRandom3(r,eta,phi);
	  par.eta1[i] = eta;
	  par.phi1[i] = phi;
	  par.r1[i] = r;
       }
       par.nhits1+=bckHits;
    }

    if (addL2Bck!=0) {
       //int bckHits = (int) (addL2Bck*gRandom->Rndm()+0.5);
       
       for (int i = par.nhits2; i < par.nhits2+bckHits; i++) {
          double eta, phi, r;
	  hLayer2Hit->GetRandom3(r,eta,phi);
	  par.eta2[i] = eta;
	  par.phi2[i] = phi;
	  par.r2[i] = r;
       }
       par.nhits2+=bckHits;
    }
    if (addL3Bck!=0) {
       //int bckHits = (int) (addL3Bck*gRandom->Rndm()+0.5);
       for (int i = par.nhits3; i < par.nhits3+bckHits; i++) {
          double eta, phi, r;
	  hLayer3Hit->GetRandom3(r,eta,phi);
	  par.eta3[i] = eta;
	  par.phi3[i] = phi;
	  par.r3[i] = r;
       }
       par.nhits3+=bckHits;
    }

    // add trackletVertex
    /*
    if (tdata12.nv == 2) tdata12.nv=3;
    if (tdata23.nv == 2) tdata23.nv=3;
    if (tdata13.nv == 2) tdata13.nv=3;
    */
    
    vector<RecoHit> layerRaw1;
    prepareHits(layerRaw1,par, cuts, 1, 0,0,0, splitProb, dropProb,  cutOnClusterSize, par.nRun, par.nLumi);
    vector<RecoHit> layerRaw2;
    prepareHits(layerRaw2,par, cuts, 2, 0,0,0, splitProb, dropProb,  cutOnClusterSize, par.nRun, par.nLumi);

    if (nPileUp!=0 ) {
       if (gRandom->Rndm()<nPileUp) {
          t->GetEntry(i+1);
          prepareHits(layerRaw1,par, cuts, 1, 0,0,0, splitProb, dropProb,  cutOnClusterSize, par.nRun, par.nLumi);
          prepareHits(layerRaw2,par, cuts, 2, 0,0,0, splitProb, dropProb,  cutOnClusterSize, par.nRun, par.nLumi);
       }
       t->GetEntry(i);
    }
    
    double trackletVertex = 0;

    // Choose KK Vertex if specified =============================================
    if (par.nhits1> vertexHitRegion|| useKKVertex) 
    { 
       // actually use pixel3vertex
       trackletVertex = par.vz[1]; 
    } else {
       trackletVertex = TrackletVertexUnbin(layerRaw1,layerRaw2,0.14,0.08);
    }
    // ===========================================================================
    
    // For particle gun =========
    //    if (i % 1000 == 0) cout <<"!!! USE GEN VERTEX (FOR PARTICLE GUN)"<<endl;
    //    trackletVertex = par.vz[0];
    //===========================

    // vz[1] is always the selected algorithm
    
   double smear =0;
   if (smearVertex!=0) {
       if (i==1) cout <<"Vertex smeared!"<<endl;
       while (smear!=0) {
          double x = gRandom->Rndm()*2-1;
          if (gRandom->Rndm()<TMath::Gaus(x,0,smearVertex,1)) {
             smear = x;
          }
       }
       trackletVertex += smear;
    }
    
    tdata12.vz[1] = trackletVertex;
    tdata23.vz[1] = trackletVertex;
    tdata13.vz[1] = trackletVertex;
    
    if(useKKVertex) {
       if (i==1) cout <<"Use Reconstructed Vertex "<<endl;
       tdata12.vx[1] = par.vx[1];
       tdata12.vy[1] = par.vy[1];
       tdata13.vx[1] = par.vx[1];
       tdata13.vy[1] = par.vy[1];
       tdata23.vx[1] = par.vx[1];
       tdata23.vy[1] = par.vy[1];
    } else {
       if (i==1)cout <<"Use Tracklet Vertex "<<endl;
    }

    if (useRandomVertex) {
       if (i==1) cout <<"Random Vertex!!!"<<endl;
           tdata12.vz[1]=gRandom->Rndm()*40-20;
           tdata13.vz[1]=tdata12.vz[1];
           tdata23.vz[1]=tdata12.vz[1];
   }
    // use trackletVertex
    if (fabs(tdata12.vz[1])>cuts.vzCut && makeVzCut == 1) continue;
    

    
    // Process hits with Vz constraint:
    vector<RecoHit> layer1;
    prepareHits(layer1,par, cuts, 1, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi);
    vector<RecoHit> layer2;
    prepareHits(layer2,par, cuts, 2, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi);
    vector<RecoHit> layer3;
    prepareHits(layer3,par, cuts, 3, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi);

    vector<RecoHit> layer1Cut;
    prepareHits(layer1Cut,par, cuts, 1, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, 1, par.nRun, par.nLumi);

    if (nPileUp!=0) {
       if (gRandom->Rndm()<nPileUp) {
          t->GetEntry(i+1);
          prepareHits(layer1,par, cuts, 1,tdata12.vx[1], tdata12.vy[1],  tdata12.vz[1], splitProb, dropProb, cutOnClusterSize);
          prepareHits(layer2,par, cuts, 2,tdata12.vx[1], tdata12.vy[1],  tdata12.vz[1], splitProb, dropProb, cutOnClusterSize);
          prepareHits(layer3,par, cuts, 3,tdata12.vx[1], tdata12.vy[1],  tdata12.vz[1], splitProb, dropProb, cutOnClusterSize);
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
    vector<Tracklet> protoTracklets12;
    vector<Tracklet> protoTracklets13;
    vector<Tracklet> protoTracklets23;
    vector<Tracklet> recoTracklets12;
    vector<Tracklet> recoTracklets13;
    vector<Tracklet> recoTracklets23;

    if (mimicPixelCounting) {
       protoTracklets12 = recoProtoTracklets(layer1,layer1);
       protoTracklets13 = recoProtoTracklets(layer2,layer2);
       protoTracklets23 = recoProtoTracklets(layer3,layer3);
    } else {
       protoTracklets12 = recoProtoTracklets(layer1,layer2);
       protoTracklets13 = recoProtoTracklets(layer1,layer3);
       protoTracklets23 = recoProtoTracklets(layer2,layer3);
    }  
    
    recoTracklets12  = cleanTracklets(protoTracklets12,0,cuts);
    recoTracklets13  = cleanTracklets(protoTracklets13,0,cuts);
    recoTracklets23  = cleanTracklets(protoTracklets23,0,cuts);
    

    // Move the Vertex back
    if (smearVertex!=0)
    {
       tdata12.vz[1] = trackletVertex-smear;
       tdata23.vz[1] = trackletVertex-smear;
       tdata13.vz[1] = trackletVertex-smear;
    }

    // Vertex Compatibility information
    float vtxQualCut=0;
    if (par.npxhits<150) {
       vtxQualCut = 1;
    } else {
       if (par.vtxqual>2) {
          vtxQualCut = 1;
       } else {
          if (par.vtxqual>0.0045*par.npxhits) {
	     vtxQualCut = 1;
 	  } else {
 	     vtxQualCut = 0;
 	  }
       }
    }
   
    // Fill Ntuple
    tdata12.nTracklet = recoTracklets12.size();
    tdata12.nhit1   = layer1.size();
    tdata12.nhit2   = layer2.size();
    tdata12.nRun    = par.nRun;
    tdata12.nEv     = par.nEv;
    tdata12.nLumi   = par.nLumi;
    tdata12.nBX     = par.nBX;
    tdata12.nHFn    = par.nHFp;
    tdata12.nHFp    = par.nHFn;
    tdata12.nHits   = layer1.size()+layer2.size()+layer3.size();
    tdata12.nHltBit = par.nHltBit;
    tdata12.nL1ABit = par.nL1ABit;
    tdata12.nL1TBit = par.nL1TBit;
    tdata12.vtxQualCut = vtxQualCut;
    tdata12.vtxqual = par.vtxqual;
    tdata12.npxhits = par.npxhits;
    
    for (int j=0;j<(int)par.nHltBit;j++)
    {
       tdata12.hltBit[j]= par.hltBit[j];
    }

    for (int j=0;j<(int)par.nL1ABit;j++)
    {
       tdata12.l1ABit[j]= par.l1ABit[j];
    }

    for (int j=0;j<(int)par.nL1TBit;j++)
    {
       tdata12.l1TBit[j]= par.l1TBit[j];
    }
    
    int ntracklet12s=0;
    int ntracklet12b=0;
    for (int j=0;j<(int)tdata12.nTracklet;j++)
    {
        tdata12.eta1[j] = recoTracklets12[j].eta1();	
        tdata12.eta2[j] = recoTracklets12[j].eta2();	
        tdata12.r1[j] = recoTracklets12[j].r1();	
        tdata12.r2[j] = recoTracklets12[j].r2();	
        tdata12.cs1[j] = recoTracklets12[j].cs1();	
        tdata12.cs2[j] = recoTracklets12[j].cs2();	
        tdata12.phi1[j] = recoTracklets12[j].phi1();	
        tdata12.phi2[j] = recoTracklets12[j].phi2();
	tdata12.deta[j] = recoTracklets12[j].deta();
   	tdata12.dphi[j] = recoTracklets12[j].dphi();
        if (fabs(tdata12.deta[j])<0.1){
	   if (fabs(tdata12.dphi[j])<1.0) {
	      ntracklet12s++;
	   } 
	   if (fabs(tdata12.dphi[j])>1.0&&fabs(tdata12.dphi[j])<2.0) {
	      ntracklet12b++;
	   }
	}
    }
    tdata12.mult=ntracklet12s-ntracklet12b;
    tdata12.mult2=layer1Cut.size();
    tdata12.npart=0;
    
    bool reWeightMultDropFlag = 0;       

    if (reweightMultiplicity) {
       reWeightMultDropFlag = 0;
       double myVz = tdata12.vz[1];       
       
       double Ratio = hRatio->GetBinContent(hRatio->FindBin(tdata12.mult));
       //cout <<Ratio<<endl;
       double x=gRandom->Rndm()*1.6;

       if (x> Ratio) reWeightMultDropFlag=1;
    }
    if (reWeightMultDropFlag) continue;

    
    double pro1=0;
    double pro2=0;
    for (int j=0;j<12;j++) tdata12.nhad[j]=0;

    for(int j=0;j<par.npart;j++)
    {
	if (fabs(par.pdg[j])==2212) {
	   double momentum = par.pt[j]*cosh(par.eta[j]);
	   if (momentum>pro2) {
	      if (momentum>pro1) {
	         pro2 = pro1;
		 pro1 = momentum;
	      } else {
	         pro2 = momentum;
	      }
//	      cout <<pro2<<" "<<pro1<<endl;
	   }
	}
        if (fabs(par.eta[j])>3||par.chg[j]==0||abs(par.pdg[j])==11||abs(par.pdg[j])==13) continue;
//        if (fabs(par.eta[j])>3) continue;
	tdata12.eta[tdata12.npart]=par.eta[j];
	tdata12.phi[tdata12.npart]=par.phi[j];
	tdata12.chg[tdata12.npart]=par.chg[j];
	tdata12.pdg[tdata12.npart]=par.pdg[j];
	if (abs(par.pdg[j]-11)<1) cout <<"Oh no???"<<endl;
	tdata12.pt[tdata12.npart]=par.pt[j];
        tdata12.npart++;
	int bin = (int)((par.eta[j]+3)*2);
	int pdg = (int)abs(par.pdg[j]);
	if (pdg== 211 || pdg == 321 || pdg == 2212 || pdg ==3122) tdata12.nhad[bin]++;
    }

    nhits->Fill(mult,layer1.size(),layer2.size());
    ntmult->Fill(mult,layer1.size(),layer2.size());

    tdata12.evtType = par.evtType;
    tdata12.pro2 = pro2;

    for (int j=0;j<par.nv;j++)
    {
       tdata12.vz[j]+=vzShift;
    }
    
    trackletTree12->Fill();

    tdata13.nTracklet = recoTracklets13.size();
    tdata13.nhit1 = layer1.size();
    tdata13.nhit2 = layer3.size();
    tdata13.nRun    = par.nRun;
    tdata13.nEv     = par.nEv;
    tdata13.nLumi   = par.nLumi;
    tdata13.nBX     = par.nBX;
    tdata13.nHFn    = par.nHFp;
    tdata13.nHFp    = par.nHFn;
    tdata13.nHits   = layer1.size()+layer2.size()+layer3.size();
    tdata13.nHltBit = par.nHltBit;
    tdata13.nL1ABit = par.nL1ABit;
    tdata13.nL1TBit = par.nL1TBit;
    tdata13.vtxQualCut = vtxQualCut;
    tdata13.vtxqual = par.vtxqual;
    tdata13.npxhits = par.npxhits;
    
    for (int j=0;j<(int)par.nHltBit;j++)
    {
       tdata13.hltBit[j]= par.hltBit[j];
    }

    for (int j=0;j<(int)par.nL1ABit;j++)
    {
       tdata13.l1ABit[j]= par.l1ABit[j];
    }

    for (int j=0;j<(int)par.nL1TBit;j++)
    {
       tdata13.l1TBit[j]= par.l1TBit[j];
    }

    int ntracklet13s=0;
    int ntracklet13b=0;
    for (int j=0;j<(int)tdata13.nTracklet;j++)
    {
        tdata13.eta1[j] = recoTracklets13[j].eta1();	
        tdata13.eta2[j] = recoTracklets13[j].eta2();	
        tdata13.r1[j] = recoTracklets13[j].r1();	
        tdata13.r2[j] = recoTracklets13[j].r2();	
        tdata13.cs1[j] = recoTracklets13[j].cs1();	
        tdata13.cs2[j] = recoTracklets13[j].cs2();	
        tdata13.phi1[j] = recoTracklets13[j].phi1();	
        tdata13.phi2[j] = recoTracklets13[j].phi2();
	tdata13.deta[j] = recoTracklets13[j].deta();
   	tdata13.dphi[j] = recoTracklets13[j].dphi();
        if (fabs(tdata13.deta[j])<0.1){
	   if (fabs(tdata13.dphi[j])<1.0) {
	      ntracklet13s++;
	   } 
	   if (fabs(tdata13.dphi[j])>1.0&&fabs(tdata13.dphi[j])<2.0) {
	      ntracklet13b++;
	   }
	}
    }
    tdata13.mult=ntracklet13s-ntracklet13b;
    tdata13.mult2=layer1Cut.size();
    tdata13.npart=0;
    
    for (int j=0;j<12;j++) tdata13.nhad[j]=0;

    for(int j=0;j<par.npart;j++)
    {
        if (fabs(par.eta[j])>3||par.chg[j]==0||fabs(par.pdg[j])==11||abs(par.pdg[j])==13) continue;
//        if (fabs(par.eta[j])>3) continue;
	tdata13.eta[tdata13.npart]=par.eta[j];
	tdata13.phi[tdata13.npart]=par.phi[j];
	tdata13.chg[tdata13.npart]=par.chg[j];
	tdata13.pdg[tdata13.npart]=par.pdg[j];
	tdata13.pt[tdata13.npart]=par.pt[j];
        tdata13.npart++;
	int bin = (int)((par.eta[j]+3)*2);
	int pdg = (int)abs(par.pdg[j]);
	if (pdg== 211 || pdg == 321 || pdg == 2213 || pdg ==3132) tdata13.nhad[bin]++;
    }

    nhits->Fill(mult,layer1.size(),layer2.size());
    ntmult->Fill(mult,layer1.size(),layer2.size());

    tdata13.evtType = par.evtType;
    tdata13.pro2 = pro2;

    for (int j=0;j<par.nv;j++)
    {
       tdata13.vz[j]+=vzShift;
    }

    trackletTree13->Fill();

    tdata23.nTracklet = recoTracklets23.size();
    tdata23.nhit1 = layer2.size();
    tdata23.nhit2 = layer3.size();
    tdata23.nRun    = par.nRun;
    tdata23.nEv     = par.nEv;
    tdata23.nLumi   = par.nLumi;
    tdata23.nBX     = par.nBX;
    tdata23.nHFn    = par.nHFp;
    tdata23.nHFp    = par.nHFn;
    tdata23.nHits   = layer1.size()+layer2.size()+layer3.size();
    tdata23.nHltBit = par.nHltBit;
    tdata23.nL1ABit = par.nL1ABit;
    tdata23.nL1TBit = par.nL1TBit;
    tdata23.vtxQualCut = vtxQualCut;
    tdata23.vtxqual = par.vtxqual;
    tdata23.npxhits = par.npxhits;
    
    for (int j=0;j<(int)par.nHltBit;j++)
    {
       tdata23.hltBit[j]= par.hltBit[j];
    }

    for (int j=0;j<(int)par.nL1ABit;j++)
    {
       tdata23.l1ABit[j]= par.l1ABit[j];
    }

    for (int j=0;j<(int)par.nL1TBit;j++)
    {
       tdata23.l1TBit[j]= par.l1TBit[j];
    }

    int ntracklet23s=0;
    int ntracklet23b=0;
    for (int j=0;j<(int)tdata23.nTracklet;j++)
    {
        tdata23.eta1[j] = recoTracklets23[j].eta1();	
        tdata23.eta2[j] = recoTracklets23[j].eta2();	
        tdata23.r1[j] = recoTracklets23[j].r1();	
        tdata23.r2[j] = recoTracklets23[j].r2();	
        tdata23.cs1[j] = recoTracklets23[j].cs1();	
        tdata23.cs2[j] = recoTracklets23[j].cs2();	
        tdata23.phi1[j] = recoTracklets23[j].phi1();	
        tdata23.phi2[j] = recoTracklets23[j].phi2();
	tdata23.deta[j] = recoTracklets23[j].deta();
   	tdata23.dphi[j] = recoTracklets23[j].dphi();
        if (fabs(tdata23.deta[j])<0.1){
	   if (fabs(tdata23.dphi[j])<1.0) {
	      ntracklet23s++;
	   } 
	   if (fabs(tdata23.dphi[j])>1.0&&fabs(tdata23.dphi[j])<2.0) {
	      ntracklet23b++;
	   }
	}
    }
    tdata23.mult=ntracklet23s-ntracklet23b;
    tdata23.mult2=layer1Cut.size();

    tdata23.npart=0;
    for (int j=0;j<12;j++) tdata23.nhad[j]=0;

    for(int j=0;j<par.npart;j++)
    {
        if (fabs(par.eta[j])>3||par.chg[j]==0||fabs(par.pdg[j])==11||abs(par.pdg[j])==13) continue;
//        if (fabs(par.eta[j])>3) continue;
	tdata23.eta[tdata23.npart]=par.eta[j];
	tdata23.phi[tdata23.npart]=par.phi[j];
	tdata23.chg[tdata23.npart]=par.chg[j];
	tdata23.pdg[tdata23.npart]=par.pdg[j];
	tdata23.pt[tdata23.npart]=par.pt[j];
        tdata23.npart++;
	int bin = (int)((par.eta[j]+3)*2);
	int pdg = (int)abs(par.pdg[j]);
	if (pdg== 211 || pdg == 321 || pdg == 2223 || pdg ==3232) tdata23.nhad[bin]++;
    }

    nhits->Fill(mult,layer1.size(),layer2.size());
    ntmult->Fill(mult,layer1.size(),layer2.size());

    tdata23.evtType = par.evtType;
    tdata23.pro2 = pro2;

    for (int j=0;j<par.nv;j++)
    {
       tdata23.vz[j]+=vzShift;
    }


    trackletTree23->Fill();
  }

  // Close outputfile ===================================================================================

  outf->Write();
  outf->Close(); 
}


