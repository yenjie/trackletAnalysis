#include "LineTrackingProducer.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h"

#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/BaseTrackerRecHit.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/Math/interface/Vector.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "../interface/Tracking.h"
#include "../interface/StraightLine.h"
#include "../interface/Vertexing.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>

#include "TVectorD.h"

using namespace std;

#define sqr(x) ((x) * (x))

/*****************************************************************************/
LineTrackingProducer::LineTrackingProducer(const edm::ParameterSet & ps)
{
  maxClusters         = ps.getParameter<int>("maxClusters");          // 6000
  maxClusterWidthDiff = ps.getParameter<int>("maxClusterWidthDiff");  // 2

  nRounds             = ps.getParameter<int>("nRounds");              // 2
  maxFirstHitRadius   = ps.getParameter<double>("maxFirstHitRadius"); // 31.
  maxSharedDets       = ps.getParameter<int>("maxSharedDets");        // 3
  maxAbsoluteZ0       = ps.getParameter<double>("maxAbsoluteZ0");     // 20.

  dMaxForVertexing    = ps.getParameter<double>("dMaxForVertexing");  // 50.
  
    produces<TrackingRecHitCollection>();
       produces<reco::TrackCollection>();
  produces<reco::TrackExtraCollection>();
     produces<reco::VertexCollection>();
}

/*****************************************************************************/
LineTrackingProducer::~LineTrackingProducer()
{
}

/*****************************************************************************/
void LineTrackingProducer::beginRun
  (edm::Run const & run, edm::EventSetup const & es)
{
  // Get tracker geometry
  edm::ESHandle<TrackerGeometry> trackerHandle;
  es.get<TrackerDigiGeometryRecord>().get(trackerHandle);
  theTracker = trackerHandle.product();
}

/*****************************************************************************/
void LineTrackingProducer::readHits(const edm::Event& ev)
{
  vector<edm::Handle<SiStripMatchedRecHit2DCollection> > stripColls;
  ev.getManyByType(stripColls);

  for(vector<edm::Handle<SiStripMatchedRecHit2DCollection> >::const_iterator
      stripColl = stripColls.begin();
      stripColl!= stripColls.end(); stripColl++)
  {
    const SiStripMatchedRecHit2DCollection* theStripHits =
                                             (*stripColl).product();

    for(SiStripMatchedRecHit2DCollection::DataContainer::const_iterator
            recHit = theStripHits->data().begin();
            recHit!= theStripHits->data().end(); recHit++)
    {
      DetId id = recHit->geographicalId();
      LocalPoint lpos = recHit->localPosition();
      GlobalPoint p = theTracker->idToDet(id)->toGlobal(lpos);

      // check if strip cluster widths are close
      if(abs(  recHit->monoHit().cluster()->amplitudes().size() -
             recHit->stereoHit().cluster()->amplitudes().size())
             <= maxClusterWidthDiff)
      {
        TVector3 hit(p.x(), p.y(), p.z());
        points.push_back(hit);

        detids.push_back(pair<unsigned long int, unsigned long int>(
            recHit->monoHit().geographicalId(),
          recHit->stereoHit().geographicalId()));
      }
    }
  }

  cerr << " read " << points.size() << " double-sided strip clusters" << endl;
}

/*****************************************************************************/
bool sortByRadius(const TVector3 & a, const TVector3 & b)
{
  return (a.Perp2() < b.Perp2());
}

/*****************************************************************************/
bool sortByAbsZ(const TVector3 & a, const TVector3 & b)
{
  return (fabs(a.z()) < fabs(b.z()));
}

/*****************************************************************************/
void LineTrackingProducer::processHits
  (vector<TrackWithRecHits> & recTracksWithHits,
   vector<reco::Vertex> & recVertexes)
{
  // Initial guess
  vector<double> z0(1,0.); // z0 = 0 cm

  if(points.size() > 0)
  for(int r = 0; r < nRounds; r++)
  {
    double dMax;
    unsigned int mHits;

    if(r == 0) { dMax = 0.10; mHits = 3; }
    if(r == 1) { dMax = 0.05; mHits = 4; }

    // Initialize clusters, at most maxVertices
    unsigned int maxVertices = points.size(); // FIXME

    //////////////////
    // fPNN
    unsigned int nOptimal;
    vector<vector<int> > lists;

    cerr << " search with z0 = (";
    for(vector<double>::const_iterator z = z0.begin();
                                       z!= z0.end(); z++)
      cerr << " " << *z;
    cerr << ") cm" << endl;

    Tracking theTracking(dMax, z0);
    theTracking.run(points, detids, nOptimal, lists, maxVertices);

    // Fill
    int K = nOptimal;
    LineCollection tracks;
    for(int k = 0; k < K; k++)
    {
      Line track;

      for(vector<int>::iterator il = lists[k].begin();
                                il!= lists[k].end(); il++)
        track.push_back(*il);

      tracks.push_back(track);
    }

    // Track filtering based on number of shared detids
    for(LineCollection::iterator track1 = tracks.begin();
                                 track1!= tracks.end(); track1++)
    for(LineCollection::iterator track2 = track1 + 1;
                                 track2!= tracks.end(); )
    {
//      int sub = -1;

      int nshared = 0;
      for(vector<int>::const_iterator i1 = track1->begin();
                                      i1!= track1->end(); i1++)
      for(vector<int>::const_iterator i2 = track2->begin();
                                      i2!= track2->end(); i2++)
        if(detids[*i1] == detids[*i2])
          nshared++;

      if(nshared > maxSharedDets) // remove track2, FIXME
      {
/*
        cerr << "  removed " << nshared << " " << sub
             << " " << track1->size() << " " << track2->size() << endl;
*/
        track2 = tracks.erase(track2);
      }
      else
        track2++;
    }

    // Look at vertices
    int n = 0;

    vector<pair<double,double> > vpoints;

//    vector<reco::TrackBaseRef> basetracks;

    for(LineCollection::const_iterator track = tracks.begin();
                                       track!= tracks.end(); track++)
    if(track->size() >= mHits)
    {
      // copy
      vector<TVector3> ps;
      for(vector<int>::const_iterator i = track->begin();
                                      i!= track->end(); i++)
        ps.push_back(points[*i]);
      sort(ps.begin(), ps.end(), sortByRadius);

      vector<TVector3> pz;
      for(vector<int>::const_iterator i = track->begin();
                                      i!= track->end(); i++)
        pz.push_back(points[*i]);
      sort(pz.begin(), pz.end(), sortByAbsZ);

      // have to have the first point close
      if(ps[0].Perp() < maxFirstHitRadius)
      {
        StraightLine theFitter;
        vector<double> x(4);
        double average_delta = theFitter.fit(ps, x);
       
        double z = x[0]; 
        double theta = x[2]; 
        double sigma_z = average_delta / sin(theta);

//        double eta = -log(tan(theta/2));
        double phi =  x[3];

        if(fabs(z) < maxAbsoluteZ0)
        {
          vpoints.push_back(pair<double,double>(z,sqr(sigma_z)));

          if(r ==  nRounds - 1)
          {
            reco::TrackBase::CovarianceMatrix cov;
  
            reco::TrackBase::TrackAlgorithm algo =
            reco::TrackBase::TrackAlgorithm::undefAlgorithm;
  
            reco::TrackBase::TrackQuality   qual =
            reco::TrackBase::TrackQuality::undefQuality;
  
            reco::Track recTrack(average_delta, // chi2
                       ps.size()-1,   // ndof
                       reco::TrackBase::Point(z,0,0),   // vertex
                       reco::TrackBase::Vector(sin(theta)*cos(phi),
                                               sin(theta)*sin(phi),
                                               cos(theta)), // momentum
                                 99,    // charge
                                 cov, algo, qual); // covar, algo, quality
  
            vector<ConstRecHitPointer> recHits;
            recTracksWithHits.push_back(TrackWithRecHits(&recTrack, recHits));
          }

          n++;
        }
      }
    }

    cerr << "  found    " << tracks.size() << " candidates" << endl;
    cerr << "  selected " << n << " tracks" << endl;

    // agglo vertexing
    if(r < nRounds - 1)
    if(n > 0)
    {
      // Initialize clusters, at most maxVertices
      const int maxVertices = 1000;
      vector<pair<TVectorD,TVectorD> > clusters;
      for(unsigned int i = 0; i <= maxVertices; i++)
      {
        TVectorD mu(i);
        TVectorD  P(i);

        clusters.push_back(pair<TVectorD,TVectorD>(mu,P));
      }

      Vertexing theVertexing(dMaxForVertexing);

      unsigned int nOptimal;
      vector<vector<int> > lists;
      theVertexing.run(vpoints, clusters, nOptimal, lists, maxVertices);

      //
      int K = nOptimal;

      TVectorD mu(K); mu = clusters[K].first;
      TVectorD P(K) ; P  = clusters[K].second;

      LineVertexCollection vertices;
      for(int k = 0; k < K; k++)
      {
        LineVertex vertex;

        double sig2 = 0;
        for(vector<int>::iterator il = lists[k].begin();
                                  il!= lists[k].end(); il++)
          sig2 += 1 / vpoints[*il].second;
    
        vertex.first = pair<double,double>(mu(k), 1/sig2);
    
        for(vector<int>::iterator il = lists[k].begin();
                                  il!= lists[k].end(); il++)
          vertex.second.push_back(*il);

        vertices.push_back(vertex);
      }

      //
      if(vertices.size() > 0)
      {
        z0.clear();

        for(LineVertexCollection::iterator vertex = vertices.begin();
                                           vertex!= vertices.end(); )
        {
          if(vertices.size() == 1 || vertex->second.size() >= 3)
          {
            z0.push_back(vertex->first.first);

            cerr << "  vertex with "
                 << vertex->second.size() << " tracks, at "
                 << vertex->first.first << " +/- "
                 << vertex->first.second << " cm" << endl;

            vertex++;
          }
          else
            vertex = vertices.erase(vertex);
        }
      }

      if(r == nRounds - 2)
      {
        cerr << " found vertices " << vertices.size() << endl;

        for(LineVertexCollection::iterator vertex = vertices.begin();
                                           vertex!= vertices.end(); vertex++)
        {
          reco::Vertex::Error err;
          err(2,2) = vertex->first.second;
          reco::Vertex recVertex(reco::Vertex::Point(0,0,vertex->first.first),
                                 err, 0, vertex->second.size(), 1);

          recVertexes.push_back(recVertex);
        }

        cerr << " vertices pushed" << endl;
      }
    }
  }
}

/*****************************************************************************/
void LineTrackingProducer::produce(edm::Event& ev, const edm::EventSetup& es)
{
  //
  auto_ptr<TrackingRecHitCollection>   recHits
      (new TrackingRecHitCollection());
  auto_ptr<reco::TrackCollection>      recTracks
      (new reco::TrackCollection());
  auto_ptr<reco::TrackExtraCollection> recTrackExtras
      (new reco::TrackExtraCollection());
  auto_ptr<reco::VertexCollection>     recVertices
      (new reco::VertexCollection());

  // clear
  points.clear();
  detids.clear();

  // read hits
  readHits(ev);

  // do it
  vector<TrackWithRecHits> recTracksWithHits;
  vector<reco::Vertex> recVertexes; // FIXME

  if(points.size() > 0 && int(points.size()) < maxClusters)
    processHits(recTracksWithHits, recVertexes);

  int cc = 0, nTracks = recTracksWithHits.size();

  // hits, track
  for(int i = 0; i < nTracks; i++)
  {
    reco::Track * recTrack = recTracksWithHits[i].first;

    for(vector<ConstRecHitPointer>::const_iterator 
           h = recTracksWithHits[i].second.begin();
           h!= recTracksWithHits[i].second.end(); h++)
    {
      TrackingRecHit * recHit = (*h)->hit()->clone();

      recTrack->appendHitPattern(*recHit);
      recHits->push_back(recHit);
    }
    recTracks->push_back(*recTrack);

//    delete recTrack;
  }

  // store
  edm::OrphanHandle<TrackingRecHitCollection> ohRH = ev.put(recHits);
  edm::RefProd<TrackingRecHitCollection> hitCollProd(ohRH);

  // track extras
  for(int i = 0; i < nTracks; i++)
  {
    reco::TrackExtra theTrackExtra{};

    // fill the TrackExtra with TrackingRecHitRef
    unsigned int nHits = recTracks->at(i).numberOfValidHits();
    theTrackExtra.setHits(hitCollProd, cc, nHits);
    cc += nHits;
    recTrackExtras->push_back(theTrackExtra);
  }

  // store
  edm::OrphanHandle<reco::TrackExtraCollection> ohTE = ev.put(recTrackExtras);

  // add ref to track extras
  for(int i = 0; i < nTracks; i++)
  {
    const reco::TrackExtraRef theTrackExtraRef(ohTE,i);
    (recTracks->at(i)).setExtra(theTrackExtraRef);
  }

  // store
  ev.put(recTracks);

  // vertices
  int nVertices = recVertexes.size();
  for(int i = 0; i < nVertices; i++)
    recVertices->push_back(recVertexes[i]);

  // store
  ev.put(recVertices);
}

