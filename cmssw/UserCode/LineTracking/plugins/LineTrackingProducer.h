#ifndef LineTrackingProducer_H
#define LineTrackingProducer_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TVector3.h"

typedef std::vector<int> Line;
typedef std::vector<Line> LineCollection;

typedef std::pair<std::pair<double,double>, std::vector<int> > LineVertex;
typedef std::vector<LineVertex> LineVertexCollection;

namespace edm { class Run; class Event; class EventSetup; }
namespace reco { class Track; class Vertex; }

class BaseTrackerRecHit;
using ConstRecHitPointer = BaseTrackerRecHit const *;

typedef std::pair<reco::Track *,
                  std::vector<ConstRecHitPointer> > TrackWithRecHits;

class TrackerGeometry;

class LineTrackingProducer : public edm::EDProducer
{
 public:
  explicit LineTrackingProducer(const edm::ParameterSet & ps);
  ~LineTrackingProducer();

  virtual void produce(edm::Event & ev, const edm::EventSetup & es);
 
 private:
  void beginRun(edm::Run const & run, edm::EventSetup const & es) override;

  void readHits(const edm::Event& ev);

  // bool sortByRadius(const TVector3 & a, const TVector3 & b);
  // bool sortByAbsZ  (const TVector3 & a, const TVector3 & b);

  void processHits(std::vector<TrackWithRecHits> & recTracksWithHits,
                   std::vector<reco::Vertex> & recVertexes);

  std::vector<TVector3> points;
  std::vector<std::pair<unsigned long int, unsigned long int> > detids;

  const TrackerGeometry * theTracker; 

  // when to stop track clustering into vertices
  int maxClusterWidthDiff;

  // when to stop track clustering into vertices
  double dMaxForVertexing;

  // number of tracking-vertexing rounds to do
  int nRounds;

  // maximal number of shared dets for two tracks
  int maxSharedDets;

  // for track hit closest to the beamline
  double maxFirstHitRadius;

  // for track z
  double maxAbsoluteZ0;

  // max number of clusters
  int maxClusters;
};
#endif
