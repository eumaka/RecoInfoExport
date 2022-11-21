#include "RecoInfoExport.h"

#include <phool/getClass.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/PHCompositeNode.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4Hit.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterv3.h>
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrClusterIterationMapv1.h>

#include <trackbase_historic/TrackSeed_v1.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/SvtxTrackSeed_v1.h>

#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack.h>

#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxTrackEval.h>
#include <g4eval/SvtxVertexEval.h>
#include <g4eval/SvtxTruthEval.h>

#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawCluster.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TDatabasePDG.h>
#include <TVector3.h>

#include <iostream>
#include <cassert>
#include <fstream>
#include <sstream>
#include <boost/format.hpp>
#include <boost/math/special_functions/sign.hpp>

using namespace std;

RecoInfoExport::RecoInfoExport(const string &name) :
    SubsysReco(name), _event(0), _calo_names(
      { "CEMC", "HCALIN", "HCALOUT" }), _tower_threshold(0), _pT_threshold(0), _min_track_hit_dist(
        0), _tower_PhiShift(0)
{
}

int
RecoInfoExport::Init(PHCompositeNode *topNode)
{

  return 0;
}

int
RecoInfoExport::process_event(PHCompositeNode *topNode)
{
  ++_event;

  stringstream fname;
  fname << _file_prefix << "event" << _event << ".json";
  fstream fdata(fname.str(), ios_base::out);
  fdata << "{" << endl;
  fdata << "\"HITS\": {" << endl;
     
  for (auto & calo_name : _calo_names)
    {
      string towernodename = "TOWER_CALIB_" + calo_name;
      RawTowerContainer* towers = findNode::getClass<RawTowerContainer>(topNode,
          towernodename.c_str());
      if (!towers)
        {
          std::cout << PHWHERE << ": Could not find node "
              << towernodename.c_str() << std::endl;
          return Fun4AllReturnCodes::ABORTRUN;
        }
      string towergeomnodename = "TOWERGEOM_" + calo_name;
      RawTowerGeomContainer *towergeom = findNode::getClass<
          RawTowerGeomContainer>(topNode, towergeomnodename.c_str());
      if (!towergeom)
        {
          cout << PHWHERE << ": Could not find node "
              << towergeomnodename.c_str() << endl;
          return Fun4AllReturnCodes::ABORTRUN;
        }

      set<const RawTower *> good_towers;
        
      RawTowerContainer::ConstRange begin_end = towers->getTowers();
      RawTowerContainer::ConstIterator rtiter;
      for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          const RawTower *tower = rtiter->second;
          assert(tower);
          if (tower->get_energy() > _tower_threshold)
            {
              good_towers.insert(tower);
            }
        }
        
      fdata <<"\""<<calo_name <<"\": [";

      bool first = true;
      const double geom_ref_radius = 225.87;
      for (const auto & tower : good_towers)
        {
          assert(tower);

          float eta = towergeom->get_etacenter(tower->get_bineta());
          float phi = towergeom->get_phicenter(tower->get_binphi());
            
          const double x(geom_ref_radius * cos(phi));
          const double y(geom_ref_radius * sin(phi));

           if (first)
            {
              first = false;
            }
          else
            fdata << ",";

           fdata
              << (boost::format(
                    "{ \"eta\": %1%, \"phi\": %2%, \"e\": %3%, \"x\": %4% , \"y\": %5%}") % eta % phi % tower->get_energy()) % x % y << endl;

        }
        fdata << "]" << endl;
        if(calo_name != "HCALOUT") fdata << "," << endl;
        fdata << endl;
    }

    {
        fdata << "}," << endl;
        fdata << "\"TRACKS\": {" << endl;
        fdata << "\"RECOTRACKS\": [" << endl;
             
        TrackSeedContainer *_track_map = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
        if(!_track_map)
        {
         std::cout << PHWHERE << "No TpcTrackSeedContainer on node tree. Can't  continue."
               << std::endl;
        }
        
        TrackSeedContainer *_track_map_silicon = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
        if(!_track_map_silicon)
        {
         std::cout << PHWHERE << "No SiliconTrackSeedContainer on node tree. Can't  continue."
               << std::endl;
        }
        
        TrkrClusterContainer *clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
        if(!clusterContainer)
        {
         std::cout << PHWHERE << "No TRKR_CLUSTER on node tree. Can't  continue."
               << std::endl;
        }
        
        ActsGeometry *geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
        if(!geometry)
        {
         std::cout << PHWHERE << "No Acts geometry on node tree. Can't  continue."
               << std::endl;
        }
        
        double _eta = 0.; double _phi = 0.; double _pt = 0.;
        double _eta_si = 0.; double _phi_si = 0.; double _pt_si = 0.;

        double x = 0.; double y = 0.; double z = 0.;
        double x_si = 0.; double y_si = 0.; double z_si = 0.;

        TVector3 pos; TVector3 pos_si;
        
        bool first = true;
        stringstream spts;
        
        int t = 1;
        int c = -1;
        
        
        for (unsigned int phtrk_iter_si = 0; phtrk_iter_si < _track_map_silicon->size(); ++phtrk_iter_si)
        {
            TrackSeed  *_tracklet_si = _track_map_silicon->get(phtrk_iter_si);
            if(!_tracklet_si) { continue; }
            _eta_si = _tracklet_si->get_eta();
            _phi_si = _tracklet_si->get_phi(clusterContainer,geometry);
            _pt_si = _tracklet_si->get_pt();
            
            for (auto iter_si = _tracklet_si->begin_cluster_keys(); iter_si != _tracklet_si->end_cluster_keys(); ++iter_si)
            {
                TrkrDefs::cluskey cluster_key_si = *iter_si;
                TrkrCluster* cluster_si = clusterContainer->findCluster(cluster_key_si);
                if(!cluster_si) continue;
                Acts::Vector3 globalClusterPosition_si = geometry->getGlobalPosition(cluster_key_si, cluster_si);
                x_si = globalClusterPosition_si(0);
                y_si = globalClusterPosition_si(1);
                z_si = globalClusterPosition_si(2);

            
            pos_si.SetX(x_si); pos_si.SetY(y_si); pos_si.SetZ(z_si);
            
            if (first)
            {
                first = false;
            }
            
            else
            spts << ",";
            spts << "[";
            spts << pos_si.x();
            spts << ",";
            spts << pos_si.y();
            spts << ",";
            spts << pos_si.z();
            spts << "]";
                
            }
            
            fdata
            << (boost::format(
              "{ \"pt\": %1%, \"t\": %2%, \"e\": %3%, \"p\": %4%, \"c\": %5%, \"pts\":[ %6% ]},")
               % _pt_si % t % _eta_si % _phi_si % c % spts.str()) << endl;
               spts.clear();
               spts.str("");
        }

        
        
        
        bool second = true;
        stringstream sptss;

        for (unsigned int phtrk_iter = 0; phtrk_iter < _track_map->size(); ++phtrk_iter)
        {
            TrackSeed *_tracklet_tpc = _track_map->get(phtrk_iter);
            if(!_tracklet_tpc) {continue;}
            _eta = _tracklet_tpc->get_eta();
            _phi = _tracklet_tpc->get_phi(clusterContainer,geometry);
            _pt = _tracklet_tpc->get_pt();
            
            for (auto iter = _tracklet_tpc->begin_cluster_keys(); iter != _tracklet_tpc->end_cluster_keys(); ++iter)
            {
                TrkrDefs::cluskey cluster_key = *iter;
                TrkrCluster* cluster = clusterContainer->findCluster(cluster_key);
                if(!cluster) continue;
                Acts::Vector3 globalClusterPosition = geometry->getGlobalPosition(cluster_key, cluster);
                x = globalClusterPosition(0);
                y = globalClusterPosition(1);
                z = globalClusterPosition(2);

            
            pos.SetX(x); pos.SetY(y); pos.SetZ(z);
            
            if (second)
            {
                second = false;
            }
            
            else
            sptss << ",";
            sptss << "[";
            sptss << pos.x();
            sptss << ",";
            sptss << pos.y();
            sptss << ",";
            sptss << pos.z();
            sptss << "]";
                
            }

            fdata
            << (boost::format(
              "{ \"pt\": %1%, \"t\": %2%, \"e\": %3%, \"p\": %4%, \"c\": %5%, \"pts\":[ %6% ]},")
               % _pt % t % _eta % _phi % c % sptss.str()) << endl;
               sptss.clear();
               sptss.str("");
        }
        
   }

   fdata << "]" << endl;
   fdata << "}}" << endl;
   fdata.close();
   return 0;
 
}


int
RecoInfoExport::End(PHCompositeNode *topNode)
{
  return 0;
}

