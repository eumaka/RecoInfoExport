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
        0)
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
      // Grab the towers
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
      for (const auto & tower : good_towers)
        {
          assert(tower);

          float eta = towergeom->get_etacenter(tower->get_bineta());
          float phi = towergeom->get_phicenter(tower->get_binphi());

          phi = atan2(cos(phi), sin(phi));

          if (first)
            {
              first = false;
            }
          else
            fdata << ",";

           fdata
              << (boost::format(
                    "{ \"eta\": %1%, \"phi\": %2%, \"e\": %3%}") % eta % phi % tower->get_energy()) << endl;

        }
        fdata << "]" << endl;
        if(calo_name != "HCALOUT") fdata << "," << endl;
        fdata << endl;
    }

    {
        fdata << "}," << endl;
        fdata << "\"TRACKS\": {" << endl;
        fdata << "\"TRUTHINFO\": [" << endl;
        
         
        stringstream spts;
        int t = 1;
        int c = -1;
        bool first = true;
        
        SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
        TrkrClusterContainer* clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
        if(!clustermap)
        clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
        ActsGeometry *tgeometry = findNode::getClass<ActsGeometry>(topNode,"ActsGeometry");
        if(!tgeometry)
        {
         std::cout << PHWHERE << "No Acts geometry on node tree. Can't  continue."
               << std::endl;
        }
        
        for (SvtxTrackMap::Iter iter = trackmap->begin(); iter != trackmap->end(); ++iter)
        {
          SvtxTrack* track = iter->second;
          std::vector<TrkrDefs::cluskey> clusters;
          auto tpcseed = track->get_tpc_seed();
          if(tpcseed)
          {
            for (auto iter = tpcseed->begin_cluster_keys(); iter != tpcseed->end_cluster_keys(); ++iter)
            {
               TrkrDefs::cluskey cluster_key = *iter;
               clusters.push_back(cluster_key);
            }
          }

        for(unsigned int iclus = 0; iclus < clusters.size(); ++iclus)
        {
          TrkrDefs::cluskey cluster_key = clusters[iclus];
          TrkrCluster* cluster = clustermap->findCluster(cluster_key);
          if(!cluster) continue;
          Acts::Vector3 glob = tgeometry->getGlobalPosition(cluster_key, cluster);
          float x = glob(0);
          float y = glob(1);
          float z = glob(2);
          TVector3 pos(x, y, z);

          float px = track->get_px();
          float py = track->get_py();
          float pz = track->get_pz();
          TVector3 mom(px, py, pz);       
         
          if (first)
          {
             first = false;
          }
            
          else
            
          spts << ",";
          spts << "[";
          spts << pos.x();
          spts << ",";
          spts << pos.y();
          spts << ",";
          spts << pos.z();
          spts << "]";
          fdata << (boost::format(
                    "{ \"pt\": %1%, \"t\": %2%, \"e\": %3%, \"p\": %4%, \"c\": %5%, \"pts\":[ %6% ]}")
                    % mom.Perp() % t % pos.Eta() % pos.Phi() % c
                    % spts.str()) << endl;
        }
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
