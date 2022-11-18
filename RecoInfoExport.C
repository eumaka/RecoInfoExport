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
        
         
        stringstream spts;
        int t = 1;
        int c = -1;
        
        SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
        TrkrClusterContainer *clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
        ActsGeometry *geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
        if(!geometry)
        {
         std::cout << PHWHERE << "No Acts geometry on node tree. Can't  continue."
               << std::endl;
        }
        
        float x = 0.; float y = 0.; float z = 0.;
        float px = 0.; float py = 0.; float pz = 0.;
        TVector3 pos; TVector3 mom;
        
        bool first = true;

        for (SvtxTrackMap::Iter iter = trackmap->begin(); iter != trackmap->end(); ++iter)
        {
          SvtxTrack* track = iter->second;
          px = track->get_px();
          py = track->get_py();
          pz = track->get_pz();
          mom.SetX(px); mom.SetY(py); mom.SetZ(pz);
            
          std::vector<TrkrDefs::cluskey> clusters;
          auto siseed = track->get_silicon_seed();
          if(siseed)
          {
            for (auto iter = siseed->begin_cluster_keys(); iter != siseed->end_cluster_keys(); ++iter)
                {
                      TrkrDefs::cluskey cluster_key = *iter;
                      clusters.push_back(cluster_key);
                }
            }
            
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
            TrkrCluster* cluster = clusterContainer->findCluster(cluster_key);
            if(!cluster) continue;
            
            Acts::Vector3 globalClusterPosition = geometry->getGlobalPosition(cluster_key, cluster);
            x = globalClusterPosition(0);
            y = globalClusterPosition(1);
            z = globalClusterPosition(2);
            pos.SetX(x); pos.SetY(y); pos.SetZ(z);
                    
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
        }
            fdata
            << (boost::format(
              "{ \"pt\": %1%, \"t\": %2%, \"e\": %3%, \"p\": %4%, \"c\": %5%, \"pts\":[ %6% ]},")
               % mom.Pt() % t % mom.PseudoRapidity() % mom.Phi() % c % spts.str()) << endl;
               spts.clear();
               spts.str("");
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

