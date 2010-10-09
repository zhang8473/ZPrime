// -*- C++ -*-

// Package:    Analyzer
// Class:      Analyzer
// 
/**\class Analyzer Analyzer.cc Analyzer/Analyzer/src/Analyzer.cc

 Description: <one line class summary>
 It will save the information about primary vertex, HLT, muons with the original cut (GlobalMuon&&pt>20Gev/c). The muon objects will be saved in descending order of pt.
 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jinzhong Zhang
//         Created:  Wed Sep  9 18:30:00 CEST 2009
// $Id: Analyzer.cc,v 1.8 2010/05/14 16:44:01 zhangjin Exp $
//

#include "../interface/Analyzer.h"

using namespace edm;
using namespace std;
using namespace reco;

#define ReportError(code,MSG) Error.ErrorCode=code;    \
  Error.Run=Info.RUN;			      \
  Error.Event=Info.EVENT;		      \
  ErrorMsg_Tree->Fill();		      \
  cerr<<MSG<<endl;			      

#define trackIdLink(TrkID)  vector<SimTrack>::iterator Trk_iter = STC.begin(); \
  for (; Trk_iter != STC.end(); ++Trk_iter )	\
    if ((int) Trk_iter->trackId() == (int) TrkID) break;		\
  if (Trk_iter==STC.end()) {ReportError(2,"parentIndex/trackId Error")}	\
  thisTrk=&*Trk_iter;

#define GenSimMomentum(ppt,peta,pphi,ppdg) Gen_pdgId->push_back(ppdg);Gen_pt->push_back(ppt);Gen_eta->push_back(peta);Gen_phi->push_back(pphi)

#define GenSimVertex(vtx,vty,vtz,vtt) Gen_vx->push_back(vtx);Gen_vy->push_back(vty);Gen_vz->push_back(vtz);Gen_vt->push_back(vtt)

#define RecordSimTrack(thisTrk) ParToSim.push_back(thisTrk);			\
  IsParInHep->push_back(false);						\
  GenSimMomentum(thisTrk->momentum().pt(),thisTrk->momentum().eta(),thisTrk->momentum().phi(),thisTrk->type());	\
  if (!thisTrk->noVertex()) {						\
    SimVertex thisVtx=SVC[thisTrk->vertIndex()];			\
    GenSimVertex(thisVtx.position().x(),thisVtx.position().y(),thisVtx.position().z(),thisVtx.position().t()*1e9); \
  }									\
  else {GenSimVertex(0,0,0,0);}

#define RecordHepMC(GenParticle) ParToHep.push_back(GenParticle);	\
  IsParInHep->push_back(true);						\
  GenSimMomentum(sqrt(GenParticle->momentum().px()*GenParticle->momentum().px()+GenParticle->momentum().py()*GenParticle->momentum().py()),GenParticle->momentum().eta(),GenParticle->momentum().phi(),GenParticle->pdg_id()); \
  HepMC::GenVertex *thisVtx=GenParticle->production_vertex();		\
  if (thisVtx) {GenSimVertex(thisVtx->position().x()/10.,thisVtx->position().y()/10.,thisVtx->position().z()/10.,thisVtx->position().t()/299.792458);} \
  else {GenSimVertex(0,0,0,0);}

// ------------ method called on each new Event  ------------

bool
DiMuonSelector::filter(edm::Event& event, const edm::EventSetup& iSetup) {
  Summarization.Total_Events++;
  ClearVecs();
  //Set Muons Handle
  Handle<reco::MuonCollection> Muon;
  event.getByLabel("muons", Muon);
  if (!Muon.isValid()) return false;
  reco::MuonCollection const & muons = *Muon;
  vector<reco::MuonCollection::const_iterator> MuonQueue;
  bool OneMuonPtPassedThreshold=false;
  for ( reco::MuonCollection::const_iterator iter = muons.begin(); iter != muons.end() ; ++iter)
    if (iter->isTrackerMuon()||iter->isGlobalMuon())
      if (iter->pt()>MuonPtSelection) {
	MuonQueue.push_back(iter);
	if (iter->pt()>MuonHighestPtSelection) OneMuonPtPassedThreshold=true;
      }
  
#ifdef HepMCGenParticle
  Handle<edm::HepMCProduct> HepMCH;
  event.getByLabel(HepMCTag, HepMCH);
  HepMC::GenEvent HepGenEvent(*(HepMCH->GetEvent()));
  //HepMC::GenEvent * HepGenEvent = new HepMC::GenEvent(*(HepMCH->GetEvent())); //Memory leak!!!
  Byte_t num_Gen_Zp_Mu=0,num_Gen_Zp_AntiMu=0;
  Byte_t num_Gen_Z_Mu=0,num_Gen_Z_AntiMu=0;
  for (HepMC::GenEvent::particle_iterator GenParticle_iter = HepGenEvent.particles_begin(); GenParticle_iter != HepGenEvent.particles_end();GenParticle_iter++) 
    if (abs((*GenParticle_iter)->pdg_id())==13) { 
      HepMC::GenVertex *thisVtx=(*GenParticle_iter)->production_vertex();
      if (thisVtx) for (HepMC::GenVertex::particles_in_const_iterator pgenD = thisVtx->particles_in_const_begin(); pgenD != thisVtx->particles_in_const_end(); ++pgenD)
	if ((*pgenD)->pdg_id()==32||(*pgenD)->pdg_id()==22||(*pgenD)->pdg_id()==23) {//if the parent particle is Z, gamma* or Z' record
	  vector<HepMC::GenParticle *>::iterator ParToHep_iter=find(ParToHep.begin(),ParToHep.end(),*pgenD);
	  if (ParToHep_iter==ParToHep.end()) {RecordHepMC((*pgenD))}
	  if ((*GenParticle_iter)->pdg_id()==13) 
	    switch ((*pgenD)->pdg_id()) {
	    case 32:
	      num_Gen_Zp_Mu++;
	      break;
	    case 23:
	      num_Gen_Z_Mu++;
	      break;
	    }
	  if ((*GenParticle_iter)->pdg_id()==-13)
	    switch ((*pgenD)->pdg_id()) {
	    case 32:
	      num_Gen_Zp_AntiMu++;
	      break;
	    case 23:
	      num_Gen_Z_AntiMu++;
	      break;
	    }
	  RecordHepMC((*GenParticle_iter))
	}
    }//search for Z'->MuMu in MonteCarlo
  if ((MuonQueue.size()<2||!OneMuonPtPassedThreshold)&&((!num_Gen_Zp_Mu)&&(!num_Gen_Zp_AntiMu))&&((!num_Gen_Z_Mu)&&(!num_Gen_Z_AntiMu))) return false;
  Info.num_ZpMuMuInMC=num_Gen_Zp_Mu<num_Gen_Zp_AntiMu?num_Gen_Zp_Mu:num_Gen_Zp_AntiMu;
  Info.num_ZMuMuInMC=num_Gen_Z_Mu<num_Gen_Z_AntiMu?num_Gen_Z_Mu:num_Gen_Z_AntiMu;
  for (unsigned int i=0;i<4;i++)
    HepMCFourVec[i]=999999.;
  HepMC::GenEvent::particle_iterator GenParticle_iter = HepGenEvent.particles_begin();
  for (;GenParticle_iter != HepGenEvent.particles_end();GenParticle_iter++) { 
    HepMC::GenVertex *thisVtx=(*GenParticle_iter)->production_vertex();
    if (thisVtx) {
      HepMCFourVec[0]=thisVtx->position().x()/10.;
      HepMCFourVec[1]=thisVtx->position().y()/10.;
      HepMCFourVec[2]=thisVtx->position().z()/10.;
      HepMCFourVec[3]=thisVtx->position().t()/299.792458;
      break;
    }
  }
  if (GenParticle_iter == HepGenEvent.particles_end()) {ReportError(7,"No HepMC(Core) Vertex Information")}
#else
#endif
  if (MuonQueue.size()<2||!OneMuonPtPassedThreshold) return false;

  for (vector<reco::MuonCollection::const_iterator>::iterator iter_queue=MuonQueue.begin();iter_queue!=MuonQueue.end(); ++iter_queue )
    for (vector<reco::MuonCollection::const_iterator>::iterator iter_queue2=iter_queue;iter_queue2!=MuonQueue.end(); ++iter_queue2 )
      if ((*iter_queue)->pt()<(*iter_queue2)->pt()) {
	reco::MuonCollection::const_iterator temp=*iter_queue2;
	*iter_queue2=*iter_queue;
	*iter_queue=temp;
      }//save from the biggest pt

  //general event information
  Info.RUN   = event.id ().run ();
  Info.EVENT = event.id ().event();
  Info.LS    = event.luminosityBlock ();
  Info.ORBIT = event.orbitNumber ();
  Info.BX = event.bunchCrossing ();
  Info.isRealData = event.isRealData();

  Handle< ValueMap<reco::DeDxData> > dEdxTrackHandle;
  event.getByLabel(dEdxTag, dEdxTrackHandle);
  const ValueMap<reco::DeDxData> dEdxTrack = *(dEdxTrackHandle.product());

  //Reco Tracks
  Handle<reco::TrackCollection> trackCollectionH;
  event.getByLabel(tracksTag,trackCollectionH);
  reco::TrackCollection tC = *trackCollectionH.product();
  Handle< View<Track> > trackCollectionHV;
  event.getByLabel(tracksTag,trackCollectionHV);

  //HLT Information
  Handle<TriggerResults> hltTriggerResults;
  event.getByLabel(TriggerResultsTag,hltTriggerResults);
  if (hltTriggerResults.isValid()) {
    if (FirstEntry) {
      string HLTTriggerNames;
      const TriggerNames &HLTNames = event.triggerNames(*hltTriggerResults);
      HLTSize = hltTriggerResults->size();
      if (HLTSize>maxHLTnum) throw cms::Exception("You need to increase the constant maxHLTnum.\n");
      HLTTriggerNames=HLTNames.triggerName(0)+"/O:";
      unsigned int i=1;
      for (; i < HLTSize-1; i++)
	HLTTriggerNames+=HLTNames.triggerName(i)+":";
      HLTTriggerNames+=HLTNames.triggerName(i);
      Muons_Tree->Branch("HLTAcceptance",HLTacceptance,HLTTriggerNames.c_str());
    }
    for (unsigned int i = 0 ; i < HLTSize ; i++) 
      HLTacceptance[i]=hltTriggerResults->accept(i);
  }
  else if (HLTSize>0) {
    ReportError(8,"Invalid TriggerResultsTag: "<<TriggerResultsTag.label())
    for (unsigned int i = 0 ; i < HLTSize ; i++) 
      HLTacceptance[i]=0;
  }
  
  //HLT Objects
  Handle<trigger::TriggerEvent> trgEvent;
  event.getByLabel(TriggerEventTag,trgEvent);
  if (trgEvent.isValid()) {
    const trigger::TriggerObjectCollection& TOC = trgEvent->getObjects();
    int Total_Filters = trgEvent->sizeFilters();
    for( unsigned int i =0; i < num_HLTsSaveObjs;i++) {
      trigger::size_type pos_Filter = trgEvent->filterIndex(HLTFilterNames[i]);
      if (pos_Filter < Total_Filters) {
	const trigger::Keys& KEYS(trgEvent->filterKeys(pos_Filter));
	unsigned int num_thisHLTObjs=KEYS.size();
	for(unsigned int j = 0; j < num_thisHLTObjs; j++) {
	  const trigger::TriggerObject& TO = TOC[KEYS[j]];
	  //printf("pt%f\teta%f\tphi%f\n",TO.pt(),TO.eta(),TO.phi());
	  HLTObj_pt[i]->push_back(TO.pt());
	  HLTObj_eta[i]->push_back(TO.eta());
	  HLTObj_phi[i]->push_back(TO.phi());
	}
      }
      else {ReportError(9,"FilterName "<<HLTFilterNames[i].label()<<" is not valid.")}
    }
  }
  else {ReportError(10,"TriggerEventTag "<<TriggerEventTag.label()<<" is not valid.")};

  //Primary Vertex
  Handle<reco::VertexCollection> recVtxs;
  event.getByLabel(PrimaryVerticesTag.c_str(),recVtxs);
  if (recVtxs.isValid())
    for(reco::VertexCollection::const_iterator v=recVtxs->begin(); v!=recVtxs->end(); ++v) {
      vx->push_back(v->x());	   vxError->push_back(v->xError());
      vy->push_back(v->y());	   vyError->push_back(v->yError());
      vz->push_back(v->z());	   vzError->push_back(v->zError());
    }
  else {ReportError(1,PrimaryVerticesTag.c_str()<<" information is not valid.")}
  
#ifdef TrackingParticles //Simulated Vertices: the vertexId() is just the position
  Handle<SimVertexContainer> SVCollectionH;
  event.getByLabel("g4SimHits", SVCollectionH);
  SVC = *SVCollectionH.product();
  
  //Simulated Tracks: the trackId() is not the position
  Handle<SimTrackContainer> STCollectionH;
  event.getByLabel("g4SimHits", STCollectionH);
  SimTrackContainer STC = *STCollectionH.product();
  
  //find daughters of each SimTrack
  Daughters.clear();
  for (vector<SimTrack>::iterator FindTrack = STC.begin(); FindTrack != STC.end(); ++FindTrack )
    if (!FindTrack->noVertex()) {
      SimVertex thisVtx=SVC[FindTrack->vertIndex()];
      if (!thisVtx.noParent()) {
	SimTrack *thisTrk;
	trackIdLink(thisVtx.parentIndex())
	Daughters[thisTrk].push_back(&*FindTrack);
      }
    }
 
  Handle<TrackingParticleCollection> TPCollectionH ;
  event.getByType(TPCollectionH);
  
  //SimToReco Tracks Association
  ESHandle<TrackAssociatorBase> AssociatorByHits;
  iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits", AssociatorByHits);
  //SimToRecoCollection SimToRecoByHits = AssociatorByHits->associateSimToReco(trackCollectionHV,TPCollectionH,&event);
  RecoToSimCollection RecoToSimByHits = AssociatorByHits->associateRecoToSim(trackCollectionHV,TPCollectionH,&event,&iSetup);
  
  //Match by chi2
  ESHandle<TrackAssociatorBase> AssociatorByChi2;
  iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByChi2", AssociatorByChi2);
  //SimToRecoCollection SimToRecoByChi2 = AssociatorByChi2->associateSimToReco(trackCollectionHV,TPCollectionH,&event);
  RecoToSimCollection RecoToSimByChi2 = AssociatorByChi2->associateRecoToSim(trackCollectionHV,TPCollectionH,&event,&iSetup);
#endif

  for (vector<reco::MuonCollection::const_iterator>::iterator iter_queue=MuonQueue.begin();iter_queue!=MuonQueue.end(); ++iter_queue ) {
    //muon loop begin
    reco::MuonCollection::const_iterator iter=*iter_queue;
    //muon basic information (pt,eta,phi,charge)
    pt->push_back(iter->pt());
    eta->push_back(iter->eta());
    phi->push_back(iter->phi());
    if (iter->charge()==-1) chargeMinus->push_back(true);
    else chargeMinus->push_back(false);
    //Muon Selectors
    if (iter->isTrackerMuon()) Summarization.Total_TrackerMuons++;
    isTrackerMu->push_back(iter->isTrackerMuon());
    if (iter->isGlobalMuon()) Summarization.Total_GlobalMuon++;
    isGlobalMu->push_back(iter->isGlobalMuon());
    if (!iter->isTrackerMuon()&&iter->isGlobalMuon()) Summarization.Total_GlobalnotTrackerMuon++;
    for (unsigned int whichcut=0;whichcut<num_Cuts;whichcut++)
      SelectorPassed[whichcut]->push_back(muon::isGoodMuon(*iter,OfficialMuonSelectors[whichcut]));
    MySelector->push_back(muon::isGoodMuon(*iter,muon::TMLastStation,1,3,3,3,3,maxChamberDist,maxChamberDistPull,MuonArbitrationType,false,true));
    //Station and Segment Matches
    StationMask->push_back(iter->stationMask(MuonArbitrationType));
    RequiredStationMask->push_back(muon::RequiredStationMask(*iter,maxChamberDist,maxChamberDistPull,MuonArbitrationType));
    for(int stationIdx = 0; stationIdx <8; ++stationIdx) {
      TrackDistToChamberEdge->push_back(iter->trackDist(stationIdx%4+1,stationIdx/4+1,MuonArbitrationType));//999999 means that there is no track
      TrackDistToChamberEdgeErr->push_back(iter->trackDistErr(stationIdx%4+1,stationIdx/4+1,MuonArbitrationType));//999999 means that there is no track
      DXTrackToSegment->push_back(999999);
      DYTrackToSegment->push_back(999999);
      DXErrTrackToSegment->push_back(999999);
      DYErrTrackToSegment->push_back(999999);
      DRTrackToSegment->push_back(999999);
      DRErrTrackToSegment->push_back(999999);
      NumberOfHitsInSegment->push_back(0);
      IsSegmentBelongsToTrackByDR->push_back(false);
      IsSegmentBelongsToTrackByCleaning->push_back(false);
    }
    for( std::vector<MuonChamberMatch>::const_iterator chamberMatch = iter->matches().begin();chamberMatch != iter->matches().end(); chamberMatch++ ) {
      if (chamberMatch->segmentMatches.empty()) continue;
      Byte_t detectorIdx=0;
      if (chamberMatch->detector()==MuonSubdetId::CSC) detectorIdx=1;
      Byte_t SavePostion=(eta->size()-1)*8+detectorIdx*4+chamberMatch->station()-1;
      for( std::vector<MuonSegmentMatch>::const_iterator segmentMatch = chamberMatch->segmentMatches.begin();segmentMatch != chamberMatch->segmentMatches.end(); segmentMatch++ ) {
	if (!segmentMatch->isMask(MuonSegmentMatch::BestInStationByDR)) continue;
	(*IsSegmentBelongsToTrackByDR)[SavePostion]=segmentMatch->isMask(MuonSegmentMatch::BelongsToTrackByDR);
	(*IsSegmentBelongsToTrackByCleaning)[SavePostion]=segmentMatch->isMask(MuonSegmentMatch::BelongsToTrackByCleaning);
	if (segmentMatch->cscSegmentRef.isNonnull()&&segmentMatch->dtSegmentRef.isNonnull()) {ReportError(12,"A segment belongs to both CSC and DT")}
	else {
	  if (segmentMatch->cscSegmentRef.isNonnull()) (*NumberOfHitsInSegment)[SavePostion]=segmentMatch->cscSegmentRef->specificRecHits().size();
	  if (segmentMatch->dtSegmentRef.isNonnull()) (*NumberOfHitsInSegment)[SavePostion]=segmentMatch->dtSegmentRef->recHits().size();
	}
	(*DXTrackToSegment)[SavePostion]=abs(segmentMatch->x-chamberMatch->x);
	(*DYTrackToSegment)[SavePostion]=abs(segmentMatch->y-chamberMatch->y);
	(*DXErrTrackToSegment)[SavePostion]=sqrt(segmentMatch->xErr*segmentMatch->xErr+chamberMatch->xErr*chamberMatch->xErr);
	(*DYErrTrackToSegment)[SavePostion]=sqrt(segmentMatch->yErr*segmentMatch->yErr+chamberMatch->yErr*chamberMatch->yErr);
	(*DRTrackToSegment)[SavePostion]=sqrt((*DXTrackToSegment)[SavePostion]*(*DXTrackToSegment)[SavePostion]+(*DYTrackToSegment)[SavePostion]*(*DYTrackToSegment)[SavePostion]);
	(*DRErrTrackToSegment)[SavePostion]=sqrt((*DXTrackToSegment)[SavePostion]*(*DXTrackToSegment)[SavePostion]*(*DXErrTrackToSegment)[SavePostion]*(*DXErrTrackToSegment)[SavePostion]+(*DYTrackToSegment)[SavePostion]*(*DYTrackToSegment)[SavePostion]*(*DYErrTrackToSegment)[SavePostion]*(*DYErrTrackToSegment)[SavePostion])/(*DRTrackToSegment)[SavePostion];
	break;
      }
    }
    //isHLTObj
    if (trgEvent.isValid()) 
      for (unsigned int i=0; i<num_HLTsSaveObjs; i++) {
	isHLTObj[i]->push_back(false);
	unsigned int num_thisHLTObjs=HLTObj_pt[i]->size();
	for (unsigned int j=0; j<num_thisHLTObjs; j++) {
	  float dR=(eta->back()-(*HLTObj_eta[i])[j])*(eta->back()-(*HLTObj_eta[i])[j])+(phi->back()-(*HLTObj_phi[i])[j])*(phi->back()-(*HLTObj_phi[i])[j]);
	  if (dR<MaxHLTObjDeviation2) {
	    isHLTObj[i]->back()=true;
	    break;
	  }
	}
      }
    //miscellany
    Vertex_x->push_back(iter->vertex().X());
    Vertex_y->push_back(iter->vertex().Y());
    Vertex_z->push_back(iter->vertex().Z());
    isoR03sumPt->push_back(iter->isolationR03().sumPt);
    isoR03emEt->push_back(iter->isolationR03().emEt);
    isoR03hadEt->push_back(iter->isolationR03().hadEt);
    isoR03hoEt->push_back(iter->isolationR03().hoEt);
    isoR03nJets->push_back(iter->isolationR03().nJets);
    isoR03nTracks->push_back(iter->isolationR03().nTracks);
    isoR05sumPt->push_back(iter->isolationR05().sumPt);
    isoR05emEt->push_back(iter->isolationR05().emEt);
    isoR05hadEt->push_back(iter->isolationR05().hadEt);
    isoR05hoEt->push_back(iter->isolationR05().hoEt);
    isoR05nJets->push_back(iter->isolationR05().nJets);
    isoR05nTracks->push_back(iter->isolationR05().nTracks);
    isoemVetoEt->push_back(iter->isolationR03().emVetoEt);
    isohadVetoEt->push_back(iter->isolationR03().hadVetoEt);
    isohoVetoEt->push_back(iter->isolationR03().hoVetoEt);
    TrkKink->push_back(iter->combinedQuality().trkKink);
    GlbKink->push_back(iter->combinedQuality().glbKink);
    TrkRelChi2->push_back(iter->combinedQuality().trkRelChi2);
    CaloE_emMax->push_back(iter->calEnergy().emMax);
    CaloE_emS9->push_back(iter->calEnergy().emS9);
    CaloE_emS25->push_back(iter->calEnergy().emS25);
    CaloE_hadMax->push_back(iter->calEnergy().hadMax);
    CaloE_hadS9->push_back(iter->calEnergy().hadS9);
    Calo_emPos_R->push_back(iter->calEnergy().ecal_position.R());
    Calo_emPos_eta->push_back(iter->calEnergy().ecal_position.eta());
    Calo_emPos_phi->push_back(iter->calEnergy().ecal_position.phi());
    Calo_hadPos_R->push_back(iter->calEnergy().hcal_position.R());
    Calo_hadPos_eta->push_back(iter->calEnergy().hcal_position.eta());
    Calo_hadPos_phi->push_back(iter->calEnergy().hcal_position.phi());
   //match to RecoTrk
    TrackRef InnerTrack = iter->innerTrack();
    dEdx->push_back(dEdxTrack[InnerTrack].dEdx());
    dEdxError->push_back(dEdxTrack[InnerTrack].dEdxError());
    dEdx_numberOfSaturatedMeasurements->push_back(dEdxTrack[InnerTrack].numberOfSaturatedMeasurements());
    dEdx_numberOfMeasurements->push_back(dEdxTrack[InnerTrack].numberOfMeasurements());
    InnerTrack_nValidTrackerHits->push_back(InnerTrack->hitPattern().numberOfValidTrackerHits());
    InnerTrack_nValidPixelHits->push_back(InnerTrack->hitPattern().numberOfValidPixelHits());
    InnerTrack_nLostPixelHits->push_back(InnerTrack->hitPattern().numberOfLostPixelHits());
    InnerTrack_nLostTrackerHits->push_back(InnerTrack->hitPattern().numberOfLostTrackerHits());
    InnerTrack_chi2->push_back(InnerTrack->chi2());
    InnerTrack_ndof->push_back(InnerTrack->ndof());
    TrackRef GlobalTrack = iter->globalTrack();
    if (GlobalTrack.isNonnull()) {
      GlobalTrack_chi2->push_back(GlobalTrack->chi2());
      GlobalTrack_ndof->push_back(GlobalTrack->ndof());
    }
    else {
           GlobalTrack_chi2->push_back(0);
	   GlobalTrack_ndof->push_back(0);
    } 

#ifdef TrackingParticles
    theSameWithMuon->push_back(-1);
    unsigned int InnerTrack_pos;  DChains->push_back(-2);
    for(InnerTrack_pos=0; InnerTrack_pos<tC.size(); InnerTrack_pos++)  //recotrk loop begin
      if (InnerTrack == TrackRef(trackCollectionH,InnerTrack_pos) ) {
	RefToBase<Track> trk(trackCollectionHV, InnerTrack_pos);
	if(InnerTrack_nValidTrackerHits->back()<minTrackHits||RecoToSimByHits.find(trk) == RecoToSimByHits.end()||RecoToSimByChi2.find(trk) == RecoToSimByChi2.end()) break;
	pair<TrackingParticleRef, double>  BestMatch=RecoToSimByHits[trk].front();
	vector<pair<TrackingParticleRef, double> > TPCByChi2=RecoToSimByChi2[trk];
	vector<pair<TrackingParticleRef, double> >::iterator TPCByChi2_iter=TPCByChi2.begin();
	for (;TPCByChi2_iter!=TPCByChi2.end();TPCByChi2_iter++)
	  if (BestMatch.first==TPCByChi2_iter->first) break;
	if (TPCByChi2_iter==TPCByChi2.end()) break;
	TrackingParticleRef tpr = BestMatch.first;
	//Simulated Tracks
	TrkParticles_pt->push_back(tpr->pt());
	TrkParticles_eta->push_back(tpr->eta());
	TrkParticles_phi->push_back(tpr->phi());
	TrkParticles_charge->push_back(tpr->charge());
	TrkParticles_pdgId->push_back(tpr->pdgId());
	SharedHitsRatio->push_back(BestMatch.second);
	MCMatchChi2->push_back(-1.*TPCByChi2_iter->second);
	//cout<<"pt:"<<track->pt()<<"vs"<<tpr->pt()<<endl<<"eta:"<<track->eta()<<"vs"<<tpr->eta()<<endl<<"phi:"<<track->phi()<<"vs"<<tpr->phi()<<endl<<"chi2:"<<-tp.begin()->second<<endl<<"-------------Next------------------"<<endl;//check it is doing correct things
	
	//Get the decay chain of this track
	MaskOut.clear();
	for (vector<SimTrack>::const_iterator g4Track_iter = tpr->g4Track_begin(); g4Track_iter != tpr->g4Track_end(); ++g4Track_iter )  {//g4Track loop begin
	  SimTrack *thisTrk;
	  trackIdLink(g4Track_iter->trackId())
	  if (find(MaskOut.begin(),MaskOut.end(),thisTrk)==MaskOut.end()) {
	    DChain.clear();SimChains.clear();
	    SimTrackDaughtersTree(thisTrk);
	    SimVertex thisVtx;
	    do {
	      if (!thisTrk->noVertex()) {
		thisVtx=SVC[thisTrk->vertIndex()];
		if (!thisVtx.noParent()) {
		  trackIdLink(thisVtx.parentIndex())
		    //add parent particle to each Chain
		  vector<SimTrack *>::iterator ParToSim_iter=find(ParToSim.begin(),ParToSim.end(),thisTrk);
		  if (ParToSim_iter==ParToSim.end()) {
		    for (vector< vector<Int_t> >::iterator SimChains_iter = SimChains.begin(); SimChains_iter !=  SimChains.end(); ++SimChains_iter )
		      SimChains_iter->insert(SimChains_iter->begin(),IsParInHep->size());
		      RecordSimTrack(thisTrk)
		  }
		  else { Int_t pos=FindSimTrackRecordingPosition(ParToSim_iter-ParToSim.begin());
		    for (vector< vector<Int_t> >::iterator SimChains_iter = SimChains.begin(); SimChains_iter !=  SimChains.end(); ++SimChains_iter )
		      SimChains_iter->insert(SimChains_iter->begin(),pos);
		    }
		}
		else break;
	      }
	      else break;
	    }
	    while(true);
	    //HepMC Particles
	    HepMCChains.clear();
	    if (!thisTrk->noGenpart()) {
	      HepMC::GenEvent::particle_iterator genPar = HepGenEvent.particles_begin();
	      for (int count=1; count<thisTrk->genpartIndex()&&genPar != HepGenEvent.particles_end(); count++ )
		genPar++;
	      if (genPar != HepGenEvent.particles_end()) HepMCParentTree(*genPar);
	      else {ReportError(6,"genpartIndex() Error or HepMC is empty")}
	    }
	    //merge the HepMC and SimTrack Decay Chains
	    for (vector< vector<Int_t> >::iterator SimChains_iter = SimChains.begin(); SimChains_iter !=  SimChains.end(); ++SimChains_iter )
	      for (vector< vector<Int_t> >::iterator HepMCChains_iter = HepMCChains.begin(); HepMCChains_iter !=  HepMCChains.end(); ++HepMCChains_iter ) {
		vector<Int_t> thisChain(HepMCChains_iter->rbegin(),HepMCChains_iter->rend());
		thisChain.insert(thisChain.end(),SimChains_iter->begin(),SimChains_iter->end());
		//see if thisChain is the same with previous muons
		int Muref=-1;
		for (vector<int>::iterator DChain_iter = DChains->begin(); DChain_iter != DChains->end()&&Muref<(int) eta->size()-1; DChain_iter++ ) {
		  if (*DChain_iter==-1) {
		    DChain_iter++;
		    vector<int>::iterator DChain_begin=DChain_iter;
		    for (; DChain_iter != DChains->end()&&*DChain_iter!=-2&&*DChain_iter!=-1; DChain_iter++ ) ;
		    if (DChain_begin!=DChain_iter) 
		      if (IstheSameDChain(thisChain,* new vector<int> (DChain_begin,DChain_iter))) {
			theSameWithMuon->back()=Muref;
			break;
		      }
		    DChain_iter--;
		  }
		  if (*DChain_iter==-2) Muref++;
		}
		DChains->push_back(-1);
		DChains->insert(DChains->end(),thisChain.begin(),thisChain.end());
	      }
	  }
	}//g4Track loop end
	break;
      }//recotrk loop end
    if (TrkParticles_eta->size()<eta->size()) {
      if (InnerTrack_pos==tC.size()) {ReportError(4,"TrkRef Error")}
      if (TrkParticles_eta->size()!=eta->size()-1) {ReportError(5,"MC Match Coding Error")}
      NumMisMatch++;
      TrkParticles_pt->push_back(0);
      TrkParticles_eta->push_back(0);
      TrkParticles_phi->push_back(0);
      TrkParticles_pdgId->push_back(0);
      TrkParticles_charge->push_back(-100.0);
      SharedHitsRatio->push_back(-100.0);
    }
#endif
  }//muon loop end

  //DiMuonInvariantMass
#define Muon_Mass 0.105658367
  TLorentzVector P1,P2;
  unsigned int num_muons=pt->size();
  for (unsigned int Mu1=0; Mu1<num_muons; Mu1++)
    for (unsigned int Mu2=Mu1+1; Mu2<num_muons; Mu2++) {
      P1.SetPtEtaPhiM((*pt)[Mu1],(*eta)[Mu1],(*phi)[Mu1],Muon_Mass);
      P2.SetPtEtaPhiM((*pt)[Mu2],(*eta)[Mu2],(*phi)[Mu2],Muon_Mass);
      DiMuonInvariantMass->push_back((P1+P2).M());
    }
  Muons_Tree->Fill();
  FirstEntry=false;
  return true;
}
#ifdef TrackingParticles
bool DiMuonSelector::IstheSameDChain(vector<int> &ThisChain,vector<int> &AnotherChain) {//whether two chains include each other
  bool ChainIncluded=false;
  vector<int>::iterator ThisChain_Particle = ThisChain.begin(),AnotherChain_Particle = AnotherChain.begin();
  for (; AnotherChain_Particle != AnotherChain.end()&&ThisChain_Particle != ThisChain.end(); AnotherChain_Particle++) {
    if (ChainIncluded&&*ThisChain_Particle!=*AnotherChain_Particle) ChainIncluded=false;
    if (ThisChain.front() == *AnotherChain_Particle) {
      ChainIncluded=true;
      ThisChain_Particle = ThisChain.begin();
    }
    if (ChainIncluded) ThisChain_Particle++;
  }
  if (!ChainIncluded) {
    AnotherChain_Particle = AnotherChain.begin(); ThisChain_Particle = ThisChain.begin();
    for (; AnotherChain_Particle != AnotherChain.end()&&ThisChain_Particle != ThisChain.end(); ThisChain_Particle++) {
      if (ChainIncluded&&*ThisChain_Particle!=*AnotherChain_Particle) ChainIncluded=false;
      if (AnotherChain.front() == *ThisChain_Particle) {
	AnotherChain_Particle = AnotherChain.begin();
	ChainIncluded=true;
      }
      if (ChainIncluded) AnotherChain_Particle++;
    }
  }
  return ChainIncluded;
}

void DiMuonSelector::SimTrackDaughtersTree(SimTrack * thisTrk)
{
  //To avoid duplicate particle saving
  vector<SimTrack *>::iterator ParToSim_iter=find(ParToSim.begin(),ParToSim.end(),thisTrk);
  if (ParToSim_iter==ParToSim.end())
    {
      DChain.push_back(IsParInHep->size());
      RecordSimTrack(thisTrk)
    }
  else DChain.push_back(FindSimTrackRecordingPosition(ParToSim_iter-ParToSim.begin()));
  MaskOut.push_back(thisTrk);
  if (Daughters[thisTrk].size()>0) 
    for (vector<SimTrack *>::iterator Daughter=Daughters[thisTrk].begin();Daughter!=Daughters[thisTrk].end();Daughter++)
      SimTrackDaughtersTree(*Daughter);
  else SimChains.push_back(DChain);
  DChain.pop_back();
}

void DiMuonSelector::HepMCParentTree(HepMC::GenParticle *genPar) {
  HepMC::GenVertex *thisVtx = genPar->production_vertex();
  bool ChainEnd=true;
  if (thisVtx) {
      for (HepMC::GenVertex::particles_in_const_iterator pgenD = thisVtx->particles_in_const_begin(); pgenD != thisVtx->particles_in_const_end(); ++pgenD)
	if ((*pgenD)->pdg_id()!=92)  {//Pythia special code for string, we only care about the particles after hadronization
	  ChainEnd=false;
	  vector<HepMC::GenParticle *>::iterator ParToHep_iter=find(ParToHep.begin(),ParToHep.end(),*pgenD);
	  if (ParToHep_iter==ParToHep.end())
	    {
	      DChain.push_back(IsParInHep->size());
	      RecordHepMC((*pgenD))
	    }
	  else DChain.push_back(FindHepMCRecordingPosition(ParToHep_iter-ParToHep.begin()));
	  HepMCParentTree(*pgenD);
	  DChain.pop_back();
	}
  }
  if (ChainEnd) HepMCChains.push_back(DChain);
}
#endif

void DiMuonSelector::ClearVecs() {
  pt->clear(); eta->clear();  phi->clear();
  chargeMinus->clear();  isGlobalMu->clear();  isTrackerMu->clear();
  Vertex_x->clear();  Vertex_y->clear();  Vertex_z->clear();
  isoR03sumPt->clear();  isoR03emEt->clear();  isoR03hadEt->clear();
  isoR03hoEt->clear();  isoR03nJets->clear();  isoR03nTracks->clear();
  isoR05sumPt->clear();  isoR05emEt->clear();  isoR05hadEt->clear();
  isoR05hoEt->clear();  isoR05nJets->clear();  isoR05nTracks->clear();
  isoemVetoEt->clear();  isohadVetoEt->clear();  isohoVetoEt->clear();
  TrkKink->clear();  GlbKink->clear();  TrkRelChi2->clear();
  CaloE_emMax->clear();  CaloE_emS9->clear();  CaloE_emS25->clear();
  CaloE_hadMax->clear();  CaloE_hadS9->clear();
  Calo_emPos_R->clear();  Calo_emPos_eta->clear();  Calo_emPos_phi->clear();
  Calo_hadPos_R->clear();  Calo_hadPos_eta->clear();  Calo_hadPos_phi->clear();
  DiMuonInvariantMass->clear();
  dEdx->clear();  dEdxError->clear();
  dEdx_numberOfSaturatedMeasurements->clear();
  dEdx_numberOfMeasurements->clear();
  InnerTrack_nValidTrackerHits->clear();
  InnerTrack_nValidPixelHits->clear();  
  InnerTrack_nLostTrackerHits->clear();
  InnerTrack_nLostPixelHits->clear(); 
  InnerTrack_chi2->clear(); InnerTrack_ndof->clear();
  GlobalTrack_chi2->clear(); GlobalTrack_ndof->clear();

  //HLT Objects
  for( unsigned int i =0; i < num_HLTsSaveObjs;i++) {
    HLTObj_pt[i]->clear();
    HLTObj_eta[i]->clear();
    HLTObj_phi[i]->clear();
    isHLTObj[i]->clear();
  }

  TrackDistToChamberEdge->clear();
  TrackDistToChamberEdgeErr->clear();
  DXTrackToSegment->clear();
  DYTrackToSegment->clear();
  DXErrTrackToSegment->clear();
  DYErrTrackToSegment->clear();
  DRTrackToSegment->clear();
  DRErrTrackToSegment->clear();
  IsSegmentBelongsToTrackByDR->clear();
  IsSegmentBelongsToTrackByCleaning->clear();
  NumberOfHitsInSegment->clear();
  StationMask->clear();RequiredStationMask->clear();

  vx->clear(); vy->clear(); vz->clear();
  vxError->clear(); vyError->clear(); vzError->clear();

  for (unsigned int whichcut=0;whichcut<num_Cuts;whichcut++)
    SelectorPassed[whichcut]->clear(); 
  MySelector->clear();

#ifdef HepMCGenParticle
  Gen_pt->clear();  Gen_eta->clear();  Gen_phi->clear(); Gen_pdgId->clear();
  Gen_vx->clear();  Gen_vy->clear();  Gen_vz->clear();  Gen_vt->clear();
  IsParInHep->clear(); ParToHep.clear();
#endif
  
#ifdef TrackingParticles
  TrkParticles_pt->clear();  TrkParticles_eta->clear(); TrkParticles_phi->clear();
  TrkParticles_pdgId->clear(); TrkParticles_charge->clear(); 
  SharedHitsRatio->clear(); MCMatchChi2->clear(); ParToSim.clear(); 
  DChains->clear(); theSameWithMuon->clear();
#endif
}

#define MakeVecBranch(Name,Var,Type) Var=new vector<Type>();Muons_Tree->Branch(Name,&Var)

DiMuonSelector::DiMuonSelector(const edm::ParameterSet& pset) {
//---- Get the input parameters
  FileName = pset.getParameter<string>("FileName");
  Summarization.CrossSection = pset.getParameter<double>("CrossSection");
  maxChamberDist = pset.getUntrackedParameter<double>("maxChamberDist",-3.);
  maxChamberDistPull = pset.getUntrackedParameter<double>("maxChamberDistPull",-3.);
  MuonPtSelection = pset.getUntrackedParameter<double>("MuonPtSelection",15.);
  MuonHighestPtSelection = pset.getUntrackedParameter<double>("MuonHighestPtSelection",40.);
  vector<string> Cuts;
  Cuts = pset.getUntrackedParameter< vector<string> >("StandardMuonCuts",Cuts);//default is empty
  string MuonArbitrationTypeStr;
  MuonArbitrationTypeStr = pset.getUntrackedParameter<string>("MuonArbitrationType","SegmentArbitration");
  MuonArbitrationType=MuonArbitrationTypeFromString(MuonArbitrationTypeStr.c_str());
  minTrackHits = pset.getUntrackedParameter<uint>("minTrackHits",3);
  PrimaryVerticesTag = pset.getUntrackedParameter<string>("PrimaryVertices","offlinePrimaryVerticesWithBS");//"offlinePrimaryVerticesWithBS": Primary vertex reconstructed using the tracks taken from the generalTracks collection, and imposing the offline beam spot as a constraint in the fit of the vertex position. Another possible tag is "offlinePrimaryVertices", which is Primary vertex reconstructed using the tracks taken from the generalTracks collection
  const InputTag tracksTag_default("generalTracks");
  tracksTag = pset.getUntrackedParameter<InputTag>("tracksTag",tracksTag_default);
  dEdxTag = pset.getUntrackedParameter<string>("dEdxTag","dedxHarmonic2");//Other options already available in RECO files are dedxMedian and dedxTruncated40. 
  const InputTag default_TriggerResultsTag("TriggerResults::HLT");
  TriggerResultsTag = pset.getUntrackedParameter<InputTag>("TriggerResultsTag",default_TriggerResultsTag);
  const InputTag default_TriggerEventTag("hltTriggerSummaryAOD","","HLT");
  TriggerEventTag = pset.getUntrackedParameter<InputTag>("triggerEventTag",default_TriggerEventTag);
  HLTFilterNames = pset.getParameter< vector<InputTag> >("hltFilterNames");
  num_HLTsSaveObjs = HLTFilterNames.size();
  if (num_HLTsSaveObjs>maxFilterObjects) throw cms::Exception("You need to increase the constant maxFilterObjects.\n");;

//---- Make a root file and plant a TTree
  file=new TFile(FileName.c_str(),"RECREATE");
  Muons_Tree = new TTree ("MuTrkCand","TrksCandidates") ;
  Summarization_Tree = new TTree ("Summary","Summary") ;
  ErrorMsg_Tree = new TTree ("Errors","Errors") ;
  Muons_Tree->SetCircular(500000);//the max events in a single root file
  ErrorMsg_Tree->SetCircular(500000);
//Build Branches
  //General event information
#ifdef HepMCGenParticle
  HepMCTag = pset.getUntrackedParameter<string>("HepMCTag","generator");
  Muons_Tree->Branch("Event_Info",&Info.RUN,"RUN/l:EVENT:LumiBlock:ORBIT:BunchCrossing:num_ZpMuMuInMC/b:num_ZMuMuInMC:isRealData/O");
#else
  Muons_Tree->Branch("Event_Info",&Info.RUN,"RUN/l:EVENT:LumiBlock:ORBIT:BunchCrossing:isRealData/O");
#endif
  //Muon Information
  MakeVecBranch("pt",pt,Float_t);  MakeVecBranch("eta",eta,Float_t);  MakeVecBranch("phi",phi,Float_t);
  MakeVecBranch("chargeMinus",chargeMinus,Bool_t); MakeVecBranch("isGlobalMu",isGlobalMu,Bool_t);  MakeVecBranch("isTrackerMu",isTrackerMu,Bool_t);
  MakeVecBranch("Vertex_x",Vertex_x,Float_t);  MakeVecBranch("Vertex_y",Vertex_y,Float_t);  MakeVecBranch("Vertex_z",Vertex_z,Float_t);
  MakeVecBranch("isoR03sumPt",isoR03sumPt,Float_t);  MakeVecBranch("isoR03emEt",isoR03emEt,Float_t);  MakeVecBranch("isoR03hadEt",isoR03hadEt,Float_t);
  MakeVecBranch("isoR03hoEt",isoR03hoEt,Float_t);  MakeVecBranch("isoR03nJets",isoR03nJets,Float_t);  MakeVecBranch("isoR03nTracks",isoR03nTracks,Float_t);
  MakeVecBranch("isoR05sumPt",isoR05sumPt,Float_t);  MakeVecBranch("isoR05emEt",isoR05emEt,Float_t);  MakeVecBranch("isoR05hadEt",isoR05hadEt,Float_t);
  MakeVecBranch("isoR05hoEt",isoR05hoEt,Float_t);  MakeVecBranch("isoR05nJets",isoR05nJets,Float_t);  MakeVecBranch("isoR05nTracks",isoR05nTracks,Float_t);
  MakeVecBranch("isoemVetoEt",isoemVetoEt,Float_t);  MakeVecBranch("isohadVetoEt",isohadVetoEt,Float_t);  MakeVecBranch("isohoVetoEt",isohoVetoEt,Float_t);
  MakeVecBranch("TrkKink",TrkKink,Float_t); MakeVecBranch("GlbKink",GlbKink,Float_t); MakeVecBranch("TrkRelChi2",TrkRelChi2,Float_t);
  MakeVecBranch("CaloE_emMax",CaloE_emMax,Float_t);  MakeVecBranch("CaloE_emS9",CaloE_emS9,Float_t);  MakeVecBranch("CaloE_emS25",CaloE_emS25,Float_t);
  MakeVecBranch("CaloE_hadMax",CaloE_hadMax,Float_t);  MakeVecBranch("CaloE_hadS9",CaloE_hadS9,Float_t);
  MakeVecBranch("Calo_emPos_R",Calo_emPos_R,Float_t);  MakeVecBranch("Calo_emPos_eta",Calo_emPos_eta,Float_t);  MakeVecBranch("Calo_emPos_phi",Calo_emPos_phi,Float_t);
  MakeVecBranch("Calo_hadPos_R",Calo_hadPos_R,Float_t);  MakeVecBranch("Calo_hadPos_eta",Calo_hadPos_eta,Float_t);  MakeVecBranch("Calo_hadPos_phi",Calo_hadPos_phi,Float_t);
  MakeVecBranch("DiMuonInvariantMass",DiMuonInvariantMass,Float_t);
  //InnerTrack
  MakeVecBranch("dEdx",dEdx,Float_t);  MakeVecBranch("dEdxError",dEdxError,Float_t);
  MakeVecBranch("dEdx_numberOfSaturatedMeasurements",dEdx_numberOfSaturatedMeasurements,Int_t);
  MakeVecBranch("dEdx_numberOfMeasurements",dEdx_numberOfMeasurements,Int_t);
  MakeVecBranch("Track_nValidTrackerHits",InnerTrack_nValidTrackerHits,UInt_t);
  MakeVecBranch("Track_nValidPixelHits",InnerTrack_nValidPixelHits,UInt_t);
  MakeVecBranch("Track_nLostTrackerHits",InnerTrack_nLostTrackerHits,UInt_t);
  MakeVecBranch("Track_nLostPixelHits",InnerTrack_nLostPixelHits,UInt_t);
  MakeVecBranch("InnerTrack_chi2",InnerTrack_chi2,Float_t);  MakeVecBranch("InnerTrack_ndof",InnerTrack_ndof,UInt_t);
  MakeVecBranch("GlobalTrack_chi2",GlobalTrack_chi2,Float_t);  MakeVecBranch("GlobalTrack_ndof",GlobalTrack_ndof,UInt_t);

  //Muon System - Segment/Chambers
  MakeVecBranch("TrackDistToChamberEdge",TrackDistToChamberEdge,Float_t);
  MakeVecBranch("TrackDistToChamberEdgeErr",TrackDistToChamberEdgeErr,Float_t);
  MakeVecBranch("DXTrackToSegment",DXTrackToSegment,Float_t);
  MakeVecBranch("DYTrackToSegment",DYTrackToSegment,Float_t);
  MakeVecBranch("DXErrTrackToSegment",DXErrTrackToSegment,Float_t);
  MakeVecBranch("DYErrTrackToSegment",DYErrTrackToSegment,Float_t);
  MakeVecBranch("DRTrackToSegment",DRTrackToSegment,Float_t);
  MakeVecBranch("DRErrTrackToSegment",DRErrTrackToSegment,Float_t);
  MakeVecBranch("IsSegmentBelongsToTrackByDR",IsSegmentBelongsToTrackByDR,Bool_t);
  MakeVecBranch("IsSegmentBelongsToTrackByCleaning",IsSegmentBelongsToTrackByCleaning,Bool_t);
  MakeVecBranch("NumberOfHitsInSegment",NumberOfHitsInSegment,UInt_t);
  MakeVecBranch("StationMask",StationMask,UInt_t);MakeVecBranch("RequiredStationMask",RequiredStationMask,UInt_t); 
  //HLT Objects
  for( unsigned int i =0; i < num_HLTsSaveObjs;i++) {
    MakeVecBranch((HLTFilterNames[i].label()+"_pt").c_str(),HLTObj_pt[i],Float_t);
    MakeVecBranch((HLTFilterNames[i].label()+"_eta").c_str(),HLTObj_eta[i],Float_t);
    MakeVecBranch((HLTFilterNames[i].label()+"_phi").c_str(),HLTObj_phi[i],Float_t);
    MakeVecBranch(("is"+HLTFilterNames[i].label()+"Obj").c_str(),isHLTObj[i],Bool_t);
  }
  //PV
  MakeVecBranch("PVx",vx,Float_t);  MakeVecBranch("PVy",vy,Float_t);  MakeVecBranch("PVz",vz,Float_t);
  MakeVecBranch("PVxError",vxError,Float_t);  MakeVecBranch("PVyError",vyError,Float_t);  MakeVecBranch("PVzError",vzError,Float_t);
  //Cuts
  num_Cuts=Cuts.size();
  for (unsigned int whichcut=0;whichcut<num_Cuts;whichcut++) 
    if (muon::selectionTypeFromString(Cuts[whichcut])!=(muon::SelectionType)-1) {
      MakeVecBranch(Cuts[whichcut].c_str(),SelectorPassed[whichcut],Bool_t); 
      OfficialMuonSelectors.push_back(muon::selectionTypeFromString(Cuts[whichcut]));
    }
  num_Cuts=OfficialMuonSelectors.size();
  MakeVecBranch("MySelector",MySelector,Bool_t);
  
#ifdef HepMCGenParticle
  //Generation and Simulation Particles
  MakeVecBranch("Gen_pt",Gen_pt,Float_t);  MakeVecBranch("Gen_eta",Gen_eta,Float_t);  MakeVecBranch("Gen_phi",Gen_phi,Float_t);   MakeVecBranch("Gen_pdgId",Gen_pdgId,Int_t);
  MakeVecBranch("Gen_vx",Gen_vx,Float_t);  MakeVecBranch("Gen_vy",Gen_vy,Float_t);  MakeVecBranch("Gen_vz",Gen_vz,Float_t);  MakeVecBranch("Gen_vt",Gen_vt,Float_t);
  MakeVecBranch("IsParInHep",IsParInHep,Bool_t);
  Muons_Tree->Branch("HepMCVertex",HepMCFourVec,"HepMCVertex[4]/F");
#endif

#ifdef TrackingParticles  
  //Simulated Tracks
  MakeVecBranch("TrkParticles_pt",TrkParticles_pt,Float_t);  MakeVecBranch("TrkParticles_eta",TrkParticles_eta,Float_t);  MakeVecBranch("TrkParticles_phi",TrkParticles_phi,Float_t);
  MakeVecBranch("TrkParticles_pdgId",TrkParticles_pdgId,Int_t); MakeVecBranch("TrkParticles_charge",TrkParticles_charge,Int_t);
  MakeVecBranch("SharedHitsRatio",SharedHitsRatio,Double_t); MakeVecBranch("MCMatchChi2",MCMatchChi2,Double_t); 
  MakeVecBranch("DChains",DChains,Int_t); MakeVecBranch("theSameWithMuon",theSameWithMuon,Int_t);
  NumMisMatch=0;
#endif

  //Summary
  ErrorMsg_Tree->Branch("ErrorMsg",&Error.ErrorCode,"ErrorCode/b:RunNum/l:EventNum");
  Summarization_Tree->Branch("Summarization",&Summarization.Total_Events,"Total_Events/l:Total_TrackerMuons:Total_GlobalMuon:Total_GlobalnotTrackerMuon:CrossSection/D");
  Summarization.Total_Events=0;Summarization.Total_TrackerMuons=0;Summarization.Total_GlobalMuon=0;Summarization.Total_GlobalnotTrackerMuon=0;
  FirstEntry=true;HLTSize=0;
}

Muon::ArbitrationType DiMuonSelector::MuonArbitrationTypeFromString( const std::string &label ) {
  struct MuonArbitrationTypeStringToEnum { const char *label; reco::Muon::ArbitrationType value; };
  static MuonArbitrationTypeStringToEnum MuonArbitrationTypeStringToEnumMap[] = {
    {"NoArbitration",Muon::NoArbitration},
    {"SegmentArbitration",Muon::SegmentArbitration},
    {"SegmentAndTrackArbitration",Muon::SegmentAndTrackArbitration},
    {"SegmentAndTrackArbitrationCleaned",Muon::SegmentAndTrackArbitrationCleaned},
    { 0, (Muon::ArbitrationType)-1 }
  };
  Muon::ArbitrationType value = (Muon::ArbitrationType)-1;
  bool found = false;
  for(int i = 0; MuonArbitrationTypeStringToEnumMap[i].label && (! found); ++i)
    if (! strcmp(label.c_str(), MuonArbitrationTypeStringToEnumMap[i].label)) {
      found = true;
      value = MuonArbitrationTypeStringToEnumMap[i].value;
    }
  
  // in case of unrecognized selection type
  if (! found) throw cms::Exception("MuonSelectorError") << label << " is not a recognized SelectionType";
  return value;
}

DiMuonSelector::~DiMuonSelector()
{
  //printf("MisRate:%f%%\n",NumMisMatch/(float) (Summarization.Total_TrackerMuons+Summarization.Total_GlobalnotTrackerMuon)*100);
  //Summarizations
  Summarization_Tree->Fill();
  file->Write();
  file->Close();
}

// Release all memory when clear vectors - Supposed to fix the memory leak problem - swap with empty vectors - does not work
/*
void DiMuonSelector::ClearVecs() {
  vector<Float_t>().swap(*pt); vector<Float_t>().swap(*eta); vector<Float_t>().swap(*phi);
  vector<Bool_t>().swap(*chargeMinus); vector<Bool_t>().swap(*isGlobalMu); vector<Bool_t>().swap(*isTrackerMu);
  vector<Float_t>().swap(*Vertex_x); vector<Float_t>().swap(* Vertex_y); vector<Float_t>().swap(*Vertex_z);
  vector<Float_t>().swap(*isoR03sumPt); vector<Float_t>().swap(*isoR03emEt); vector<Float_t>().swap(*isoR03hadEt);
  vector<Float_t>().swap(*isoR03hoEt); vector<Float_t>().swap(*isoR03nJets); vector<Float_t>().swap(*isoR03nTracks);
  vector<Float_t>().swap(*isoR05sumPt); vector<Float_t>().swap(*isoR05emEt); vector<Float_t>().swap(*isoR05hadEt);
  vector<Float_t>().swap(*isoR05hoEt); vector<Float_t>().swap(*isoR05nJets); vector<Float_t>().swap(*isoR05nTracks);
  vector<Float_t>().swap(*isoemVetoEt); vector<Float_t>().swap(*isohadVetoEt); vector<Float_t>().swap(*isohoVetoEt);
  vector<Float_t>().swap(*TrkKink); vector<Float_t>().swap(*GlbKink); vector<Float_t>().swap(*TrkRelChi2);
  vector<Float_t>().swap(*CaloE_emMax); vector<Float_t>().swap(*CaloE_emS9); vector<Float_t>().swap(*CaloE_emS25);
  vector<Float_t>().swap(*CaloE_hadMax); vector<Float_t>().swap(*CaloE_hadS9);
  vector<Float_t>().swap(*Calo_emPos_R); vector<Float_t>().swap(*Calo_emPos_eta); vector<Float_t>().swap(*Calo_emPos_phi);
  vector<Float_t>().swap(*Calo_hadPos_R); vector<Float_t>().swap(*Calo_hadPos_eta); vector<Float_t>().swap(*Calo_hadPos_phi);
  vector<Float_t>().swap(*DiMuonInvariantMass);
  vector<Float_t>().swap(*dEdx); vector<Float_t>().swap(*dEdxError);
  vector<Int_t>().swap(*dEdx_numberOfSaturatedMeasurements);
  vector<Int_t>().swap(*dEdx_numberOfMeasurements);
  vector<UInt_t>().swap(*InnerTrack_nValidTrackerHits);
  vector<UInt_t>().swap(*InnerTrack_nValidPixelHits); 
  vector<UInt_t>().swap(*InnerTrack_nLostTrackerHits);
  vector<UInt_t>().swap(*InnerTrack_nLostPixelHits);
  vector<Float_t>().swap(*InnerTrack_chi2); vector<UInt_t>().swap(*InnerTrack_ndof);
  vector<Float_t>().swap(*GlobalTrack_chi2); vector<UInt_t>().swap(*GlobalTrack_ndof);

  //HLT Objects
  for( unsigned int i =0; i < num_HLTsSaveObjs;i++) {
    vector<Float_t>().swap(*HLTObj_pt[i]);
    vector<Float_t>().swap(*HLTObj_eta[i]);
    vector<Float_t>().swap(*HLTObj_phi[i]);
    vector<Bool_t>().swap(*isHLTObj[i]);
  }

  vector<Float_t>().swap(*TrackDistToChamberEdge);
  vector<Float_t>().swap(*TrackDistToChamberEdgeErr);
  vector<Float_t>().swap(*DXTrackToSegment);
  vector<Float_t>().swap(*DYTrackToSegment);
  vector<Float_t>().swap(*DXErrTrackToSegment);
  vector<Float_t>().swap(*DYErrTrackToSegment);
  vector<Float_t>().swap(*DRTrackToSegment);
  vector<Float_t>().swap(*DRErrTrackToSegment);
  vector<Bool_t>().swap(*IsSegmentBelongsToTrackByDR);
  vector<Bool_t>().swap(*IsSegmentBelongsToTrackByCleaning);
  vector<UInt_t>().swap(*NumberOfHitsInSegment);
  vector<UInt_t>().swap(*StationMask);vector<UInt_t>().swap(*RequiredStationMask);

  vector<Float_t>().swap(*vx); vector<Float_t>().swap(*vy); vector<Float_t>().swap(*vz);
  vector<Float_t>().swap(*vxError); vector<Float_t>().swap(*vyError); vector<Float_t>().swap(*vzError);

  for (unsigned int whichcut=0;whichcut<num_Cuts;whichcut++)
    vector<Bool_t>().swap(*SelectorPassed[whichcut]);
  vector<Bool_t>().swap(*MySelector);

#ifdef HepMCGenParticle
  vector<Float_t>().swap(*Gen_pt); vector<Float_t>().swap(*Gen_eta); vector<Float_t>().swap(*Gen_phi); vector<Int_t>().swap(*Gen_pdgId);
  vector<Float_t>().swap(*Gen_vx); vector<Float_t>().swap(*Gen_vy); vector<Float_t>().swap(*Gen_vz); vector<Float_t>().swap(*Gen_vt);
  vector<Bool_t>().swap(*IsParInHep);
#endif
  
#ifdef TrackingParticles
  vector<Float_t>().swap(*TrkParticles_pt); vector<Float_t>().swap(*TrkParticles_eta); vector<Float_t>().swap(*TrkParticles_phi);
  vector<Int_t>().swap(*TrkParticles_pdgId); vector<Int_t>().swap(*TrkParticles_charge); 
  vector<Double_t>().swap(*SharedHitsRatio); vector<Double_t>().swap(*MCMatchChi2);
  vector<SimTrack *>().swap(ParToSim); vector<HepMC::GenParticle *>().swap(*ParToHep);
  vector<Int_t>().swap(*DChains); vector<Int_t>().swap(*theSameWithMuon);
#endif
}*/

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonSelector);
