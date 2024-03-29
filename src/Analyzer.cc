// -*- C++ -*-

// Package:    Analyzer
// Clas:      Analyzer
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
// $Id: Analyzer.cc,v 1.10 2010/10/19 00:39:38 zhangjin Exp $
//

#include "../interface/Analyzer.h"
#define Muon_Mass 0.105658367

#define trackIdLink(TrkID)  vector<SimTrack>::iterator Trk_iter = STC.begin(); \
  for (; Trk_iter != STC.end(); ++Trk_iter )	\
    if ((int) Trk_iter->trackId() == (int) TrkID) break;		\
  if (Trk_iter==STC.end()) LogWarning("RefNULL")<<"parentIndex/trackId Error"<<endl; \
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

#ifdef TrackingParticles
#define RecordHepMC(GenParticle) ParToHep.push_back(GenParticle);	\
  IsParInHep->push_back(true);						\
  GenSimMomentum(sqrt(GenParticle->momentum().px()*GenParticle->momentum().px()+GenParticle->momentum().py()*GenParticle->momentum().py()),GenParticle->momentum().eta(),GenParticle->momentum().phi(),GenParticle->pdg_id()); \
  HepMC::GenVertex *thisVtx=GenParticle->production_vertex();		\
  if (thisVtx) {GenSimVertex(thisVtx->position().x()/10.,thisVtx->position().y()/10.,thisVtx->position().z()/10.,thisVtx->position().t()/299.792458);} \
  else {GenSimVertex(0,0,0,0);}
#else
#define RecordHepMC(GenParticle) ParToHep.push_back(GenParticle);	\
  GenSimMomentum(sqrt(GenParticle->momentum().px()*GenParticle->momentum().px()+GenParticle->momentum().py()*GenParticle->momentum().py()),GenParticle->momentum().eta(),GenParticle->momentum().phi(),GenParticle->pdg_id()); \
  HepMC::GenVertex *thisVtx=GenParticle->production_vertex();		\
  if (thisVtx) {GenSimVertex(thisVtx->position().x()/10.,thisVtx->position().y()/10.,thisVtx->position().z()/10.,thisVtx->position().t()/299.792458);} \
  else {GenSimVertex(0,0,0,0);}
#endif

// ------------ method called on each new Event  ------------

void MuESelector::analyze(const edm::Event& event, const edm::EventSetup& iSetup) {
  Bool_t SaveEvent=true;
//########## preselect muons according to optimized pt() #####################
  Handle<reco::MuonCollection> Muon;
  event.getByLabel("muons", Muon);
  typedef pair<reco::MuonCollection::const_iterator,double> MuonPreSelection;
  vector<MuonPreSelection> MuonQueue;
  reco::MuonCollection const & muons = *Muon;
  Handle <reco::TrackToTrackMap> tevMapH1,tevMapH2,tevMapH3;
  event.getByLabel("tevMuons", "default", tevMapH1);
  const reco::TrackToTrackMap tevMap1 = *(tevMapH1.product());
  event.getByLabel("tevMuons", "firstHit", tevMapH2);
  const reco::TrackToTrackMap tevMap2 = *(tevMapH2.product());
  event.getByLabel("tevMuons", "picky", tevMapH3);
  const reco::TrackToTrackMap tevMap3 = *(tevMapH3.product());
  if (Muon.failedToGet()) SaveEvent=false;
  else {
    for ( reco::MuonCollection::const_iterator iter = muons.begin(); iter != muons.end() ; ++iter) {
      Float_t pt;
      if ( iter->globalTrack().isNonnull() ) pt=muon::tevOptimized(*iter,tevMap1,tevMap2,tevMap3)->pt();
      else pt=iter->pt();
      if (pt>MuonPtCut) MuonQueue.push_back(make_pair(iter,pt));
    }
    if (MuonQueue.size()<2) SaveEvent=false;
  }
//########## Generation Information #####################
//IsMC? based on event content
  if (FirstEvent) {
    Handle<edm::HepMCProduct> HepMCH;
    event.getByLabel(HepMCTag, HepMCH);
    if (HepMCH.failedToGet()){
      LogInfo("DataFormat")<<"It is real data. HepMC collection is disabled"<<endl;
      if ( !event.isRealData() ) LogWarning("DataLost")<<"Corrupted MC samples! HepMC collection is missing"<<endl;
      IsMC=false;
    }
    FirstEvent=false;
  }
//Save all generated Z/gamma*, Z' events for efficiency calculation
  HepMC::GenEvent *HepGenEvent=NULL;
  if (IsMC) {
    ClearVecs_HepMC();
    Handle<edm::HepMCProduct> HepMCH;
    event.getByLabel(HepMCTag, HepMCH);
    HepGenEvent=new HepMC::GenEvent(*(HepMCH->GetEvent()));
    for (HepMC::GenEvent::particle_iterator GenParticle_iter = HepGenEvent->particles_begin(); GenParticle_iter != HepGenEvent->particles_end();GenParticle_iter++) 
      if (abs((*GenParticle_iter)->pdg_id())==13) { 
	HepMC::GenVertex *thisVtx=(*GenParticle_iter)->production_vertex();
	if (thisVtx) for (HepMC::GenVertex::particles_in_const_iterator pgenD = thisVtx->particles_in_const_begin(); pgenD != thisVtx->particles_in_const_end(); ++pgenD)
	  if ((*pgenD)->pdg_id()==32||(*pgenD)->pdg_id()==23) {//if the parent particle is Z or Z'
	    vector<HepMC::GenParticle *>::iterator ParToHep_iter=find(ParToHep.begin(),ParToHep.end(),*pgenD);
	    if (ParToHep_iter==ParToHep.end()) {RecordHepMC((*pgenD))}
	    RecordHepMC((*GenParticle_iter))
	      SaveEvent=true;
	  }
      }//search for Z'->MuMu in MonteCarlo
    for (unsigned int i=0;i<4;i++)
      HepMCFourVec[i]=999999.;
    HepMC::GenEvent::particle_iterator GenParticle_iter = HepGenEvent->particles_begin();
    for (;GenParticle_iter != HepGenEvent->particles_end();GenParticle_iter++) { 
      HepMC::GenVertex *thisVtx=(*GenParticle_iter)->production_vertex();
      if (thisVtx) {
	HepMCFourVec[0]=thisVtx->position().x()/10.;
	HepMCFourVec[1]=thisVtx->position().y()/10.;
	HepMCFourVec[2]=thisVtx->position().z()/10.;
      HepMCFourVec[3]=thisVtx->position().t()/299.792458;
      break;
      }
    }
    if (GenParticle_iter == HepGenEvent->particles_end()) LogWarning("DataLost")<<"No HepMC(Core) Vertex Information"<<endl;
    Handle<GenEventInfoProduct> hEvtInfo;
    event.getByLabel("generator", hEvtInfo);
    GenEventWeight = hEvtInfo->weight();
  }
  else GenEventWeight = 1.;

  if (!SaveEvent) return;
//############### standard filters ############################
  isHLTTriggerred=HLTFilter->TriggerFilter(event,iSetup);
  isGoodVertex=VertexFilter->filter(event,iSetup);
  isNoScrapping=ScrapingFilter->filter(event);
//############### sort Muon according to Tuned pt ############
  ClearVecs_RECO();
  if (MuonQueue.size()>1)
    for (vector<MuonPreSelection>::iterator iter_queue=MuonQueue.begin();iter_queue!=MuonQueue.end()-1; ++iter_queue )
      for (vector<MuonPreSelection>::iterator iter_queue2=iter_queue+1;iter_queue2!=MuonQueue.end(); ++iter_queue2 )
	if (iter_queue->second<iter_queue2->second) swap(*iter_queue,*iter_queue2);

//############### event information ##########################
  Info.RUN   = event.id ().run ();
  Info.EVENT = event.id ().event();
  Info.LS    = event.luminosityBlock ();
  Info.ORBIT = event.orbitNumber ();
  Info.BX = event.bunchCrossing ();
  isRealData = event.isRealData();

  //Beamspot
  Handle<reco::BeamSpot> beamSpotHandle;
  if (!event.getByLabel(InputTag("offlineBeamSpot"), beamSpotHandle)) LogWarning("DataLost")<<"BeamSpot Information is not found."<<endl;

  //Primary Vertex
  Handle<reco::VertexCollection> recVtxs;
  event.getByLabel(PrimaryVerticesTag.c_str(),recVtxs);
  if (recVtxs.isValid())
    for(reco::VertexCollection::const_iterator v=recVtxs->begin(); v!=recVtxs->end(); ++v) {
      vx->push_back(v->x());	   vxError->push_back(v->xError());
      vy->push_back(v->y());	   vyError->push_back(v->yError());
      vz->push_back(v->z());	   vzError->push_back(v->zError());
    }
  else LogWarning("InvalidTag")<<PrimaryVerticesTag.c_str()<<" information is not valid."<<endl;
  
  // this is needed by the IPTools methods from the tracking group
  ESHandle<TransientTrackBuilder> trackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);
  //Pileup
  if (IsMC) {
    edm::InputTag PileupSrc_("addPileupInfo");
    Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    if ( !event.getByLabel(PileupSrc_, PupInfo) ) LogWarning("DataLost")<<"Pileup Information is not found."<<endl;
    else {
      numberOfPUVertices=0;
      for(std::vector<PileupSummaryInfo>::const_iterator PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI)
	numberOfPUVertices+=PVI->getPU_NumInteractions();
    }
  }

  //HLT
  Handle<TriggerResults> hltTriggerResults;
  event.getByLabel(TriggerResultsTag,hltTriggerResults);
  if (hltTriggerResults.isValid())
    for (UInt_t i=0; i < HLTSize; i++)
      HLTacceptance->push_back(hltTriggerResults->accept(i));
  else LogWarning("InvalidTag")<<"TriggerResultsTag: "<<TriggerResultsTag.label()<<"is not valid"<<endl;

  //HLT Objects
  Handle<trigger::TriggerEvent> trgEvent;
  event.getByLabel(TriggerEventTag,trgEvent);
  if (trgEvent.isValid()) {
    const trigger::TriggerObjectCollection& TOC = trgEvent->getObjects();
    int Total_Filters = trgEvent->sizeFilters();
    for( unsigned int i =0; i < num_HLTsSaveObjs;i++) 
      if (HLTFilterNamesAcceptenceIndex[i]>=0) {
	trigger::size_type pos_Filter = trgEvent->filterIndex(HLT_ModuleNames[i]);
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
	else if ((*HLTacceptance)[HLTFilterNamesAcceptenceIndex[i]])
	  LogWarning("DataLost")<<(*HLTNamesSet)[HLTFilterNamesAcceptenceIndex[i]]<<" passed but "<<HLT_ModuleNames[i].label()<<" is not found."<<endl;
      }
  }
  else LogWarning("InvalidTag")<<"TriggerEventTag "<<TriggerEventTag.label()<<" is not valid."<<endl;

//################# Tracks Information #####################
#ifdef TrackingParticles
  //Reco Tracks
  Handle<reco::TrackCollection> trackCollectionH;
  event.getByLabel(tracksTag,trackCollectionH);
  reco::TrackCollection tC = *trackCollectionH.product();
  Handle< View<Track> > trackCollectionHV;
  event.getByLabel(tracksTag,trackCollectionHV);

  Handle<SimVertexContainer> SVCollectionH;
  event.getByLabel("g4SimHits", SVCollectionH);
  SVC = *SVCollectionH.product();//Simulated Vertices: the vertexId() is just the position
  
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
  //DEDX
  Handle< ValueMap<reco::DeDxData> > dEdxTrackHandle;
  event.getByLabel(dEdxTag, dEdxTrackHandle);
  const ValueMap<reco::DeDxData> dEdxTrack = *(dEdxTrackHandle.product());
//################# Main Muon LOOP #########################
  for (vector<MuonPreSelection>::iterator iter_queue=MuonQueue.begin();iter_queue!=MuonQueue.end(); ++iter_queue ) {//muon loop begin
    reco::MuonCollection::const_iterator iter=iter_queue->first;
    //muon basic information (pt,eta,phi,charge)
    //Global Muon uses Cocktail track selection, Track Muon uses tracker track
    TrackRef InnerTrack = iter->innerTrack();
    TrackRef GlobalTrack = iter->globalTrack();
    if ( InnerTrack.isNull() ) continue;
    if ( iter->isGlobalMuon() ) {
      GlobalTrack = muon::tevOptimized(GlobalTrack, InnerTrack, tevMap1, tevMap2, tevMap3);//change the global track to the optimized track (TuneP algo)
      pt->push_back(GlobalTrack->pt());
      eta->push_back(GlobalTrack->eta());
      phi->push_back(GlobalTrack->phi());
      isTrackerMu->push_back(iter->isTrackerMuon());
      isGlobalMu->push_back(true);
    }
    else {
      if ( iter->isTrackerMuon() ) {
	pt->push_back(iter->pt());
	eta->push_back(iter->eta());
	phi->push_back(iter->phi());
	isTrackerMu->push_back(true);
	isGlobalMu->push_back(false);
      }
      else continue;
    }
    if (GlobalTrack.isNonnull()) {
      numberOfValidMuonHits->push_back( GlobalTrack->hitPattern().numberOfValidMuonHits() );
      GlobalTrack_chi2->push_back(GlobalTrack->chi2());
      GlobalTrack_ndof->push_back(GlobalTrack->ndof());
    }
    else {
      numberOfValidMuonHits->push_back(0);
      GlobalTrack_chi2->push_back(0);
      GlobalTrack_ndof->push_back(0);
    } 
    if (iter->charge()==-1) chargeMinus->push_back(true);
    else chargeMinus->push_back(false);
    //Muon Selectors
    for (unsigned int whichcut=0;whichcut<num_Cuts;whichcut++)
      SelectorPassed[whichcut]->push_back(muon::isGoodMuon(*iter,OfficialMuonSelectors[whichcut]));
    MySelector->push_back(muon::isGoodMuon(*iter,muon::TMLastStation,1,3,3,3,3,maxChamberDist,maxChamberDistPull,MuonArbitrationType,false,true));
    //Station and Segment Matches
    StationMask->push_back(iter->stationMask(MuonArbitrationType));
    RequiredStationMask->push_back(muon::RequiredStationMask(*iter,maxChamberDist,maxChamberDistPull,MuonArbitrationType));
    //isHLTObj
    for (unsigned int i=0; i<num_HLTsSaveObjs; i++)
      isHLTObj[i]->push_back(false);
    if ( trgEvent.isValid() )
      for (unsigned int i=0; i<num_HLTsSaveObjs; i++) {
	unsigned int num_thisHLTObjs=HLTObj_pt[i]->size();
	for (unsigned int j=0; j<num_thisHLTObjs; j++) {
	  float dR=(eta->back()-(*HLTObj_eta[i])[j])*(eta->back()-(*HLTObj_eta[i])[j])+(phi->back()-(*HLTObj_phi[i])[j])*(phi->back()-(*HLTObj_phi[i])[j]);
	  if (dR<MaxHLTObjDeviation2) {
	    isHLTObj[i]->back()=true;
	    break;
	  }
	}
      }
    //Muon Vertex
    Vertex_x->push_back(iter->vertex().X());
    Vertex_y->push_back(iter->vertex().Y());
    Vertex_z->push_back(iter->vertex().Z());
    //numberofmatches
    numberOfMatchedSegments->push_back( iter->numberOfMatches() );
    numberOfMatchedStations->push_back( iter->numberOfMatchedStations() );
    //Isolation
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
    //Kinks
    TrkKink->push_back(iter->combinedQuality().trkKink);
    GlbKink->push_back(iter->combinedQuality().glbKink);
    //Combined track quality
    TrkRelChi2->push_back(iter->combinedQuality().trkRelChi2);

   //match to RecoTrk
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
    //DXY using beamspot
    DXYwtBS->push_back( InnerTrack->dxy(beamSpotHandle->position()) );
    //DXY using first PV
    std::pair<bool,Measurement1D> result = IPTools::absoluteTransverseImpactParameter(trackBuilder->build(InnerTrack),recVtxs->at(0));
    DXYwtPV->push_back( result.second.value() );
    DXYErrwtPV->push_back( result.second.error() );
    
//chambers being considered
    vector<const MuonChamberMatch *> Chambers;
    for( vector<MuonChamberMatch>::const_iterator chamberMatch = iter->matches().begin();chamberMatch != iter->matches().end(); chamberMatch++ )
      if ( chamberMatch->detector()!=MuonSubdetId::RPC ) Chambers.push_back( &(*chamberMatch) );
    //align save position from DT to CSC, from - to +, from inner to outter,from smaller chamber # to big
    if (Chambers.size()>1)
      for( vector<const MuonChamberMatch *>::iterator chamberMatch1 = Chambers.begin();chamberMatch1 != Chambers.end()-1; chamberMatch1++ ) {
	Bool_t IsCSC1=(*chamberMatch1)->detector()==MuonSubdetId::CSC;
	Int_t wheelstation1=0,stationring1=0,sectorchamber1=0;
	if ( IsCSC1 ) {
	  const CSCDetId ChamberID( (*chamberMatch1)->id.rawId() );
	  wheelstation1=ChamberID.endcap()==1?-ChamberID.station():ChamberID.station();
	  stationring1=ChamberID.ring()%4;
	  sectorchamber1=ChamberID.chamber();
	}
	else {
	  const DTChamberId ChamberID( (*chamberMatch1)->id.rawId() );
	  wheelstation1=ChamberID.station();
	  stationring1=ChamberID.wheel();
	  sectorchamber1=ChamberID.sector();
	}
	for( vector<const MuonChamberMatch *>::iterator chamberMatch2 = chamberMatch1+1;chamberMatch2 != Chambers.end(); chamberMatch2++ ) {
	  Bool_t IsCSC2=(*chamberMatch2)->detector()==MuonSubdetId::CSC;
	  if ( IsCSC2&&! IsCSC1 ) continue;
	  Int_t wheelstation2=0,stationring2=0,sectorchamber2=0;
	  if ( IsCSC2 ) {
	    const CSCDetId ChamberID( (*chamberMatch2)->id.rawId() );
	    wheelstation2=ChamberID.endcap()==1?-ChamberID.station():ChamberID.station();
	    stationring2=ChamberID.ring()%4;
	    sectorchamber2=ChamberID.chamber();
	  }
	  else {
	    const DTChamberId ChamberID( (*chamberMatch2)->id.rawId() );
	    wheelstation2=ChamberID.station();
	    stationring2=ChamberID.wheel();
	    sectorchamber2=ChamberID.sector();
	  }
	  if ( IsCSC2==IsCSC1 ) {
	    if ( wheelstation2<wheelstation1 ) continue;
	    if ( wheelstation2==wheelstation1 ) {
	      if ( stationring2>stationring1 ) continue;
	      if ( ( stationring2==stationring1 )&&( sectorchamber2>sectorchamber1 ) ) continue;
	    }
	  }
	  swap(*chamberMatch1,*chamberMatch2);
	}
      }
    //save the chambers and the segments in the chambers information 
    for( vector<const MuonChamberMatch *>::iterator chamberMatch = Chambers.begin();chamberMatch != Chambers.end(); chamberMatch++ ) {//chamberloop begin
      std::vector<MuonSegmentMatch>::const_iterator segmentMatch = (*chamberMatch)->segmentMatches.begin();
      for( ;segmentMatch != (*chamberMatch)->segmentMatches.end(); segmentMatch++ ) 
	if (segmentMatch->isMask(MuonSegmentMatch::BestInStationByDR)) break;
      if ( segmentMatch == (*chamberMatch)->segmentMatches.end() &&  (*chamberMatch)->dist()==999999 ) continue;
      Bool_t IsCSC=(*chamberMatch)->detector()==MuonSubdetId::CSC;
      //Chamber Information
      IsCSCChamber->push_back(IsCSC);
      if ( IsCSC ) {
	const CSCDetId ChamberID( (*chamberMatch)->id.rawId() );
	StationRing->push_back( ChamberID.endcap()==1?ChamberID.station()*10+ChamberID.ring():-ChamberID.station()*10-ChamberID.ring() );
	SectorChamber->push_back( ChamberID.chamber() );
      }
      else {
	const DTChamberId ChamberID( (*chamberMatch)->id.rawId() );
	StationRing->push_back( ChamberID.wheel()>0?ChamberID.wheel()*10+ChamberID.station():ChamberID.wheel()*10-ChamberID.station() );
	SectorChamber->push_back( ChamberID.sector() );
      }
      MuonIndex->push_back( eta->size()-1 );
      TrackDistToChamberEdge->push_back( (*chamberMatch)->dist() );
      TrackDistToChamberEdgeErr->push_back( (*chamberMatch)->distErr() );
      
      XTrack->push_back( (*chamberMatch)->x );
      YTrack->push_back( (*chamberMatch)->y );
      XErrTrack->push_back( (*chamberMatch)->xErr );
      YErrTrack->push_back( (*chamberMatch)->yErr );
      XSegment->push_back(999999);
      YSegment->push_back(999999);
      XErrSegment->push_back(999999);
      YErrSegment->push_back(999999);
      NumberOfHitsInSegment->push_back(0);
      IsSegmentOwnedExclusively->push_back( false );
      IsSegmentBelongsToTrackByDR->push_back( false );
      IsSegmentBelongsToTrackByCleaning->push_back( false );
      if ( segmentMatch != (*chamberMatch)->segmentMatches.end() ) {
	if (IsCSC) 
	  if (segmentMatch->cscSegmentRef.isNonnull()) NumberOfHitsInSegment->back() = segmentMatch->cscSegmentRef->nRecHits();
	  else LogWarning("RefNull")<<"ME "<<StationRing->back()<<", Segment Ref not found"<<endl;
	else
	  if (segmentMatch->dtSegmentRef.isNonnull()) NumberOfHitsInSegment->back() = segmentMatch->dtSegmentRef->recHits().size();
	  else LogWarning("RefNull")<<"MB "<<StationRing->back()<<", Segment Ref not found"<<endl;	   
	XSegment->back() = segmentMatch->x;
	YSegment->back() = segmentMatch->y;
	XErrSegment->back() = segmentMatch->xErr;
	YErrSegment->back() = segmentMatch->yErr;
	UInt_t CurrentChamber=XTrack->size()-1;
	IsSegmentOwnedExclusively->back()=true;
	for (UInt_t chamber=0; chamber<CurrentChamber&&MuonIndex->at(chamber)!=MuonIndex->back(); chamber++)
	  if ( IsTheSameSegment(chamber,CurrentChamber) ) {
	    IsSegmentOwnedExclusively->back()=false;
	    IsSegmentOwnedExclusively->at(chamber)=false;
	  }
	IsSegmentBelongsToTrackByDR->back() = segmentMatch->isMask(MuonSegmentMatch::BelongsToTrackByDR);
	IsSegmentBelongsToTrackByCleaning->back() = segmentMatch->isMask(MuonSegmentMatch::BelongsToTrackByCleaning);
      }
    }//chamber loop end
#ifdef TrackingParticles
    //Match to Tracking Particle Collection
    SavedTP.clear();
    TrkParticles_pt->push_back(0);
    TrkParticles_eta->push_back(0);
    TrkParticles_phi->push_back(0);
    TrkParticles_pdgId->push_back(0);
    TrkParticles_charge->push_back(-100.0);
    SharedHitsRatio->push_back(-100.0);
    MCMatchChi2->push_back(-100.0);
    TTTruthDChains->push_back(-2);
    TTTruthMuType->push_back(0);
    SegTruthDChains->push_back(-2);
    SegTruthMuType->push_back(0);
    NumMisMatch++;
    ChamberSimHits.clear();
    RefToBase<Track> trk(InnerTrack);
    if (InnerTrack_nValidTrackerHits->back()>=minTrackHits&&RecoToSimByHits.find(trk) != RecoToSimByHits.end()&&RecoToSimByChi2.find(trk) != RecoToSimByChi2.end()) {
      pair<TrackingParticleRef, double>  BestMatch=RecoToSimByHits[trk].front();
      vector<pair<TrackingParticleRef, double> > TPCByChi2=RecoToSimByChi2[trk];
      vector<pair<TrackingParticleRef, double> >::iterator TPCByChi2_iter=TPCByChi2.begin();
      for (;TPCByChi2_iter!=TPCByChi2.end();TPCByChi2_iter++)
	if (BestMatch.first==TPCByChi2_iter->first) break;
      if (TPCByChi2_iter!=TPCByChi2.end()) {
	TrackingParticleRef tpr = BestMatch.first;
	//Simulated Tracks
	TrkParticles_pt->back()=tpr->pt();
	TrkParticles_eta->back()=tpr->eta();
	TrkParticles_phi->back()=tpr->phi();
	TrkParticles_charge->back()=tpr->charge();
	TrkParticles_pdgId->back()=tpr->pdgId();
	SharedHitsRatio->back()=BestMatch.second;
	MCMatchChi2->back()=-1.*TPCByChi2_iter->second;
	NumMisMatch--;
	//Get the decay chain of this track
	vector <TheMuonType> type;
	GetDecayChains(tpr,TTTruthDChains,type,HepGenEvent);
	Long64_t TypeRecord=0;
	for (vector <TheMuonType>::iterator type_iter=type.begin();type_iter!=type.end();type_iter++)
	  if (*type_iter) TypeRecord=TypeRecord*100+Long64_t(*type_iter);
	TTTruthMuType->back()=TypeRecord;
      }
    }
#endif
  }//muon loop end
  //DiMuonInvariantMass
  TLorentzVector P1,P2,Pc;
  unsigned int num_muons=pt->size();
  for (unsigned int Mu1=0; Mu1<num_muons; Mu1++)
    for (unsigned int Mu2=Mu1+1; Mu2<num_muons; Mu2++) {
      P1.SetPtEtaPhiM((*pt)[Mu1],(*eta)[Mu1],(*phi)[Mu1],Muon_Mass);
      P2.SetPtEtaPhiM((*pt)[Mu2],(*eta)[Mu2],(*phi)[Mu2],Muon_Mass);
      Float_t  InvarMass=(P1+P2).M();
      DiMuonInvariantMass->push_back(InvarMass);
      Pc=P1+P2;
      Float_t Gamma=Pc.Gamma();
      TLorentzRotation l;
      l.Boost(Pc.Px()/Gamma/InvarMass,Pc.Py()/Gamma/InvarMass,Pc.Pz()/Gamma/InvarMass);
      l.Invert();
      P1.Transform(l);
      CosThetaStar->push_back(P1.CosTheta());
      AngleBetweenDiMuon->push_back( acos( -(P1.Px()*P2.Px()+P1.Py()*P2.Py()+P1.Pz()*P2.Pz())/P1.P()/P2.P() ) ); 
    }
  Muons_Tree->Fill();
  delete HepGenEvent;
  return;
}

void MuESelector::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup) {
  Bool_t changed=true;
  if (HLTConfig.init(iRun,iSetup,HLTProc,changed)) {
    // if init returns TRUE, initialisation has succeeded!
    if (changed) {
      // The HLT config has actually changed wrt the previous Run, hence rebook your
      // histograms or do anything else dependent on the revised HLT config
      Info.HLT_FirstRunNum =iRun.id ().run ();
      HLTTableName->assign(HLTConfig.tableName());
      *HLTNamesSet=HLTConfig.triggerNames();
      HLTNames_Tree->Fill();
      LogInfo("DataFormat")<<"HLT table is changed to"<<endl;
      HLTConfig.dump("TableName");
      for (UInt_t j=0;j<num_HLTsSaveObjs;j++) {
	UInt_t i,length=HLTObj_HLTNames[j].rfind("_v");
	if (length<6) length=HLTObj_HLTNames[j].size();
	HLTSize=HLTNamesSet->size();
	for (i=0;i<HLTSize;i++)
	  if ( HLTObj_HLTNames[j].compare(0,length,HLTNamesSet->at(i),0,length)==0 ) {HLTFilterNamesAcceptenceIndex[j]=i;break;}//ignore the versions of the HLT name
	if (i==HLTSize) {
	  LogWarning("DataLost")<<"the HLT object request:"<<HLTObj_HLTNames[j].substr(0,length)<<" does not exist in the current HLT scope"<<endl;
	  HLTFilterNamesAcceptenceIndex[j]=-1;
	}
	else //for (vector<string>::const_iterator iter=HLTConfig.moduleLabels(HLTFilterNamesAcceptenceIndex[j]).begin();iter!=HLTConfig.moduleLabels(HLTFilterNamesAcceptenceIndex[j]).end();iter++)
	  HLT_ModuleNames[j]=InputTag((HLTConfig.moduleLabels(HLTFilterNamesAcceptenceIndex[j]).end()-2)->c_str(),"",HLTProc.c_str());
      }
    }
  } else {
    // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
    // with the file and/or code and needs to be investigated!
    LogWarning("DataLost") << " HLT config extraction failure with process name " <<HLTProc<<endl;
    // In this case, all access methods will return empty values!
  }
}

#define MakeVecBranch(Name,Var,Type) Var=new vector<Type>();Muons_Tree->Branch(Name,&Var)
MuESelector::MuESelector(const edm::ParameterSet& pset) {

//###############Get the input parameters###################
//Standard filters
  HLTFilter=new TriggerResultsFilter(pset.getParameter<edm::ParameterSet>("HLTFilterPSet"));
  VertexFilter=new GoodVertexFilter(pset.getParameter<edm::ParameterSet>("VertexFilterPSet"));
  ScrapingFilter=new FilterOutScraping(pset.getParameter<edm::ParameterSet>("noscrapingPSet"));
//Output filename
  FileName = pset.getParameter<string>("FileName");
//Muon Cut
  MuonPtCut = pset.getUntrackedParameter<double>("MuonPtCut",0.);
  //DiMuonInvarMassCut = pset.getUntrackedParameter<double>("DiMuonInvarMassCut",0.);
//chamber match
  maxChamberDist = pset.getUntrackedParameter<double>("maxChamberDist",-3.);
  maxChamberDistPull = pset.getUntrackedParameter<double>("maxChamberDistPull",-3.);
  //  HasTrackingParticle = pset.getUntrackedParameter<bool>("HasTrackingParticle",false);
//standard muon types e.g. TMLastStation......
  vector<string> Cuts;
  Cuts = pset.getUntrackedParameter< vector<string> >("StandardMuonCuts",Cuts);//default is empty
  string MuonArbitrationTypeStr;
  MuonArbitrationTypeStr = pset.getUntrackedParameter<string>("MuonArbitrationType","SegmentArbitration");
  MuonArbitrationType=MuonArbitrationTypeFromString(MuonArbitrationTypeStr.c_str());
//PV
  PrimaryVerticesTag = pset.getUntrackedParameter<string>("PrimaryVertices","offlinePrimaryVertices");//"offlinePrimaryVerticesWithBS": Primary vertex reconstructed using the tracks taken from the generalTracks collection, and imposing the offline beam spot as a constraint in the fit of the vertex position. Another possible tag is "offlinePrimaryVertices", which is Primary vertex reconstructed using the tracks taken from the generalTracks collection
//tracks
  tracksTag = pset.getUntrackedParameter<InputTag>("tracksTag",InputTag("generalTracks"));
  dEdxTag = pset.getUntrackedParameter<string>("dEdxTag","dedxHarmonic2");//Other options already available in RECO files are dedxMedian and dedxTruncated40. 
//Triggers
  TriggerResultsTag = pset.getParameter<InputTag>("TriggerResultsTag");
  HLTProc=TriggerResultsTag.process();
  TriggerEventTag = pset.getUntrackedParameter<InputTag>("triggerEventTag",InputTag("hltTriggerSummaryAOD","",HLTProc.c_str()));
  HLTObj_HLTNames = pset.getParameter< vector<string> >("HLTObj_HLTNames");
  num_HLTsSaveObjs = HLTObj_HLTNames.size();
  HLT_ModuleNames.assign(num_HLTsSaveObjs,InputTag(""));
  if (num_HLTsSaveObjs>maxFilterObjects) throw cms::Exception("You need to increase the constant maxFilterObjects.\n");;
//Generation Information
  HepMCTag = pset.getUntrackedParameter<string>("HepMCTag","generator");

//#################### Make a root file and plant TTrees #########################
  file=new TFile(FileName.c_str(),"RECREATE");
  Muons_Tree = new TTree ("MuTrkCand","TrksCandidates") ;
  HLTNames_Tree = new TTree ("HLTNames","HLTNames") ;
  Muons_Tree->SetCircular(500000);
  HLTNames_Tree->SetCircular(500000);//the max events in a single root file
//#################### Build Branches #########################
  //General event information
  Muons_Tree->Branch("Event_Info",&Info.RUN,"RUN/l:EVENT:LumiBlock:ORBIT:BunchCrossing:HLT_FirstRunNum");
  Muons_Tree->Branch("isRealData",&isRealData,"isRealData/O");
  Muons_Tree->Branch("isHLTTriggerred",&isHLTTriggerred,"isHLTTriggerred/O");
  Muons_Tree->Branch("isGoodVertex",&isGoodVertex,"isGoodVertex/O");
  Muons_Tree->Branch("isNoScrapping",&isNoScrapping,"isNoScrapping/O");
  Muons_Tree->Branch("GenEventWeight",&GenEventWeight,"GenEventWeight/D");
  numberOfPUVertices=0;
  Muons_Tree->Branch("numberOfPUVertices",&numberOfPUVertices,"numberOfPUVertices/I");
  //Muon Basic Information
  MakeVecBranch("pt",pt,Float_t);  MakeVecBranch("eta",eta,Float_t);  MakeVecBranch("phi",phi,Float_t);
  MakeVecBranch("chargeMinus",chargeMinus,Bool_t); MakeVecBranch("isGlobalMu",isGlobalMu,Bool_t);  MakeVecBranch("isTrackerMu",isTrackerMu,Bool_t);
  //Muon Vertex
  MakeVecBranch("Vertex_x",Vertex_x,Float_t);  MakeVecBranch("Vertex_y",Vertex_y,Float_t);  MakeVecBranch("Vertex_z",Vertex_z,Float_t);
  MakeVecBranch("DXYwtBS",DXYwtBS,Float_t); MakeVecBranch("DXYwtPV",DXYwtPV,Float_t); MakeVecBranch("DXYErrwtPV",DXYErrwtPV,Float_t);
  //Muon Isolation
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
  //trackquality
  MakeVecBranch("numberOfMatchedSegments",numberOfMatchedSegments,Int_t);
  MakeVecBranch("numberOfMatchedStations",numberOfMatchedStations,Int_t);
  MakeVecBranch("dEdx",dEdx,Float_t);  MakeVecBranch("dEdxError",dEdxError,Float_t);
  MakeVecBranch("dEdx_numberOfSaturatedMeasurements",dEdx_numberOfSaturatedMeasurements,Int_t);
  MakeVecBranch("dEdx_numberOfMeasurements",dEdx_numberOfMeasurements,Int_t);
  MakeVecBranch("Track_nValidTrackerHits",InnerTrack_nValidTrackerHits,UInt_t);
  MakeVecBranch("Track_nValidPixelHits",InnerTrack_nValidPixelHits,UInt_t);
  MakeVecBranch("Track_nLostTrackerHits",InnerTrack_nLostTrackerHits,UInt_t);
  MakeVecBranch("Track_nLostPixelHits",InnerTrack_nLostPixelHits,UInt_t);
  MakeVecBranch("numberOfValidMuonHits",numberOfValidMuonHits,UInt_t);
  MakeVecBranch("InnerTrack_chi2",InnerTrack_chi2,Float_t);  MakeVecBranch("InnerTrack_ndof",InnerTrack_ndof,UInt_t);
  MakeVecBranch("GlobalTrack_chi2",GlobalTrack_chi2,Float_t);  MakeVecBranch("GlobalTrack_ndof",GlobalTrack_ndof,UInt_t);
  //Information of two muons
  MakeVecBranch("DiMuonInvariantMass",DiMuonInvariantMass,Float_t);
  MakeVecBranch("CosThetaStar",CosThetaStar,Float_t);
  MakeVecBranch("AngleBetweenDiMuon",AngleBetweenDiMuon,Float_t);
  //Muon System - Segment/Chambers
  MakeVecBranch("TrackDistToChamberEdge",TrackDistToChamberEdge,Float_t);
  MakeVecBranch("TrackDistToChamberEdgeErr",TrackDistToChamberEdgeErr,Float_t);
  MakeVecBranch("XTrack",XTrack,Float_t);
  MakeVecBranch("YTrack",YTrack,Float_t);
  MakeVecBranch("XErrTrack",XErrTrack,Float_t);
  MakeVecBranch("YErrTrack",YErrTrack,Float_t);
  MakeVecBranch("XSegment",XSegment,Float_t);
  MakeVecBranch("YSegment",YSegment,Float_t);
  MakeVecBranch("XErrSegment",XErrSegment,Float_t);
  MakeVecBranch("YErrSegment",YErrSegment,Float_t);
  MakeVecBranch("IsCSCChamber",IsCSCChamber,Bool_t);
  MakeVecBranch("SectorChamber",SectorChamber,Byte_t);
  MakeVecBranch("StationRing",StationRing,Char_t);
  MakeVecBranch("IsSegmentOwnedExclusively",IsSegmentOwnedExclusively,Bool_t);
  MakeVecBranch("IsSegmentBelongsToTrackByDR",IsSegmentBelongsToTrackByDR,Bool_t);
  MakeVecBranch("IsSegmentBelongsToTrackByCleaning",IsSegmentBelongsToTrackByCleaning,Bool_t);
  MakeVecBranch("NumberOfHitsInSegment",NumberOfHitsInSegment,Byte_t);
  MakeVecBranch("MuonIndex",MuonIndex,Byte_t);
  MakeVecBranch("StationMask",StationMask,UInt_t);MakeVecBranch("RequiredStationMask",RequiredStationMask,UInt_t);
  //HLT
  MakeVecBranch("HLTacceptance",HLTacceptance,Bool_t);
  for( unsigned int i =0; i < num_HLTsSaveObjs;i++) {
    MakeVecBranch((HLTObj_HLTNames[i]+"_pt").c_str(),HLTObj_pt[i],Float_t);
    MakeVecBranch((HLTObj_HLTNames[i]+"_eta").c_str(),HLTObj_eta[i],Float_t);
    MakeVecBranch((HLTObj_HLTNames[i]+"_phi").c_str(),HLTObj_phi[i],Float_t);
    MakeVecBranch(("is"+HLTObj_HLTNames[i]+"Obj").c_str(),isHLTObj[i],Bool_t);
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
  
  //Generation and Simulation Particles
  MakeVecBranch("Gen_pt",Gen_pt,Float_t);  MakeVecBranch("Gen_eta",Gen_eta,Float_t);  MakeVecBranch("Gen_phi",Gen_phi,Float_t);   MakeVecBranch("Gen_pdgId",Gen_pdgId,Int_t);
  MakeVecBranch("Gen_vx",Gen_vx,Float_t);  MakeVecBranch("Gen_vy",Gen_vy,Float_t);  MakeVecBranch("Gen_vz",Gen_vz,Float_t);  MakeVecBranch("Gen_vt",Gen_vt,Float_t);
  Muons_Tree->Branch("HepMCVertex",HepMCFourVec,"HepMCVertex[4]/F");

#ifdef TrackingParticles  
  //Simulated Tracks
  minTrackHits = pset.getUntrackedParameter<uint>("minTrackHits",3);
  MakeVecBranch("IsParInHep",IsParInHep,Bool_t);
  MakeVecBranch("TrkParticles_pt",TrkParticles_pt,Float_t);  MakeVecBranch("TrkParticles_eta",TrkParticles_eta,Float_t);  MakeVecBranch("TrkParticles_phi",TrkParticles_phi,Float_t);
  MakeVecBranch("TrkParticles_pdgId",TrkParticles_pdgId,Int_t); MakeVecBranch("TrkParticles_charge",TrkParticles_charge,Int_t);
  MakeVecBranch("SharedHitsRatio",SharedHitsRatio,Double_t); MakeVecBranch("MCMatchChi2",MCMatchChi2,Double_t); 
  MakeVecBranch("DChains",DChains,Int_t); MakeVecBranch("theSameWithMuon",theSameWithMuon,Int_t);
  NumMisMatch=0;
#endif

  //Summary
  HLTTableName = new string();
  HLTNames_Tree->Branch("HLTScopeFirstRun",&Info.HLT_FirstRunNum,"HLT_FirstRunNum/l");
  HLTNames_Tree->Branch("TableName",&HLTTableName);
  HLTNamesSet = new vector<string>();
  HLTNames_Tree->Branch("HLTNamesSet",&HLTNamesSet);
  FirstEvent=true;IsMC=true;
}

Muon::ArbitrationType MuESelector::MuonArbitrationTypeFromString( const std::string &label ) {
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

void MuESelector::ClearVecs_HepMC() {
  Gen_pt->clear();  Gen_eta->clear();  Gen_phi->clear(); Gen_pdgId->clear();
  Gen_vx->clear();  Gen_vy->clear();  Gen_vz->clear();  Gen_vt->clear();
  ParToHep.clear();
}

void MuESelector::ClearVecs_RECO() {
  //basic information
  pt->clear(); eta->clear();  phi->clear();
  chargeMinus->clear();  isGlobalMu->clear();  isTrackerMu->clear();
  //muon vertex
  Vertex_x->clear();  Vertex_y->clear();  Vertex_z->clear();
  DXYwtBS->clear();  DXYwtPV->clear(); DXYErrwtPV->clear();
  //isolation
  isoR03sumPt->clear();  isoR03emEt->clear();  isoR03hadEt->clear();
  isoR03hoEt->clear();  isoR03nJets->clear();  isoR03nTracks->clear();
  isoR05sumPt->clear();  isoR05emEt->clear();  isoR05hadEt->clear();
  isoR05hoEt->clear();  isoR05nJets->clear();  isoR05nTracks->clear();
  isoemVetoEt->clear();  isohadVetoEt->clear();  isohoVetoEt->clear();
  CaloE_emMax->clear();  CaloE_emS9->clear();  CaloE_emS25->clear();
  CaloE_hadMax->clear();  CaloE_hadS9->clear();
  Calo_emPos_R->clear();  Calo_emPos_eta->clear();  Calo_emPos_phi->clear();
  Calo_hadPos_R->clear();  Calo_hadPos_eta->clear();  Calo_hadPos_phi->clear();
  //track quality
  TrkKink->clear();  GlbKink->clear();  TrkRelChi2->clear();
  dEdx->clear();  dEdxError->clear();
  dEdx_numberOfSaturatedMeasurements->clear();
  dEdx_numberOfMeasurements->clear();
  InnerTrack_nValidTrackerHits->clear();
  InnerTrack_nValidPixelHits->clear();  
  InnerTrack_nLostTrackerHits->clear();
  InnerTrack_nLostPixelHits->clear(); 
  InnerTrack_chi2->clear(); InnerTrack_ndof->clear();
  numberOfValidMuonHits->clear();
  numberOfMatchedSegments->clear();
  numberOfMatchedStations->clear();
  GlobalTrack_chi2->clear(); GlobalTrack_ndof->clear();
  //Information of two muons
  DiMuonInvariantMass->clear();  CosThetaStar->clear(); AngleBetweenDiMuon->clear();

  //HLT
  HLTacceptance->clear();
  for( unsigned int i =0; i < num_HLTsSaveObjs;i++) {
    HLTObj_pt[i]->clear();
    HLTObj_eta[i]->clear();
    HLTObj_phi[i]->clear();
    isHLTObj[i]->clear();
  }

  //Chamber
  TrackDistToChamberEdge->clear();  TrackDistToChamberEdgeErr->clear();
  XTrack->clear();YTrack->clear();  XErrTrack->clear();YErrTrack->clear();
  IsCSCChamber->clear(); StationRing->clear();SectorChamber->clear();MuonIndex->clear();

  //Segment
  XSegment->clear();  YSegment->clear();
  XErrSegment->clear();  YErrSegment->clear();
  IsSegmentOwnedExclusively->clear();
  IsSegmentBelongsToTrackByDR->clear();  IsSegmentBelongsToTrackByCleaning->clear();
  NumberOfHitsInSegment->clear();
  StationMask->clear();RequiredStationMask->clear();

  vx->clear(); vy->clear(); vz->clear();
  vxError->clear(); vyError->clear(); vzError->clear();

  for (unsigned int whichcut=0;whichcut<num_Cuts;whichcut++)
    SelectorPassed[whichcut]->clear(); 
  MySelector->clear();

#ifdef TrackingParticles
  TrkParticles_pt->clear();  TrkParticles_eta->clear(); TrkParticles_phi->clear();
  TrkParticles_pdgId->clear(); TrkParticles_charge->clear(); 
  SharedHitsRatio->clear(); MCMatchChi2->clear(); ParToSim.clear(); 
  DChains->clear(); theSameWithMuon->clear();
  IsParInHep->clear();
#endif
}

// ------- method called once each job just before starting event loop  ------------
void MuESelector::beginJob(){}

// -------- method called once each job just after ending the event loop  ------------
void MuESelector::endJob() {
  file->Write();
  file->Close();
}

MuESelector::~MuESelector() {
  //standard filters
  HLTFilter->~TriggerResultsFilter();
  VertexFilter->~GoodVertexFilter();
  ScrapingFilter->~FilterOutScraping();
  //Generation Information
  delete Gen_pt; delete Gen_eta; delete Gen_phi; delete Gen_pdgId;
  delete Gen_vx; delete Gen_vy;  delete Gen_vz;  delete Gen_vt;
 //basic information
  delete pt; delete eta; delete  phi;
  delete chargeMinus;
  delete  isGlobalMu; delete  isTrackerMu;
  //muon vertex
  delete Vertex_x; delete  Vertex_y; delete  Vertex_z;
  delete DXYwtBS; delete  DXYwtPV; delete DXYErrwtPV;
  //isolation
  delete isoR03sumPt; delete isoR03emEt; delete isoR03hadEt;
  delete isoR03hoEt; delete isoR03nJets; delete isoR03nTracks;
  delete isoR05sumPt; delete isoR05emEt; delete isoR05hadEt;
  delete isoR05hoEt; delete isoR05nJets; delete isoR05nTracks;
  delete isoemVetoEt; delete isohadVetoEt; delete isohoVetoEt;
  delete CaloE_emMax; delete CaloE_emS9; delete CaloE_emS25;
  delete CaloE_hadMax; delete CaloE_hadS9; delete Calo_emPos_R;
  delete Calo_emPos_eta; delete Calo_emPos_phi; delete Calo_hadPos_R;
  delete Calo_hadPos_eta; delete Calo_hadPos_phi;
  //track quality
  delete TrkKink; delete GlbKink; delete TrkRelChi2;
  delete dEdx; delete  dEdxError;
  delete dEdx_numberOfSaturatedMeasurements;
  delete dEdx_numberOfMeasurements;
  delete InnerTrack_nValidTrackerHits;
  delete InnerTrack_nValidPixelHits;
  delete InnerTrack_nLostTrackerHits;
  delete InnerTrack_nLostPixelHits;
  delete InnerTrack_chi2; delete InnerTrack_ndof;
  delete numberOfValidMuonHits;
  delete numberOfMatchedSegments;
  delete numberOfMatchedStations;
  delete GlobalTrack_chi2; delete GlobalTrack_ndof;
  //Information of two muons
  delete DiMuonInvariantMass;
  delete  CosThetaStar; delete AngleBetweenDiMuon;

  //HLT
  delete HLTTableName;
  delete HLTacceptance;
  for( unsigned int i =0; i < num_HLTsSaveObjs;i++) {
    delete HLTObj_pt[i];
    delete HLTObj_eta[i];
    delete HLTObj_phi[i];
    delete isHLTObj[i];
  }

  //Chamber
  delete TrackDistToChamberEdge; delete TrackDistToChamberEdgeErr;
  delete XTrack; delete YTrack; delete XErrTrack; delete YErrTrack;
  delete IsCSCChamber; delete StationRing; delete SectorChamber;
  delete MuonIndex;

  //Segment
  delete XSegment; delete  YSegment;
  delete XErrSegment; delete  YErrSegment;
  delete IsSegmentOwnedExclusively; delete IsSegmentBelongsToTrackByDR;
  delete IsSegmentBelongsToTrackByCleaning; delete NumberOfHitsInSegment;
  delete StationMask; delete RequiredStationMask;
  
  delete vx; delete vy; delete vz;
  delete vxError; delete vyError; delete vzError;

  for (unsigned int whichcut=0;whichcut<num_Cuts;whichcut++)
    delete SelectorPassed[whichcut];
  delete MySelector;
#ifdef TrackingParticles
  delete TrkParticles_pt; delete TrkParticles_eta; delete TrkParticles_phi;
  delete TrkParticles_pdgId; delete  TrkParticles_charge;
  delete SharedHitsRatio; delete MCMatchChi2; 
  delete DChains; delete theSameWithMuon;
  delete IsParInHep; 
#endif
}

#ifdef TrackingParticles
void MuESelector::GetDecayChains(TrackingParticleRef tpr,vector<int> *DChains, vector <TheMuonType> &type, HepMC::GenEvent &HepGenEvent) {
  for (vector<SimTrack>::const_iterator g4Track_iter = tpr->g4Track_begin(); g4Track_iter != tpr->g4Track_end(); ++g4Track_iter )  {//g4Track loop begin
    DChain.clear();SimChains.clear();
    const SimTrack *thisTrk=&(*g4Track_iter);
    SimTrackDaughtersTree( thisTrk,tpr );
    Bool_t ChainEnd; TrackingParticleRef tpr_tmp=tpr;
    do {
      ChainEnd=true;
      if ( !thisTrk->noVertex() ) {
	TrackingVertexRef tvr=tpr_tmp->parentVertex();
	for ( TrackingParticleRefVector::iterator parenttp=tvr->sourceTracks_begin();parenttp!=tvr->sourceTracks_end();parenttp++ ) {
	  for (vector<SimTrack>::const_iterator g4Trk_iter = (*parenttp)->g4Track_begin() ; g4Trk_iter != (*parenttp)->g4Track_end(); ++g4Trk_iter )
	    if ( SVC[thisTrk->vertIndex()].parentIndex()==Int_t(g4Trk_iter->trackId()) && g4Trk_iter->eventId().rawId()==thisTrk->eventId().rawId()) { 
	      thisTrk=&(*g4Trk_iter);tpr_tmp=*parenttp;
	      Bool_t HasMuonHits=false;
	      for ( vector<PSimHit>::const_iterator g4Hit_iter=tpr_tmp->pSimHit_begin();g4Hit_iter!=tpr_tmp->pSimHit_end();g4Hit_iter++ )
		if ( g4Hit_iter->trackId()==thisTrk->trackId() && g4Hit_iter->eventId().rawId()==thisTrk->eventId().rawId() && g4Hit_iter->particleType() == thisTrk->type() )      HasMuonHits=SimHitToSegment( *g4Hit_iter);
	      vector<const SimTrack *>::iterator SavedSimTrk_iter=find(SavedSimTrk.begin(),SavedSimTrk.end(),thisTrk);
	      if (SavedSimTrk_iter==SavedSimTrk.end())
		{
		  DChain.push_back(IsParInHep->size());
		  IsParHasMuonHits->push_back(HasMuonHits);
		  RecordSimTrack(thisTrk)
		    }
	      else DChain.push_back(FindSimTrackRecordingPosition(SavedSimTrk_iter-SavedSimTrk.begin()));
	      SavedTP.push_back(tpr_tmp);
	      ChainEnd=false;
	      break;
	    }
	  if (!ChainEnd) break;
	}
      }
    }while (!ChainEnd);
    for (vector< vector<Int_t> >::iterator SimChains_iter = SimChains.begin(); SimChains_iter !=  SimChains.end(); ++SimChains_iter )
      SimChains_iter->insert(SimChains_iter->begin(),DChain.rbegin(),DChain.rend());
    DChain.clear();
    
    //HepMC Particles
    HepMCChains.clear();
    if (!thisTrk->noGenpart()) {
      vector<const SimTrack *>::iterator SavedSimTrk_iter=find(SavedSimTrk.begin(),SavedSimTrk.end(),thisTrk);
      if (SavedSimTrk_iter!=SavedSimTrk.end()) (*IsParInHep)[FindSimTrackRecordingPosition(SavedSimTrk_iter-SavedSimTrk.begin())]=true;
      else LogError("CodeWrong")<<"Cannot find the simulated track in saved sim tracks."<<endl;
      HepMC::GenEvent::particle_iterator genPar = HepGenEvent->particles_begin();
      for (int count=1; count<thisTrk->genpartIndex()&&genPar != HepGenEvent->particles_end(); count++ )
	genPar++;
      if (genPar != HepGenEvent->particles_end()) HepMCParentTree(*genPar);
      else LogWarning("RefNull")<<"genpartIndex() Error or HepMC is empty"<<endl;
    }
    
#define SaveAndClassifyDChain(thisChain) int Muref=-1;			\
    for (vector<int>::iterator DChain_iter = DChains->begin(); DChain_iter != DChains->end()&&Muref<(int) eta->size()-1; DChain_iter++ ) { \
      if (*DChain_iter==-1) {						\
	DChain_iter++;							\
	vector<int>::iterator DChain_begin=DChain_iter;			\
	for (; DChain_iter != DChains->end()&&*DChain_iter!=-2&&*DChain_iter!=-1; DChain_iter++ ) ; \
	if (DChain_begin!=DChain_iter)					\
	  if (IstheSameDChain(thisChain,* new vector<int> (DChain_begin,DChain_iter))) { \
	    theSameWithMuon->back()=Muref;				\
	    break;							\
	  }								\
	DChain_iter--;							\
      }									\
      if (*DChain_iter==-2) Muref++;					\
    }									\
    DChains->push_back(-1);						\
    DChains->insert(DChains->end(),thisChain.begin(),thisChain.end());	\
    TheMuonType newtype=classification(thisChain);			\
    if ( type.empty() ) type.push_back(newtype);			\
    else if ( type.size()==1 && type.back()==Others ) type.back()=newtype; \
    else if ( type.size()==1 && type.back()==NoMuSysHit && newtype!=Others ) type.back()=newtype; \
    else if (find(type.begin(),type.end(),newtype)==type.end() && newtype!=Others && newtype!=NoMuSysHit ) type.push_back(newtype)
    
    //merge the HepMC and SimTrack Decay Chains
    for (vector< vector<Int_t> >::iterator SimChains_iter = SimChains.begin(); SimChains_iter !=  SimChains.end(); ++SimChains_iter )
      if ( !HepMCChains.empty() )
	for (vector< vector<Int_t> >::iterator HepMCChains_iter = HepMCChains.begin(); HepMCChains_iter !=  HepMCChains.end(); ++HepMCChains_iter ) {
	  vector<Int_t> thisChain(HepMCChains_iter->rbegin(),HepMCChains_iter->rend());
	  thisChain.insert(thisChain.end(),SimChains_iter->begin(),SimChains_iter->end());
	  SaveAndClassifyDChain(thisChain);
	}
      else {SaveAndClassifyDChain((*SimChains_iter));}
  }//g4Track loop end
}

bool MuESelector::IstheSameDChain(vector<int> &ThisChain,vector<int> &AnotherChain) {//whether two chains include each other
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

void MuESelector::SimTrackDaughtersTree(SimTrack * thisTrk)
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

void MuESelector::HepMCParentTree(HepMC::GenParticle *genPar) {
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
//define this as a plug-in
DEFINE_FWK_MODULE(MuESelector);
