#include "HGCal/Reco/plugins/HGCalTBRecHitProducer.h"
#include "HGCal/Reco/interface/PulseFitter.h"
#include "HGCal/Reco/interface/CommonMode.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"

#include <iostream>


const static int SENSORSIZE = 128;

HGCalTBRecHitProducer::HGCalTBRecHitProducer(const edm::ParameterSet& cfg) : 
  m_CommonModeNoiseCollectionName(cfg.getUntrackedParameter<std::string>("CommonModeNoiseCollectionName", "HGCALTBCOMMONMODENOISEMAP")),
  m_outputCollectionName(cfg.getParameter<std::string>("OutputCollectionName")),
  m_electronicMap(cfg.getUntrackedParameter<std::string>("ElectronicsMap","HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers_AllFlipped.txt")),
  m_detectorLayoutFile(cfg.getUntrackedParameter<std::string>("DetectorLayout","HGCal/CondObjects/data/layerGeom_oct2017_h2_17layers.txt")),
  m_adcCalibrationsFile(cfg.getUntrackedParameter<std::string>("ADCCalibrations","HGCal/CondObjects/data/hgcal_calibration.txt")),
  m_timeSample3ADCCut(cfg.getUntrackedParameter<double>("TimeSample3ADCCut",15)),
  m_maskNoisyChannels(cfg.getUntrackedParameter<bool>("MaskNoisyChannels",false)),
  m_channelsToMask_filename(cfg.getUntrackedParameter<std::string>("ChannelsToMaskFileName","HGCal/CondObjects/data/noisyChannels.txt")),
  m_NHexaBoards(cfg.getUntrackedParameter<int>("NHexaBoards", 10)) { 
  
  m_HGCalTBRawHitCollection = consumes<HGCalTBRawHitCollection>(cfg.getParameter<edm::InputTag>("InputCollection"));
  RunDataToken= consumes<RunData>(cfg.getParameter<edm::InputTag>("RUNDATA"));
  HGCalTBGlobalTimestampsToken= consumes<HGCalTBGlobalTimestamps>(cfg.getParameter<edm::InputTag>("GlobalTimestampCollectionName"));
  investigatePulseShape = cfg.getUntrackedParameter<bool>("investigatePulseShape", false);
  std::cout << cfg.dump() << std::endl;

  produces <std::map<int, commonModeNoise> >(m_CommonModeNoiseCollectionName);
  produces <HGCalTBRecHitCollection>(m_outputCollectionName);
}

void HGCalTBRecHitProducer::beginJob()
{
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_electronicMap);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  };

  edm::Service<TFileService> fs;


  std::ostringstream os( std::ostringstream::ate );
  TH2F* htmp2;
  for(int ib = 0; ib<m_NHexaBoards; ib++) {
    for( size_t iski=0; iski<HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; iski++ ){
      os.str("");os<<"HexaBoard"<<ib<<"_Skiroc"<<iski;
      TFileDirectory dir = fs->mkdir( os.str().c_str() );
      for( size_t ichan=0; ichan<HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; ichan++ ){
        if ((ichan % 2) == 1) continue;

        int key = ib * 10000 + iski * 100 + ichan;

        os.str("");os<<"Channel"<<ichan<<"__LGShape";
        shapesLG[key] = dir.make<TH2F>(os.str().c_str(),os.str().c_str(), 100, 0, 225, 100, -0.2, 1.);
        os.str("");os<<"Channel"<<ichan<<"__HGShape";
        shapesHG[key] = dir.make<TH2F>(os.str().c_str(),os.str().c_str(), 100, 0, 225, 100, -0.2, 1.);
        
        os.str("");os<<"Channel"<<ichan<<"__ToARiseVsTMaxLG";
        ToARisevsTMaxLG[key] = dir.make<TH2F>(os.str().c_str(),os.str().c_str(), 100, 50., 150., 100, 4., 3500.);
        os.str("");os<<"Channel"<<ichan<<"__ToARiseVsTMaxHG";
        ToARisevsTMaxHG[key] = dir.make<TH2F>(os.str().c_str(),os.str().c_str(), 100, 50., 150., 100, 4., 3500.);

        os.str("");os<<"Channel"<<ichan<<"__ToAFallVsTMaxLG";
        ToAFallvsTMaxLG[key] = dir.make<TH2F>(os.str().c_str(),os.str().c_str(), 100, 50., 150., 100, 4., 3500.);
        os.str("");os<<"Channel"<<ichan<<"__ToAFallVsTMaxHG";
        ToAFallvsTMaxHG[key] = dir.make<TH2F>(os.str().c_str(),os.str().c_str(), 100, 50., 150., 100, 4., 3500.);
        
        
        os.str("");os<<"Channel"<<ichan<<"__TMaxHGVsTMaxLG";
        TMaxHGvsTMaxLG[key] = dir.make<TH2F>(os.str().c_str(),os.str().c_str(), 100, 50., 150., 100, 50., 150.);
      }
      TFileDirectory gaindir = dir.mkdir( "Gains" );
      htmp2=gaindir.make<TH2F>("HighGainVsLowGainAmpl","HighGainVsLowGainAmpl",200,-500.,1500,200,-500.,3500);
      htmp2->GetXaxis()->SetTitle("Low Gain Amplitude [ADC]"); 
      htmp2->GetYaxis()->SetTitle("High Gain Amplitude [ADC]"); 
      m_h_HighVsLowGainAmpl.insert( std::pair<int,TH2F*>(ib*10+iski, htmp2) );  

      htmp2=gaindir.make<TH2F>("LowGainVsTOTAmpl","LowGainVsTOTAmpl",150, 0.,1500,175,-500.,3000);
      htmp2->GetXaxis()->SetTitle("TOT [ADC]");   
      htmp2->GetYaxis()->SetTitle("Low Gain Amplitude [ADC]"); 
      m_h_LowGainVsTOTAmpl.insert( std::pair<int,TH2F*>(ib*10+iski, htmp2) ); 
    }
  }

  outtree = fs->make<TTree>("timeInvestigationTree", "timeInvestigationTree");

  outtree->Branch("event", &event_for_tree);
  outtree->Branch("run", &run_for_tree);
  outtree->Branch("pdgID", &pdgID_for_tree);
  outtree->Branch("beamEnergy", &beamEnergy_for_tree);
  outtree->Branch("board", &board_for_tree);
  outtree->Branch("skiroc", &skiroc_for_tree);
  outtree->Branch("channel", &channel_for_tree);
  outtree->Branch("globalTimestamp", &globalTimestamp_for_tree);
  outtree->Branch("MIPEnergy", &MIPEnergy_for_tree);
  outtree->Branch("HG_max", &HG_max_for_tree);
  outtree->Branch("LG_max", &LG_max_for_tree);
  outtree->Branch("TMax_HG", &TMax_HG_for_tree);
  outtree->Branch("TMax_LG", &TMax_LG_for_tree);
  outtree->Branch("TOT", &TOT_for_tree);
  outtree->Branch("TOA_rise", &TOA_rise_for_tree);
  outtree->Branch("TOA_fall", &TOA_fall_for_tree);
  outtree->Branch("chi2_HG", &chi2_HG_for_tree);
  outtree->Branch("chi2_LG", &chi2_LG_for_tree);
  outtree->Branch("trise_HG", &trise_HG_for_tree);
  outtree->Branch("trise_LG", &trise_LG_for_tree);
  outtree->Branch("errortmax_HG", &errortmax_HG_for_tree);
  outtree->Branch("errortmax_LG", &errortmax_LG_for_tree);
  outtree->Branch("erroramplitude_HG", &erroramplitude_HG_for_tree);
  outtree->Branch("erroramplitude_LG", &erroramplitude_LG_for_tree);

  if( m_maskNoisyChannels ){
    FILE* file;
    char buffer[300];
    //edm::FileInPath fip();
    file = fopen (m_channelsToMask_filename.c_str() , "r");
    if (file == NULL){
      perror ("Error opening noisy channels file"); exit(1); 
    } else{
    
      while ( ! feof (file) ){
        if ( fgets (buffer , 300 , file) == NULL ) break;
        const char* index = buffer;
        int layer,skiroc,channel,ptr,nval;
        nval=sscanf( index, "%d %d %d %n",&layer,&skiroc,&channel,&ptr );
        int skiId=HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA*layer+(HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA-skiroc)%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA+1;
        if( nval==3 ){
          HGCalTBElectronicsId eid(skiId,channel);      
          if (essource_.emap_.existsEId(eid.rawId()))
            m_noisyChannels.push_back(eid.rawId());
        } else continue;
      }
    }
    fclose (file);
  }

  fip=edm::FileInPath(m_detectorLayoutFile);
  if (!io.load(fip.fullPath(), essource_.layout_)) {
    throw cms::Exception("Unable to load detector layout file");
  };
  for( auto layer : essource_.layout_.layers() )
    layer.print();
  
  fip=edm::FileInPath(m_adcCalibrationsFile);
  if (!io.load(fip.fullPath(), essource_.adccalibmap_)) {
    throw cms::Exception("Unable to load ADC conversions map");
  };
}

void HGCalTBRecHitProducer::produce(edm::Event& event, const edm::EventSetup& iSetup)
{
  edm::Handle<RunData> rd;
  //get the relevant event information
  event.getByToken(RunDataToken, rd);


  edm::Handle<HGCalTBGlobalTimestamps> gts;
  event.getByToken(HGCalTBGlobalTimestampsToken, gts);

  std::auto_ptr<HGCalTBRecHitCollection> rechits(new HGCalTBRecHitCollection);

  edm::Handle<HGCalTBRawHitCollection> rawhits;
  event.getByToken(m_HGCalTBRawHitCollection, rawhits);

  //compute the common mode noise 
  CommonMode cm(essource_.emap_, true, true, -1); //default is common mode per chip using the median
  cm.Evaluate( rawhits );

  //..and at it to the edm event
  std::auto_ptr<std::map<int, commonModeNoise> > cmMap(new std::map<int, commonModeNoise>);
  (*cmMap)=cm.CommonModeNoiseMap();
  

  std::vector<std::pair<double, double> > CellXY;


  PulseFitter fitter(0,150);
  PulseFitterResult fitresultLG;
  PulseFitterResult fitresultHG;

  for( auto rawhit : *rawhits ){

  event_for_tree = rd->event;
  run_for_tree = rd->run;
  pdgID_for_tree = rd->pdgID;
  beamEnergy_for_tree = rd->energy;



    HGCalTBElectronicsId eid( essource_.emap_.detId2eid(rawhit.detid().rawId()) );
    if( !essource_.emap_.existsEId(eid.rawId()) || std::find(m_noisyChannels.begin(),m_noisyChannels.end(),eid.rawId())!=m_noisyChannels.end() )
      continue;
    
    int iski=rawhit.skiroc();
    int iboard=iski/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
    int ichannel=eid.ichan();
    board_for_tree = iboard;
    skiroc_for_tree = iski;
    channel_for_tree = ichannel;
    int key = iboard * 10000 + (iski % 4) * 100 + ichannel;

    uint32_t globalTimestamp = gts->skiroc_to_timestamps.at(iski);

    globalTimestamp_for_tree = globalTimestamp;

    MIPEnergy_for_tree =  HG_max_for_tree =  LG_max_for_tree =  TMax_HG_for_tree =  TMax_LG_for_tree =  TOT_for_tree =  TOA_rise_for_tree =  TOA_fall_for_tree =  chi2_HG_for_tree =  chi2_LG_for_tree =  trise_HG_for_tree =  trise_LG_for_tree =  errortmax_HG_for_tree =  errortmax_LG_for_tree =  erroramplitude_HG_for_tree =  erroramplitude_LG_for_tree = -999;

    std::vector<double> sampleHG, sampleLG, sampleT;

    float highGain(0), lowGain(0), totGain(0),toaRise(0), toaFall(0), energy(0), _time(-1);

    float subHG[NUMBER_OF_TIME_SAMPLES],subLG[NUMBER_OF_TIME_SAMPLES];

    totGain = rawhit.totSlow();
    toaRise = rawhit.toaRise();
    toaFall = rawhit.toaFall();

    switch ( rawhit.detid().cellType() ){
    default :
      for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	subHG[it]=0;
      	subLG[it]=0;
      }
    case 0 : 
      for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	subHG[it]=cmMap->at(iski).fullHG[it]; 
      	subLG[it]=cmMap->at(iski).fullLG[it]; 
      }
      break;
    case 2 : 
      for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	subHG[it]=cmMap->at(iski).halfHG[it]; 
      	subLG[it]=cmMap->at(iski).halfLG[it]; 
      }
      break;
    case 3 : 
      for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	subHG[it]=cmMap->at(iski).mouseBiteHG[it]; 
      	subLG[it]=cmMap->at(iski).mouseBiteLG[it]; 
      }
      break;
    case 4 : for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      	subHG[it]=cmMap->at(iski).outerHG[it]; 
      	subLG[it]=cmMap->at(iski).outerLG[it]; 
      }
      break;
    case 5 : for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
        subHG[it]=cmMap->at(iski).mergedHG[it]; 
        subLG[it]=cmMap->at(iski).mergedLG[it]; 
      }
      break;
    }

    
    for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      sampleHG.push_back(rawhit.highGainADC(it)-subHG[it]);
      sampleLG.push_back((rawhit.lowGainADC(it)-subLG[it]));
      if( highGain<sampleHG[it] && it>1 && it<5 ) highGain=sampleHG[it];
      if( lowGain<sampleLG[it] && it>1 && it<5 ) lowGain=sampleLG[it];
      sampleT.push_back(25*it+12.5);
    }

    //this is a just try to isolate hits with signal
    float en1=sampleHG[1];
    float en2=sampleHG[2];
    float en3=sampleHG[3];
    float en4=sampleHG[4];
    float en6=sampleHG[6];
    if( en1<en3 && en3>en6 && (en4>en6||en2>en6) && en3>m_timeSample3ADCCut){
      HGCalTBRecHit recHit(rawhit.detid(), energy, lowGain, highGain, totGain, _time);
      if( rawhit.isUnderSaturationForHighGain() ) recHit.setUnderSaturationForHighGain();
      if( rawhit.isUnderSaturationForLowGain() ) recHit.setUnderSaturationForLowGain();
      HGCalTBDetId detid = rawhit.detid();
      CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(detid.layer(), detid.sensorIU(), detid.sensorIV(), detid.iu(), detid.iv(), SENSORSIZE );
      double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.first - HGCAL_TB_GEOMETRY::DELTA);
      double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.second - HGCAL_TB_GEOMETRY::DELTA);
      recHit.setCellCenterCoordinate(iux, iuy);
        
      HGCalTBLayer layer= essource_.layout_.at(rawhit.detid().layer()-1);
      int moduleId= layer.at( recHit.id().sensorIU(),recHit.id().sensorIV() ).moduleID();
      iski = rawhit.skiroc()%4;
      ASIC_ADC_Conversions adcConv=essource_.adccalibmap_.getASICConversions(moduleId,iski);

      if( rawhit.lowGainADC(3) > adcConv.TOT_lowGain_transition() ){
        energy = totGain * adcConv.TOT_to_lowGain() * adcConv.lowGain_to_highGain();
        recHit.setEnergyTOT(totGain);
        recHit.setFlag(HGCalTBRecHit::kLowGainSaturated);
        recHit.setFlag(HGCalTBRecHit::kGood);
      } else {
        
        fitter.run(sampleT, sampleLG, fitresultLG);
        fitter.run(sampleT, sampleHG, fitresultHG);

        if (fitresultLG.status==0) {
          m_h_LowGainVsTOTAmpl[10*iboard + iski%4]->Fill(totGain, fitresultLG.amplitude);   
          if (fitresultHG.status==0)
            m_h_HighVsLowGainAmpl[10*iboard + iski%4]->Fill(fitresultLG.amplitude, fitresultHG.amplitude);   
          
        }

        if( rawhit.highGainADC(3) > adcConv.lowGain_highGain_transition() ){
          recHit.setFlag(HGCalTBRecHit::kHighGainSaturated);
          if( fitresultLG.status==0 ){
            recHit.setEnergyLow(fitresultLG.amplitude);
            energy = fitresultLG.amplitude * adcConv.lowGain_to_highGain();
            _time = fitresultLG.tmax - fitresultLG.trise;
            recHit.setFlag(HGCalTBRecHit::kGood);

            LG_max_for_tree = fitresultLG.amplitude;
            TMax_LG_for_tree = fitresultLG.tmax;
            chi2_LG_for_tree = fitresultLG.chi2;
            trise_LG_for_tree = fitresultLG.trise;
            errortmax_LG_for_tree = fitresultLG.errortmax;
            erroramplitude_LG_for_tree = fitresultLG.erroramplitude;

            if (investigatePulseShape) {
              for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++) {
                shapesLG[key]->Fill(25*it+12.5-(fitresultLG.tmax - fitresultLG.trise), sampleLG[it]/fitresultLG.amplitude);
              }
              ToARisevsTMaxLG[key]->Fill(fitresultLG.tmax, toaRise);
              ToAFallvsTMaxLG[key]->Fill(fitresultLG.tmax, toaFall);
            }
          }
        } else{
          if( fitresultHG.status==0 ){
            recHit.setEnergyHigh(fitresultHG.amplitude);
            energy = fitresultHG.amplitude;
            _time = fitresultHG.tmax - fitresultHG.trise;
            recHit.setFlag(HGCalTBRecHit::kGood);


            HG_max_for_tree = fitresultHG.amplitude;
            TMax_HG_for_tree = fitresultHG.tmax;
            chi2_HG_for_tree = fitresultHG.chi2;
            trise_HG_for_tree = fitresultHG.trise;
            errortmax_HG_for_tree = fitresultHG.errortmax;
            erroramplitude_HG_for_tree = fitresultHG.erroramplitude;

            if (investigatePulseShape) {
              for( int it=0; it<NUMBER_OF_TIME_SAMPLES; it++) {
                shapesHG[key]->Fill(25*it+12.5-(fitresultHG.tmax - fitresultHG.trise), sampleHG[it]/fitresultHG.amplitude);
              }
              ToARisevsTMaxHG[key]->Fill(fitresultHG.tmax, toaRise);
              ToAFallvsTMaxHG[key]->Fill(fitresultHG.tmax, toaFall);
            }
          }
        }
      }

      //std::cout<<"Setting energy and time of: "<<detid.layer()<<"  "<<detid.iu()+7<<"  "<<detid.iv()+7<<"  "<<energy*adcConv.adc_to_MIP()<<"  "<<_time<<std::endl;
      recHit.setEnergy(energy*adcConv.adc_to_MIP());
      MIPEnergy_for_tree = energy*adcConv.adc_to_MIP();
      TOA_rise_for_tree = toaRise;
      TOA_fall_for_tree = toaFall;
      TOT_for_tree = totGain;   //is TOT slow
      //recHit.setTime(_time);
      recHit.setTime(toaRise);


      outtree->Fill();
      rechits->push_back(recHit);
    }

  }
  //event.put(cmMap, m_CommonModeNoiseCollectionName);
  event.put(rechits, m_outputCollectionName);
  #ifdef DEBUG
    eventCounter++;
  #endif
}

DEFINE_FWK_MODULE(HGCalTBRecHitProducer);
