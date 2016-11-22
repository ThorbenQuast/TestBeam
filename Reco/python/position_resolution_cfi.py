import FWCore.ParameterSet.Config as cms

position_resolution_analyzer = cms.EDAnalyzer("Position_Resolution_Analyzer",
                                weightingMethod = cms.string('squaredWeighting'),
                                fittingMethod = cms.string('lineTGraphErrors'),
                                pedestalThreshold = cms.double(30.),   #-99999 for no threshold-->subtracts the average
                                Layer_Z_Positions = cms.vdouble([1.2, 2., 3.5, 4.3, 5.8, 6.3, 8.7, 9.5, 11.4, 12.2, 13.8, 14.6, 16.6, 17.4, 20., 20.8]),
                                nLayers = cms.int32(8),
                                SensorSize = cms.int32(128),
                                EventsFor2DGraphs = cms.vint32([1]),
                                RUNDATA = cms.InputTag("source","RunData","unpack" ), 
                                HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" )#,
                              )