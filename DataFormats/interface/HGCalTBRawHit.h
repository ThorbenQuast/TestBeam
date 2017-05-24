#ifndef HGCALTBRAWHIT_H_INCLUDED
#define HGCALTBRAWHIT_H_INCLUDED 1

#include <vector>
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"

const static int NUMBER_OF_TIME_SAMPLES = 11;

class HGCalTBRawHit
{
 public:

  HGCalTBRawHit( ){;}
 HGCalTBRawHit( unsigned int rawid, std::vector<float> &adcHigh, std::vector<float> &adcLow) :
  m_rawid(rawid),
    m_adcHigh(adcHigh),
    m_adcLow(adcLow)
    {;}

  HGCalTBDetId detid() const {return HGCalTBDetId(m_rawid);}
  void setRawId(unsigned int rawid){m_rawid=rawid;}
  void setHighGainADCs(std::vector<float> vec){m_adcHigh=vec;}
  void setLowGainADCs(std::vector<float> vec){m_adcLow=vec;}
  float highGainADC(int timeSample){return m_adcHigh.at(timeSample);}
  float lowGainADC(int timeSample){return m_adcHigh.at(timeSample);}

 private:
  unsigned int m_rawid; //for some reasons (I don't know) root does not allow saving HGCalTBDetId because of some dictionary issue. 
  std::vector<float> m_adcHigh;
  std::vector<float> m_adcLow;
};

std::ostream& operator<<(std::ostream&, HGCalTBRawHit&);

#endif