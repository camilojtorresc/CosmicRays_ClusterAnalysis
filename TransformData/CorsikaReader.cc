/* Tanguy Pierog and Ralf Ulrich  18.12.2020
This is a simple template to read a CORSIKA DAT File and print part of the input
using COAST libraries. It is meant as a simple example that user should adapt
to their use. Most of the variables have been indicated but more can be found
in COAST doxygen documentation (in particular for the missing (here) longitudinal
profiles, Cherenkov or muon information)*/

/* Modified by Oliver Ruiz (Puebla, Mexico) October 2023
This script read the binary file and write two ASCII
files, one with the output of the particles information and
another with the summary of primary particle information
*/

#include <crsRead/MCorsikaReader.h>

#include <crs/TSubBlock.h>
#include <crs/MRunHeader.h>
#include <crs/MEventHeader.h>
#include <crs/MEventEnd.h>
#include <crs/MParticleBlock.h>
#include <crs/MLongitudinalBlock.h>
#include <crs/MParticle.h>

#include <TFile.h>

#include <iostream>
#include <sstream>
#include <map>
using namespace std;


// to hold data of one observation level
struct ObsLevel {
  double x;
  double y;
  double x2;
  double y2;
  double w;
};
          



int
main (int argc, char **argv) 
{  
  if (argc<2) {
    cout << "  please specify Corsika file " << endl;
    return 1;
  }
    
  string fname(argv[1]);
  crsRead::MCorsikaReader cr(fname, 3);
  
  string inputFile;
  if (fname.rfind('/')==string::npos) {
    inputFile = fname;
  } else {
    inputFile = fname.substr(fname.rfind('/')+1, 
                             fname.size()-fname.rfind('/'));
  }
  
  
  int ShowerCounter = 0;
  
  crs::MRunHeader Run;
  while (cr.GetRun (Run)) {

    ofstream MyFile(inputFile + "_particles.txt");                              

          MyFile << "sh\tid\tx\ty\tt\tpx\tpy\tpz\tPsq\tek\tzha\taza\n";

    ofstream SummaryFile(inputFile + "_showers.txt");
          
          SummaryFile << "Shower\tEnergy\tZfirst\tTheta\tPhi\tParticles\n";
          
    crs::MEventHeader Shower;
    while (cr.GetShower(Shower)) {
      
      ++ShowerCounter;
      
      const int nObsLevel = Shower.GetNObservationLevels();
      map<int, ObsLevel> obsLevel;
      
      for (int iObsLevel=1; iObsLevel<=nObsLevel; ++iObsLevel) { 
        
        double height = Shower.GetObservationHeight(iObsLevel-1);            
        ObsLevel emptyLevel;
        ostringstream tTitle, tName;
        tTitle << "Data at level " << iObsLevel;
        tName << "data_" << iObsLevel;
        
          emptyLevel.x  = 0;
          emptyLevel.y  = 0;
          emptyLevel.w  = 0;
          emptyLevel.x2 = 0;
          emptyLevel.y2 = 0;
        
        obsLevel[iObsLevel] = emptyLevel;

      } // end loop observation levels

      const double Energy = Shower.GetEnergy();
      const double zenith = Shower.GetTheta();
      const double azimuth = Shower.GetPhi();
      const double Zfirst = Shower.GetZFirst();
     
      crs::TSubBlock Data;
      while (cr.GetData (Data)) {

        switch (Data.GetBlockType ()) {
          
            case crs::TSubBlock::ePARTDATA:
            {

              const crs::MParticleBlock& ParticleData = Data;
              crs::MParticleBlock::ParticleListConstIterator iEntry;
              for (iEntry = ParticleData.FirstParticle();
                   iEntry != ParticleData.LastParticle();
                   ++iEntry) {

		// DUMP
		//iEntry->Dump();

                if (iEntry->IsParticle()) {
                  
                  crs::MParticle iPart(*iEntry);
                  
                  const int id    = iPart.GetParticleID();
                  const int level = iPart.GetObservationLevel();
                  const double w  = iPart.GetWeight();
                  const double ek  = iPart.GetKinEnergy();
                  const double px = iPart.GetPx();
                  const double py = iPart.GetPy();
                  const double pz = iPart.GetPz();
                  const double x  = iPart.GetX();
                  const double y  = iPart.GetY();
                  const double t  = iPart.GetTime();
                  const double Psq = iPart.GetPSquared();
                  const double zha = iPart.GetTheta();
                  const double aza = atan2(py,px);

                  /* Full list of variables available : 
CREAL 	GetWeight () const 
virtual std::string GetParticleName () const 
int 	GetParticleID () const 
int 	GetParticleId () const 
bool 	IsParticle () const 
bool 	IsNucleus () const 
bool 	IsCherenkov () const 
bool 	IsMuonProductionInfo () const 
bool 	IsEmpty () const 
ParticleType 	GetType () const 
virtual int 	GetObservationLevel () const 
virtual int 	GetHadronicGeneration () const 
virtual CREAL 	GetPx () const 
virtual CREAL 	GetPy () const 
virtual CREAL 	GetPz () const 
virtual CREAL 	GetX () const 
virtual CREAL 	GetY () const 
virtual CREAL 	GetTime () const 
double 	GetMass () const 	mass in GeV
int 	GetPDGCode () const 
double 	GetKinEnergy () const 	kin. energy in GeV
double 	GetTheta () const 	zenith angle in rad
double 	GetPSquared () const 	squared of momentum in GeV
double 	GetE () const 
                  */
                
		  if (obsLevel.count(level)==0) {
		    cout << " detected new obs-level " << level 
			 << ". Possibly INCLIN-Option " << endl;
		    ObsLevel emptyLevel;
		    ostringstream tTitle, tName;
		    tTitle << "Data at level " << level;
		    tName << "data_" << level;
		    obsLevel[level] = emptyLevel;
		  }
		  
                  obsLevel[level].x  += x*w;
                  obsLevel[level].y  += y*w;
                  obsLevel[level].x2 += x*x*w;
                  obsLevel[level].y2 += y*y*w;
                  obsLevel[level].w  += w;
                  

          MyFile << ShowerCounter
          << "\t" << id 
          << "\t" << x 
          << "\t" << y 
          << "\t" << t 
          << "\t" << px
          << "\t" << py 
          << "\t" << pz
          << "\t" << Psq
          << "\t" << ek
          << "\t" << zha
          << "\t" << aza
          << "\n";

                }

              } // end particle loop

              break;
            }
            
            case crs::TSubBlock::eLONG:
              break;
              
            default:
              break;
        } // end data block

      } // loop data


      for (map<int, ObsLevel>::iterator iLevel = obsLevel.begin();
           iLevel != obsLevel.end();
           ++iLevel) {

        double npart=iLevel->second.w;
        if(npart>0)
        SummaryFile << ShowerCounter 
        << "\t" << Energy
        << "\t" << Zfirst           
        << "\t" << zenith
        << "\t" << azimuth        
        << "\t" << npart << "\n";
        
      } // loop observation levels
      
      
    } // loop shower

    SummaryFile.close();
    
    MyFile.close();                                                                  

 } // loop runs (usually just 1)    
  
  cout << " Read " << ShowerCounter << " showers from file " << endl;
  
  return 0;
}


  
