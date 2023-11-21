// mclust_7.C    6.junio.2022
// Version _6.C + Incluyo un nuevo fichero de salida con resumen de particulas individuales
// Version _5.C. Incluyo energia de protones, Energia de electrones aislados, dif. de tiempos electrones-muones, densidad de energia y un fichero con la altura media de produccion
// 23.4.22: Redefino algunos intervalos
// Nueva version de _3.C incluyendo distribucion radial de tiempos de llegada
// Nueva version del _2.C con nuevas variables de la salida para NNetwork analysis
// Atencion!!! Cambio el orden de algunas variables en algun fichero de salida !!!
// Nueva version mclust_ce con nuevas variables en el output file de showers
// Version con estimacion de la energia de los clusters de electrones y con informacion de la altura de primera interaccion (mclust_he.C)
// ... with many, many, many improvements !!!!  Feb14.22
//    Cambios relevantes en la informacion de los clusters y en los resumenes
// El programa graba ficheros con la siguiente informacion:
// ========================
// *** OutputFiles con informacion detallada. Se recomienda un archivo por run.
//of_clushc.txt: Informacion detallada de cada cluster.
//of_clushs.txt: Informacion detallada de cada shower.
// *** OutputFiles con informacion resumida. Se recomienda grabar la cabecera en la primera pasada y despues comentarla para producir una unica tabla de datos.
//of_clushp.txt: Informacion muy resumida del numero de particulas por run
//of_clushn.txt: Informacion seleccionada para su uso en analisis mediante NeuralNetworks
//of_clushr.txt: Informacion radial de clusters EM, MU y MX
// =========================
#define mclust_cxx
#include "mclust.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <fstream>
#include <stdlib.h>
// Secondary masses
#define mele   0.000511 // GeV
#define mmu    0.10566  // GeV
#define mgam   0.0      // GeV
#define mp     0.9383   // GeV
#define mn     0.9396   // GeV
#define mother 0.14     // GeV. we assume they are mainly pions
// Cluster size analyzed
#define dsize  1.8   // m^2. Size of the detector
//#define mxdst 1.5 //m. Distancia maxima entre partículas 2R = 2 * sqrt(1.8/pi)
//#define sigrmx 1.1 //m. Dispersion radial maxima del cluster ~mxdst/0.7
// ------------------------------------------------------------------------------
// ------------------------------------------------------------- Printing options
#define iprhd 1     // Index for printing headers. 0=No Print.  1=Print
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
#define fout0 "xmdatp_p_E30b_2m2.txt" // Particle output
#define fout1 "xmdatc_p_E30b_2m2.txt" // Clusters output
#define fout2 "xmdats_p_E30b_2m2.txt" // Showers output
#define fout3 "xmdatm_p_E30b_2m2.txt" // Main Summary
#define fout4 "xmdatr_p_E30b_2m2.txt" // Cluster Radial distribution
#define fout5 "xmdath_p_E30b_2m2.txt" // Mean Height estimation
// -----------------------------
// Primary cosmic ray properties
#define mpcr 1       // A mass of primary cosmic ray: 1, 4, 12, 56
// f0/m2·s·sr·GeV/nucleon
#define f0 1.8E4     // Intercepts H, He, C, Fe: 1.8E4, 5.0E3, 1.0E3, 1.5E3
#define spindex 2.7  // Spectral index
#define zhmin 0     // Zenith angle min 0, 25, 37, 46, 53
#define zhmax 25    // Zenith angle max
#define dene 0.25    // Width of the LogEnergy interval generated
// ------------------------------
//  H:  8 Height Intervals: 0, 12, 16, 18, 20, 22, 28, 40, 150 km
//  He: 8 Height Intervals: 0, 16, 20, 22, 24, 26, 30, 36, 150 km
//  C:  8 Height Intervals: 0, 24, 30, 36, 40, 44, 48, 64, 150 km
//  Fe: 8 Height Intervals: 0, 32, 40, 44, 48, 52, 56, 64, 150 km
// ------------------------------
#define hmin 0.
#define hmax 150.
//
Int_t nshana = 100000; // Max. N.Showers to analize
//
Int_t mxshop = 50;   // nshana/50; // 100; // Mx. N.Showers to print
Int_t mxclup = 50;   // 10;  // Mx. N.Clusters/Shw to print
Int_t mxparp = 500;  // 10; // Mx. N.Particles/Energy to print
Int_t eperiod = 0, enemod=0;     // shower period and module to print particles
Int_t itsana=0, ishowp=0, iclusp=0; // iCounters
// ------------------------------
// Data Files
// Simulation with 4 samples
// nh: n alturas. an: num. de intervalo angular
Int_t nsamp = 1;     // Number of samples per shower.
Int_t isami = 1;     // N. Initial Sample
Float_t gm1  = spindex - 1;  // pcr energy distribution decrease factor
//
void mclust::Loop()
{
// Analisis de cascadas atmosfericas de rayos cosmicos
/*
 *  Copyright LabCAF. IGFAE/USC. All rights reserved.
 *  Created by Juan A. Garzon on 23/01/13
 *  Modified by G.Kornakov on  11/03/2013
 *  Modified by Yanis Fontenla Barba on 23/05/2017 for Corsika
 *  Modified by JAGarzon on Apr.2021 for cluster analysis
*/
    
if (fChain == 0) return;
Long64_t nshows = fChain->GetEntriesFast();     // no. of showers in file

    if(iprhd == 1){cout <<  endl << endl << "***   Print Header"   << endl << endl << endl;}
    
cout << "***************************** " << endl;
cout << "- nShowers in file:  " << nshows << endl;
if(nshana > nshows) {nshana = nshows;}
cout << "- nShowers analyzed: " << nshana << endl;
Float_t denes = dene/nshana;  // ~ Energy Interval Width associated to each Shower
cout << "***************************** " << endl;
    
eperiod = nshana/ mxparp;
Long64_t nbytes = 0, nb = 0;
Long64_t fbytes = 0, fb = 0;
Long64_t ebytes = 0, eb = 0;
Long64_t hbytes = 0, hb = 0;
Long64_t tbytes = 0, tb = 0;
    
// ===========================================
// Identificadores de particulas en Corsika:
//  1        gammas
//  2 3      e- e+
//  5 6      mu- mu+
//  13 14    n p
// ============================================
Float_t rad2g;
rad2g = 180/TMath::Pi();
// -----------------------------
Float_t mxdst, sigrmx;
mxdst  = 2 * sqrt(dsize/TMath::Pi());  // diameter of a dsize area
sigrmx = mxdst / sqrt(2);
// -------------------------------
Int_t tag;
Float_t gpcut = 0.1;    // gammas minimum momentum cut
Int_t   zero=0, iconta=0, icontb=0, icontc=0; // auxiliar counters for cout's
Int_t   ncont=0, n=0, nsecs, isep=0, ntsep=0; // ntsec;
Long64_t icgam=1, icele=1000, icmu=100000, icpt=10000000, icnt=100000000;  // codes
Int_t   pid, pidf, pidfr1=-1, pidfr2=-1, pidfr3=-1, pidfr4=-1, pidfr5=-1, pidfr6=-1;
Long64_t pic=0, pic1, pic2, picfp, piclp;  // id code of particles, first, last...
Long64_t clid,  sclid;  // ClusterId, ShorClusterId
Int_t   ngamt=0, nelt=0, nmut=0, nnt=0, npt=0, nopt=0;  // N. of particles
Int_t   ngams=0, nels=0, nmus=0, nns=0, nps=0, nops=0;  // N. of particles in shower
Int_t   ngamc=0, nelec=0, nmuc=0;  // N. of particles in cluster
Int_t   ncluss, nclust=0, cmult, cmulp1,
        nelcls=0, nmucls=0, nmxcls=0, negcls=0, ngmcls=0, notcls=0, // cluster count in showers
        nelclt=0, nmuclt=0, nmxclt=0, negclt=0, ngmclt=0, notclt=0; // total cluster count
//
Float_t pcrene, emin=1.E20, emax=0., lemin, lemax, flint, flints, lflint, lflins; // ELimits,
Float_t ene, gmten=0, elten=0;
Float_t elment=0, mument=0, elmens=0, mumens=0, emrats=0; // Mean Energies
Float_t rx, ry, rr, dxc, dyc, dtc, edc, x1, x2, y1, y2, drc, rclmn;
Float_t xfpc=0, xlpc=0, yfpc=0, ylpc=0, drcsq=0, cltdt=0, clang=0, cltst=0, cltsz=0;
    
Float_t eclmlt=0, eclen=0, ecldt=0, eclth=0, eclsz=0, ecled=0, eclem=0, eclet=0, eclest=0;
Float_t egcmlt=0, egcen=0, egcdt=0, egcth=0, egcsz=0, egced=0, egcem=0, egcet=0, egcest=0;
Float_t gmcmlt=0, gmcen=0, gmcdt=0, gmcth=0, gmcsz=0, gmced=0, gmcem=0, gmcet=0, gmcest=0;
    
Float_t dt, dtsq, dxsq, dysq,  sigrc;
Float_t pxelc=0, pyelc=0, pzelc=0, pxegc=0, pyegc=0, pzegc=0, pxgmc=0, pygmc=0, pzgmc=0;
Float_t rmgms=0, rmgmt=0, rmels=0, rmelt=0, rmmus=0, rmmut=0,
        rmnts=0, rmntt=0, rmpts=0, rmptt=0, rmots=0, rmott=0;
Float_t reclsq=0, recl=0, rmecls=0, rmeclt=0, rmclsq=0, rmcl=0, rmmcls=0, rmmclt=0,
    rxclsq=0, rxcl=0, rmxcls=0, rmxclt=0;
Float_t regcsq=0, regc=0, rmegcs=0, rmegct=0, rgclsq=0, rgcl=0, rmgcls=0, rmgclt=0, roclsq=0, rocl=0, rmocls=0, rmoclt=0;
Float_t zha=0, zen1=0, zen2=0, aza=0, azh1=0, azh2=0;
Float_t t0, time =0, t1=0, t2=0, tfg=0, tlg=0, tfe=0, tle=0, tfm=0, tlm=0,
        tfpc=0, tlpc = 0;
Float_t zhfp=0, zhlp=0, azfp=0, azlp=0,
        zhfe=0, zhle=0, zhfm=0, zhlm=0, azfe=0, azle=0, azfm=0, azlm=0;
Float_t px=0, py=0, pz=0, px1=0, py1=0, pz1=0, px2=0, py2=0, pz2=0,
        pmod, pm1, pm2, psq, pfpc, plpc; // particle momenta
Float_t xclmn,   yclmn,   tclmn, xcmno,   ycmno,   tcmno;
Float_t sigx=0, sigxsq=0, sigy=0, sigysq=0, sigtc=0, sigtsq=0;
Float_t ssx=0, ssy=0, sst=0, ssxo, ssyo, ssto; // auxiliar variable for variance
Float_t nxfp, nyfp, nzfp, nxlp, nylp, nzlp;
Float_t finhgt, mnhgt=0;  // first height, mean hght
Float_t lpcren, lpcrem=0; // E0 energy, energy interval, mean energy
//
Int_t jr, jt, ja; // indexes for distances, times and angles
// Some counter intervals
Int_t   ir1   = 100,  ir2 = 200,    ir3 = 500,    ir4 = 1000,   ir5 = 2500;   // Radial Rings
Float_t igmen1= 0.15, igmen2 =0.25, igmen3 = 0.4, igmen4 = 0.6, igmen5 = 1.;  // GamEnergy
Float_t ielen1 = 0.1, ielen2 = 0.2, ielen3 = 0.3, ielen4 = 0.6, ielen5 = 1.;  // EleEnergy
Int_t   imuen1 = 2,   imuen2 = 4,   imuen3 = 6,   imuen4 = 10,  imuen5 = 20;  // MuEnergy
Int_t   inten1 = 2,   inten2 = 3,   inten3 = 5,   inten4 = 8,   inten5 = 15;  // NeutEnergy
Int_t   ipren1 = 2,   ipren2 = 3,   ipren3 = 5,   ipren4 = 8,   ipren5 = 15;  // ProtEnergy
Int_t   ieclm1 = 1,   ieclm2 = 2,   ieclm3 = 3,   ieclm4 = 4,   ieclm5 = 5;   // ECl.Mult.
Float_t iecen1 = 0.1, iecen2 = 0.2, iecen3 = 0.3, iecen4 = 0.6, iecen5 = 1.;  // ECl Energy
Float_t iecdt1 = 0.3, iecdt2 = 0.6, iecdt3 = 1.0, iecdt4 = 1.5, iecdt5 = 2.5; // ECl dTime
Int_t   iecth1 = 5,   iecth2 = 10,  iecth3 = 15,  iecth4 = 20,  iecth5 = 25;  // ECl.Therm.
Int_t   iecan1 = 2,   iecan2 = 4,   iecan3 = 8,   iecan4 = 20,  iecan5 = 40;  // ECl AnglInt.
Float_t iecsz1 = 0.2, iecsz2 = 0.4, iecsz3 = 0.6, iecsz4 = 0.8, iecsz5 = 1.0; // ECl Size
Float_t ieced1 = 0.5, ieced2 = 1,   ieced3 = 2,   ieced4 = 5,   ieced5 = 20;  // ECl EnDens
Float_t iecem1 = 0.1, iecem2 = 0.2, iecem3 = 0.3, iecem4 = 0.4, iecem5 = 0.8; // ECl En/Mult
Float_t iecet1 = 0.4, iecet2 = 0.8, iecet3 = 1.5, iecet4 = 3,   iecet5 = 6;   // ECl En/dTime
Float_t ieest1 = 0.2, ieest2 = 0.4, ieest3 = 0.8, ieest4 = 2,   ieest5 = 5;   // ECl EnDens/dT
Float_t iegen1 = 0.3, iegen2 = 0.4, iegen3 = 0.6, iegen4 = 1,   iegen5 = 1.5; // EGCl. Energy
Float_t igcen1 = 0.3, igcen2 = 0.4, igcen3 = 0.6, igcen4 = 1,   igcen5 = 1.5; // GmCl Energy

// Counters
Int_t   ngmen1=0, ngmen2=0, ngmen3=0, ngmen4=0, ngmen5=0, ngmen6=0; // GamEnergy Count.
Int_t   ngmes1=0, ngmes2=0, ngmes3=0, ngmes4=0, ngmes5=0, ngmes6=0; // GamEnergy CountShw
Int_t   nelen1=0, nelen2=0, nelen3=0, nelen4=0, nelen5=0, nelen6=0; // EleEnergy Count.
Int_t   neles1=0, neles2=0, neles3=0, neles4=0, neles5=0, neles6=0; // EleEnergy CountShw
Int_t   nmuen1=0, nmuen2=0, nmuen3=0, nmuen4=0, nmuen5=0, nmuen6=0; // MuEnergy Count
Int_t   nmues1=0, nmues2=0, nmues3=0, nmues4=0, nmues5=0, nmues6=0; // MuEnergy CountShow
Int_t   nnten1=0, nnten2=0, nnten3=0, nnten4=0, nnten5=0, nnten6=0; // NeutEnergy Count.
Int_t   nntes1=0, nntes2=0, nntes3=0, nntes4=0, nntes5=0, nntes6=0; // NeutEnergy CountShw
Int_t   npren1=0, npren2=0, npren3=0, npren4=0, npren5=0, npren6=0; // ProtEnergy Count.
Int_t   npres1=0, npres2=0, npres3=0, npres4=0, npres5=0, npres6=0; // ProtEnergy CountShw
Int_t   neclm1=0, neclm2=0, neclm3=0, neclm4=0, neclm5=0, neclm6=0; // EClMult. Count
Int_t   necms1=0, necms2=0, necms3=0, necms4=0, necms5=0, necms6=0; // EClMult. CountShow
Int_t   necen1=0, necen2=0, necen3=0, necen4=0, necen5=0, necen6=0; // EClEnergy Count
Int_t   neces1=0, neces2=0, neces3=0, neces4=0, neces5=0, neces6=0; // EClEnergy CountShw
Int_t   necdt1=0, necdt2=0, necdt3=0, necdt4=0, necdt5=0, necdt6=0; // ECldT Counters
Int_t   nedts1=0, nedts2=0, nedts3=0, nedts4=0, nedts5=0, nedts6=0; // ECldT COuntShow
Int_t   necth1=0, necth2=0, necth3=0, necth4=0, necth5=0, necth6=0; // EClTherm Count
Int_t   neths1=0, neths2=0, neths3=0, neths4=0, neths5=0, neths6=0; // EClTherm CountSh
Int_t   necsz1=0, necsz2=0, necsz3=0, necsz4=0, necsz5=0, necsz6=0; // EClSize Count
Int_t   neszs1=0, neszs2=0, neszs3=0, neszs4=0, neszs5=0, neszs6=0; // EClSize CountShow
Int_t   neced1=0, neced2=0, neced3=0, neced4=0, neced5=0, neced6=0; // EClEnDens Count.
Int_t   needs1=0, needs2=0, needs3=0, needs4=0, needs5=0, needs6=0; // EClEnDens COuntShow
Int_t   necem1=0, necem2=0, necem3=0, necem4=0, necem5=0, necem6=0; // EClEn/Mult Count.
Int_t   neems1=0, neems2=0, neems3=0, neems4=0, neems5=0, neems6=0; // EClEn/Mult COuntShow
Int_t   necet1=0, necet2=0, necet3=0, necet4=0, necet5=0, necet6=0; // EClEn/dT Count.
Int_t   neets1=0, neets2=0, neets3=0, neets4=0, neets5=0, neets6=0; // EClEn/dT COuntShow
Int_t   neest1=0, neest2=0, neest3=0, neest4=0, neest5=0, neest6=0; // EClEn/Sz/dT Count.
Int_t   nests1=0, nests2=0, nests3=0, nests4=0, nests5=0, nests6=0; // EClEn/dT COuntShow
Int_t   negen1=0, negen2=0, negen3=0, negen4=0, negen5=0, negen6=0; // EGClEnergy Count.
Int_t   neges1=0, neges2=0, neges3=0, neges4=0, neges5=0, neges6=0; // EGClEnergy Count.
Int_t   ngcen1=0, ngcen2=0, ngcen3=0, ngcen4=0, ngcen5=0, ngcen6=0; // EGClEnergy Count.
Int_t   ngces1=0, ngces2=0, ngces3=0, ngces4=0, ngces5=0, ngces6=0; // EGClEnergy Count.
// Mean Heigth
Float_t hgmen1=0, hgmen2=0, hgmen3=0, hgmen4=0, hgmen5=0, hgmen6=0; // GamEnerg WCount
Float_t helen1=0, helen2=0, helen3=0, helen4=0, helen5=0, helen6=0; // EleEnerg WCount
Float_t heclm1=0, heclm2=0, heclm3=0, heclm4=0, heclm5=0, heclm6=0; // EClMult. WCount
Float_t hecen1=0, hecen2=0, hecen3=0, hecen4=0, hecen5=0, hecen6=0; // EClEnerg WCount
Float_t hecdt1=0, hecdt2=0, hecdt3=0, hecdt4=0, hecdt5=0, hecdt6=0; // ECldT WCount
Float_t hecth1=0, hecth2=0, hecth3=0, hecth4=0, hecth5=0, hecth6=0; // EClTherm WCount
Float_t hecsz1=0, hecsz2=0, hecsz3=0, hecsz4=0, hecsz5=0, hecsz6=0; // EClSize WCount
Float_t heced1=0, heced2=0, heced3=0, heced4=0, heced5=0, heced6=0; // EClEnDens WCount
Float_t hecem1=0, hecem2=0, hecem3=0, hecem4=0, hecem5=0, hecem6=0; // EClEn/Mult WCount
Float_t hecet1=0, hecet2=0, hecet3=0, hecet4=0, hecet5=0, hecet6=0; // EClEn/dTime WCount
Float_t heets1=0, heets2=0, heets3=0, heets4=0, heets5=0, heets6=0; // EClEn/dTime WCount
Float_t hegen1=0, hegen2=0, hegen3=0, hegen4=0, hegen5=0, hegen6=0; // EGClEnergy WCount
Float_t hgcen1=0, hgcen2=0, hgcen3=0, hgcen4=0, hgcen5=0, hgcen6=0; // GmClEnergy WCount
Float_t hmuen1=0, hmuen2=0, hmuen3=0, hmuen4=0, hmuen5=0, hmuen6=0; // MuEnergy WCount
Float_t hnten1=0, hnten2=0, hnten3=0, hnten4=0, hnten5=0, hnten6=0; // NeutEnerg WCount
Float_t hpren1=0, hpren2=0, hpren3=0, hpren4=0, hpren5=0, hpren6=0; // ProtEnerg WCount
Float_t mhgme1=0, mhgme2=0, mhgme3=0, mhgme4=0, mhgme5=0, mhgme6=0; // EleEnerg
Float_t mhele1=0, mhele2=0, mhele3=0, mhele4=0, mhele5=0, mhele6=0; // EleEnerg
Float_t mhecm1=0, mhecm2=0, mhecm3=0, mhecm4=0, mhecm5=0, mhecm6=0; // EClMult.
Float_t mhece1=0, mhece2=0, mhece3=0, mhece4=0, mhece5=0, mhece6=0; // EClEnerg
Float_t mhedt1=0, mhedt2=0, mhedt3=0, mhedt4=0, mhedt5=0, mhedt6=0; // ECldT
Float_t mheth1=0, mheth2=0, mheth3=0, mheth4=0, mheth5=0, mheth6=0; // EClTherm
Float_t mhcsz1=0, mhcsz2=0, mhcsz3=0, mhcsz4=0, mhcsz5=0, mhcsz6=0; // EClSize
Float_t mheed1=0, mheed2=0, mheed3=0, mheed4=0, mheed5=0, mheed6=0; // EClEnDens
Float_t mheem1=0, mheem2=0, mheem3=0, mheem4=0, mheem5=0, mheem6=0; // EClEnDens
Float_t mheet1=0, mheet2=0, mheet3=0, mheet4=0, mheet5=0, mheet6=0; // EClEnDens
Float_t mhets1=0, mhets2=0, mhets3=0, mhets4=0, mhets5=0, mhets6=0; // EClEnDens
Float_t mhege1=0, mhege2=0, mhege3=0, mhege4=0, mhege5=0, mhege6=0; // EGClEnergy
Float_t mhgce1=0, mhgce2=0, mhgce3=0, mhgce4=0, mhgce5=0, mhgce6=0; // GmClEnergy
Float_t mhmue1=0, mhmue2=0, mhmue3=0, mhmue4=0, mhmue5=0, mhmue6=0; // MuEnergy
Float_t mhnte1=0, mhnte2=0, mhnte3=0, mhnte4=0, mhnte5=0, mhnte6=0; // NeutEnerg
Float_t mhpre1=0, mhpre2=0, mhpre3=0, mhpre4=0, mhpre5=0, mhpre6=0; // ProtEnerg
//Mean Distance (radius)
Float_t recen1=0, recen2=0, recen3=0, recen4=0, recen5=0, recen6=0; // EClEnerg RCount
Float_t mrgam=0,  mrele=0,  mrmuon=0, mrneut=0, mrprot=0, mroth=0; // Particles RMean
Float_t mrece1=0, mrece2=0, mrece3=0, mrece4=0, mrece5=0, mrece6=0; // EClEnerg RMean
// Mean dTime
Float_t tecen1=0, tecen2=0, tecen3=0, tecen4=0, tecen5=0, tecen6=0; // EClEnerg dTCount
Float_t mtece1=0, mtece2=0, mtece3=0, mtece4=0, mtece5=0, mtece6=0; // EClEnerg dTMean
// Relative rates for different selections
Float_t rrecm1, rrecm2, rrecm3, rrecm4, rrecm5, rrecm6; // EClMult.
Float_t rrece1, rrece2, rrece3, rrece4, rrece5, rrece6; // EClEnerg
Float_t rredt1, rredt2, rredt3, rredt4, rredt5, rredt6; // ECldT
Float_t rreth1, rreth2, rreth3, rreth4, rreth5, rreth6; // EClTherm
Float_t rrcsz1, rrcsz2, rrcsz3, rrcsz4, rrcsz5, rrcsz6; // EClSize
Float_t rreed1, rreed2, rreed3, rreed4, rreed5, rreed6; // EClEnDens
Float_t rrmue1, rrmue2, rrmue3, rrmue4, rrmue5, rrmue6; // MuEnergy
// Radial analysis
Float_t t0r1=0,   t0r2=0,   t0r3=0,   t0r4=0,   t0r5=0,   t0r6=0;   // Radial t0's
Int_t   nelr1 =0, nelr2 =0, nelr3 =0, nelr4 =0, nelr5 =0, nelr6 =0; // ElRadial Count
Int_t   nelrs1=0, nelrs2=0, nelrs3=0, nelrs4=0, nelrs5=0, nelrs6=0; // ElRadial CountShow
Int_t   nmur1 =0, nmur2 =0, nmur3 =0, nmur4 =0, nmur5 =0, nmur6 =0; // MuRadial Count
Int_t   nmurs1=0, nmurs2=0, nmurs3=0, nmurs4=0, nmurs5=0, nmurs6=0; // MuRadial CountShow
Int_t   neclr1=0, neclr2=0, neclr3=0, neclr4=0, neclr5=0, neclr6=0; // EClRadial Count
Int_t   necrs1=0, necrs2=0, necrs3=0, necrs4=0, necrs5=0, necrs6=0; // EClRadial CountShow
Float_t elmer1=0, elmer2=0, elmer3=0, elmer4=0, elmer5=0, elmer6=0; // El mEnergy in R slices
Float_t emers1=0, emers2=0, emers3=0, emers4=0, emers5=0, emers6=0; // El mEnergy in show
Float_t ecmer1=0, ecmer2=0, ecmer3=0, ecmer4=0, ecmer5=0, ecmer6=0; // ECl mEnergy
Float_t elmtr1=0, elmtr2=0, elmtr3=0, elmtr4=0, elmtr5=0, elmtr6=0; // El mTime in R slices
Float_t emtrs1=0, emtrs2=0, emtrs3=0, emtrs4=0, emtrs5=0, emtrs6=0; // El mTime in R slices
Float_t clmer1=0, clmer2=0, clmer3=0, clmer4=0, clmer5=0, clmer6=0; // ECl mTime in R slices
Float_t clmtr1=0, clmtr2=0, clmtr3=0, clmtr4=0, clmtr5=0, clmtr6=0; // ECl mTime in R slices
Float_t mumer1=0, mumer2=0, mumer3=0, mumer4=0, mumer5=0, mumer6=0; // Mu mEne in Rslices
Float_t mmers1=0, mmers2=0, mmers3=0, mmers4=0, mmers5=0, mmers6=0; // Mu mEne Rslices in Shws
Float_t mumtr1=0, mumtr2=0, mumtr3=0, mumtr4=0, mumtr5=0, mumtr6=0; // Mu mTime
Float_t mmtrs1=0, mmtrs2=0, mmtrs3=0, mmtrs4=0, mmtrs5=0, mmtrs6=0; // Mu mTime in Shows
//
Float_t emdtr1=0, emdtr2=0, emdtr3=0, emdtr4=0, emdtr5=0, emdtr6=0; // El&Mu dTime
Float_t emdts1=0, emdts2=0, emdts3=0, emdts4=0, emdts5=0, emdts6=0; // El&Mu dTim in Show
Float_t emrat1=0, emrat2=0, emrat3=0, emrat4=0, emrat5=0, emrat6=0; // El/Mu Rate
Float_t emrts1=0, emrts2=0, emrts3=0, emrts4=0, emrts5=0, emrts6=0; // El/Mu Rate in Show
// -------------------------------------------------------------
fstream file0; file0.open(fout0, fstream::app); //- Particle Summary
fstream file1; file1.open(fout1, fstream::app); //- Clusters Summary  (or fstream::out);
fstream file2; file2.open(fout2, fstream::app); //- Shower Summary
fstream file3; file3.open(fout3, fstream::app); //- Main Summary
fstream file4; file4.open(fout4, fstream::app); //- Radial distributions Summary
fstream file5; file5.open(fout5, fstream::app); //- Mean height estimation
// --------------------------------------------------------------
// -------------------------------------------------------------- Print Main headers
if(iprhd==1){
file0 << "# Particle Summary: DetMxD, DetMSg, PCRZhMin, PCRZhMx: " <<"\t"<<
        mxdst <<"\t"<< sigrmx <<"\t"<< zhmin <<"\t"<< zhmax <<"\t"<<
    " HMin, HMax: "<<"\t"<< hmin <<"\t"<< hmax << endl;
//    "Height intervals:  0, 12, 16, 18, 20, 22, 28, 40, 150" << endl;
//file0 << "---------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
file0 <<   // - Particle Summary header
    "#PrCR"  <<"\t"<< "LPCREn" <<"\t"<< "HghPCR" <<"\t"<<
    "IShow"  <<"\t"<<
    "PId"    <<"\t"<< "RPart"  <<"\t"<< "enert" <<"\t"<<
    "AzAPar" <<"\t"<< "ZhAPar" <<"\t"<<
    endl;
file1 << "# Clusters Summary: DetMxD, DetMSg, PCRZhMin, PCRZhMx: " <<"\t"<<
        mxdst <<"\t"<< sigrmx <<"\t"<< zhmin <<"\t"<< zhmax <<"\t" <<
    " HMin, HMax: "<<"\t"<< hmin <<"\t"<< hmax << endl;
//    "Height intervals:  0, 12, 16, 18, 20, 22, 28, 40, 150" << endl;
//file1 << "---------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
file1 <<   // - Clusters summary header
    "#PrCR"   <<"\t"<< "LPCREn"  <<"\t"<< "HghPCR" <<"\t"<<
    "IShow"   <<"\t"<< "iClust"  <<"\t"<<
    "ClsIdC"  <<"\t"<< "sClIdC"  <<"\t"<<  // ClId & Short ClId Code: 0, 1, 2, 3, 5, 6
    "NpartC"  <<"\t"<< "NGamC"   <<"\t"<< "NEleC"  <<"\t"<< "NMuCl"    <<"\t"<<
    "XmClst"  <<"\t"<< "YmClst"  <<"\t"<< "RmClst" <<"\t"<< "SigRCl"   <<"\t"<<
    "TmClst"  <<"\t"<< "dTClst"  <<"\t"<< "sTClst" <<"\t"<<
    "FstPID"  <<"\t"<< "FstPCX"  <<"\t"<< "FstPCY" <<"\t"<<  "FstPCT"  <<"\t"<<
    "FstPZh"  <<"\t"<< "FstPAz"  <<"\t"<< "FstPPm" <<"\t"<<
    "LstPID"  <<"\t"<< "LstPCX"  <<"\t"<< "LstPCY" <<"\t"<<  "LstPCT"  <<"\t"<<
    "LstPZh"  <<"\t"<< "LstPAz"  <<"\t"<< "LstPPm" <<"\t"<<
    "EClMlt"  <<"\t"<< "ECltEn"  <<"\t"<< "ECltdT" <<"\t"<<  "ECltTh" <<"\t"<<
    "ECltSz"  <<"\t"<< "ECltED"  <<"\t"<< "ECltEM" <<"\t"<<  "ECltET" <<"\t"<<
    "EClETS"  <<"\t"<<
    "EGmMlt"  <<"\t"<< "EGmCEn"  <<"\t"<< "EGmCdT" <<"\t"<<  "EGmCTh" <<"\t"<<
    "EGmCSz"  <<"\t"<< "EGmCED"  <<"\t"<< "EGmCEM" <<"\t"<<  "EGmCET" <<"\t"<<
    "EGmETS"  <<"\t"<<
    "GmCMlt"  <<"\t"<< "GmClEn"  <<"\t"<< "GmCldT" <<"\t"<<  "GmClTh" <<"\t"<<
    "GmClSz"  <<"\t"<< "GmClED"  <<"\t"<< "GmClEM" <<"\t"<<  "GmClET" <<"\t"<<
    "GmCETS"  <<"\t"<<
    endl;     // cluster openning angle
file2 << "# Showers Summary: DetMxD, DetMSg, PCRZhMin, PCRZhMx: " <<"\t"<<
        mxdst <<"\t"<< sigrmx <<"\t"<< zhmin <<"\t"<< zhmax <<"\t"<<
        " HMin, HMax: "<<"\t"<< hmin <<"\t"<< hmax << endl;
//        "Height intervals:  0, 12, 16, 18, 20, 22, 28, 40, 150" << endl;
file2 << "# \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t "
    << 0 <<"\t"<< igmen1 <<"\t"<< igmen2 <<"\t"<< igmen3 <<"\t"<< igmen4 <<"\t"<< igmen5 <<"\t"
    << 0 <<"\t"<< ielen1 <<"\t"<< ielen2 <<"\t"<< ielen3 <<"\t"<< ielen4 <<"\t"<< ielen5 <<"\t"
    << 0 <<"\t"<< imuen1 <<"\t"<< imuen2 <<"\t"<< imuen3 <<"\t"<< imuen4 <<"\t"<< imuen5 <<"\t"
    << ieclm1 <<"\t"<< ieclm2 <<"\t"<< ieclm3 <<"\t"<< ieclm4 <<"\t"<< ieclm5 <<"\t"<< ieclm5+1 <<"\t"
    << 0 <<"\t"<< iecen1 <<"\t"<< iecen2 <<"\t"<< iecen3 <<"\t"<< iecen4 <<"\t"<< iecen5 <<"\t"
    << 0 <<"\t"<< iecdt1 <<"\t"<< iecdt2 <<"\t"<< iecdt3 <<"\t"<< iecdt4 <<"\t"<< iecdt5 <<"\t"
    << 0 <<"\t"<< iecth1 <<"\t"<< iecth2 <<"\t"<< iecth3 <<"\t"<< iecth4 <<"\t"<< iecth5 <<"\t"
    << 0 <<"\t"<< iecsz1 <<"\t"<< iecsz2 <<"\t"<< iecsz3 <<"\t"<< iecsz4 <<"\t"<< iecsz5 <<"\t"
    << 0 <<"\t"<< iegen1 <<"\t"<< iegen2 <<"\t"<< iegen3 <<"\t"<< iegen4 <<"\t"<< iegen5 <<"\t"
    << 0 <<"\t"<< igcen1 <<"\t"<< igcen2 <<"\t"<< igcen3 <<"\t"<< igcen4 <<"\t"<< igcen5 <<"\t"
    << 0 <<"\t"<< ir1    <<"\t"<< ir2    <<"\t"<< ir3    <<"\t"<< ir4    <<"\t"<< ir5 <<"\t"
    << 0 <<"\t"<< ir1    <<"\t"<< ir2    <<"\t"<< ir3    <<"\t"<< ir4    <<"\t"<< ir5 <<"\t"
    << 0 <<"\t"<< ir1    <<"\t"<< ir2    <<"\t"<< ir3    <<"\t"<< ir4    <<"\t"<< ir5 <<"\t"
    << endl;
file2 <<  // - Shower summary
    "#PrimCR"<<"\t"<< "LPCREn"  <<"\t"<< "HghPCR"  <<"\t"<< "FlinSh" <<"\t"<<
    "IShow"  <<"\t"<< "NSecP"   <<"\t"<< "NClst"   <<"\t"<<
    "NGamS"  <<"\t"<< "NEleS"   <<"\t"<< "NMuS"    <<"\t"<< "NNeutS" <<"\t"<<
    "NProtS" <<"\t"<< "NOthrS"  <<"\t"<<
    "NElClS" <<"\t"<< "NMuClS"  <<"\t"<< "NMxClS"  <<"\t"<< "NEGClS" <<"\t"<<
    "NGmClS" <<"\t"<< "NOtClS"  <<"\t"<<
    "ElMEnS" <<"\t"<< "MuMEnS"  <<"\t"<<
    "GmMnR"  <<"\t"<< "ElMnR"   <<"\t"<< "MuMnR"   <<"\t"<< "NtMnR"  <<"\t"<<  "PtMnR"   <<"\t"
    "OtMnR"  <<"\t"
    "ElCMnR" <<"\t"<< "MuCMnR"  <<"\t"<<  "MxCMnR"  <<"\t"
    "EGCMnR" <<"\t"<< "GmCMnR"  <<"\t"<<  "OtCMnR"  <<"\t"
    "NGmEn1" <<"\t"<< "NGmEn2"  <<"\t"<<  "NGmEn3"  <<"\t"<< "NGmEn4" <<"\t"<<
    "NGmEn5" <<"\t"<< "NGmEn6"  <<"\t"<<
    "NElEn1" <<"\t"<< "NElEn2"  <<"\t"<<  "NElEn3"  <<"\t"<< "NElEn4" <<"\t"<<
    "NElEn5" <<"\t"<< "NElEn6"  <<"\t"<<
    "NMuEn1" <<"\t"<< "NMuEn2"  <<"\t"<<  "NMuEn3"  <<"\t"<< "NMuEn4" <<"\t"<<
    "NMuEn5" <<"\t"<< "NMuEn6"  <<"\t"<<
    "NEClM1" <<"\t"<< "NEClM2"  <<"\t"<<  "NEClM3"  <<"\t"<< "NEClM4" <<"\t"<<
    "NEClM5" <<"\t"<< "NEClM6"  <<"\t"<<
    "NEClE1" <<"\t"<< "NEClE2"  <<"\t"<<  "NEClE3"  <<"\t"<< "NEClE4" <<"\t"<<
    "NEClE5" <<"\t"<< "NEClE6"  <<"\t"<<
    "NECdT1" <<"\t"<< "NECdT2"  <<"\t"<<  "NECdT3"  <<"\t"<< "NECdT4" <<"\t"<<
    "NECdT5" <<"\t"<< "NECdT6"  <<"\t"<<
    "NECTh1" <<"\t"<< "NECTh2"  <<"\t"<<  "NECTh3"  <<"\t"<< "NECTh4" <<"\t"<<
    "NECTh5" <<"\t"<< "NECTh6"  <<"\t"<<
    "NECSz1" <<"\t"<< "NECSz2"  <<"\t"<<  "NECSz3"  <<"\t"<< "NECSz4" <<"\t"<<
    "NECSz5" <<"\t"<< "NECSz6"  <<"\t"<<
    "NEGEn1" <<"\t"<< "NEGEn2"  <<"\t"<<  "NEGEn3"  <<"\t"<< "NEGEn4" <<"\t"<<
    "NEGEn5" <<"\t"<< "NEGEn6"  <<"\t"<<
    "NGCEn1" <<"\t"<< "NGCEn2"  <<"\t"<<  "NGCEn3"  <<"\t"<< "NGCEn4" <<"\t"<<
    "NGCEn5" <<"\t"<< "NGCEn6"  <<"\t"<<
    "NEleR1" <<"\t"<< "NEleR2"  <<"\t"<<  "NEleR3"  <<"\t"<< "NEleR4" <<"\t"<<
    "NEleR5" <<"\t"<< "NEleR6"  <<"\t"<<
    "NMunR1" <<"\t"<< "NMunR2"  <<"\t"<<  "NMunR3"  <<"\t"<< "NMunR4"  <<"\t"<<
    "NMunR5" <<"\t"<< "NMunR6"  <<"\t"<<
    "NEClR1" <<"\t"<< "NEClR2"  <<"\t"<<  "NEClR3"  <<"\t"<< "NEClR4" <<"\t"<<
    "NEClR5" <<"\t"<< "NEClR6"  <<"\t"<<
    endl ;
//
file3 << "# Main Summary: DetMxD, DetMSg, PCRZhMin, PCRZhMx: " <<"\t"<<
        mxdst <<"\t"<< sigrmx <<"\t"<< zhmin <<"\t"<< zhmax <<"\t"<<
        " HMin, HMax: "<<"\t"<< hmin <<"\t"<< hmax << endl;
 //       "Height intervals:  0, 12, 16, 18, 20, 22, 28, 40, 150" << endl;
file3 << "# \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t  \t \t \t \t \t \t \t \t"
    << 0 <<"\t"<< igmen1 <<"\t"<< igmen2 <<"\t"<< igmen3 <<"\t"<< igmen4 <<"\t"<< igmen5 <<"\t"
    << 0 <<"\t"<< ielen1 <<"\t"<< ielen2 <<"\t"<< ielen3 <<"\t"<< ielen4 <<"\t"<< ielen5 <<"\t"
    << 0 <<"\t"<< imuen1 <<"\t"<< imuen2 <<"\t"<< imuen3 <<"\t"<< imuen4 <<"\t"<< imuen5 <<"\t"
    << 0 <<"\t"<< inten1 <<"\t"<< inten2 <<"\t"<< inten3 <<"\t"<< inten4 <<"\t"<< inten5 <<"\t"
    << 0 <<"\t"<< ipren1 <<"\t"<< ipren2 <<"\t"<< ipren3 <<"\t"<< ipren4 <<"\t"<< ipren5 <<"\t"
    << ieclm1 <<"\t"<< ieclm2 <<"\t"<< ieclm3 <<"\t"<< ieclm4 <<"\t"<< ieclm5 <<"\t"<< ieclm5+1 <<"\t"
    << 0 <<"\t"<< iecen1 <<"\t"<< iecen2 <<"\t"<< iecen3 <<"\t"<< iecen4 <<"\t"<< iecen5 <<"\t"
    << 0 <<"\t"<< iecdt1 <<"\t"<< iecdt2 <<"\t"<< iecdt3 <<"\t"<< iecdt4 <<"\t"<< iecdt5 <<"\t"
    << 0 <<"\t"<< iecth1 <<"\t"<< iecth2 <<"\t"<< iecth3 <<"\t"<< iecth4 <<"\t"<< iecth5 <<"\t"
    << 0 <<"\t"<< iecsz1 <<"\t"<< iecsz2 <<"\t"<< iecsz3 <<"\t"<< iecsz4 <<"\t"<< iecsz5 <<"\t"
    << 0 <<"\t"<< ieced1 <<"\t"<< ieced2 <<"\t"<< ieced3 <<"\t"<< ieced4 <<"\t"<< ieced5 <<"\t"
    << 0 <<"\t"<< iecem1 <<"\t"<< iecem2 <<"\t"<< iecem3 <<"\t"<< iecem4 <<"\t"<< iecem5 <<"\t"
    << 0 <<"\t"<< iecet1 <<"\t"<< iecet2 <<"\t"<< iecet3 <<"\t"<< iecet4 <<"\t"<< iecet5 <<"\t"
    << 0 <<"\t"<< ieest1 <<"\t"<< ieest2 <<"\t"<< ieest3 <<"\t"<< ieest4 <<"\t"<< ieest5 <<"\t"
    << 0 <<"\t"<< iegen1 <<"\t"<< iegen2 <<"\t"<< iegen3 <<"\t"<< iegen4 <<"\t"<< iegen5 <<"\t"
    << 0 <<"\t"<< igcen1 <<"\t"<< igcen2 <<"\t"<< igcen3 <<"\t"<< igcen4 <<"\t"<< igcen5 <<"\t"
    << endl;
file3 <<   // -Main summary
    "#PrimCR"<<"\t"<< "LPCREn" <<"\t"<< "MnHPCR" <<"\t"<< "FlxInt" <<"\t"<<
    "NShAna" <<"\t"<< "NtSePs" <<"\t"<< "NtClts" <<"\t"<<
    "NtGam"  <<"\t"<< "NtEle"  <<"\t"<< "NtMu"   <<"\t"<< "NtNeu"  <<"\t"<<
    "NtPrt"  <<"\t"<< "NtOth"  <<"\t"<<
    "NElCl"  <<"\t"<< "NMuCl"  <<"\t"<< "NMxCl"  <<"\t"<< "NEGCl"  <<"\t"<<
    "NGmCl"  <<"\t"<< "NOtCl"  <<"\t"<<
    "GmMnR"  <<"\t"<< "ElMnR"  <<"\t"<< "MuMnR"  <<"\t"<< "NtMnR"  <<"\t"<<  "PtMnR"   <<"\t"
    "OtMnR"  <<"\t"<<
    "ElCMnR" <<"\t"<< "MuCMnR" <<"\t"<< "MxCMnR" <<"\t"
    "EGCMnR" <<"\t"<< "GmCMnR" <<"\t"<< "OtCMnR" <<"\t"
    "NGmEn1" <<"\t"<< "NGmEn2" <<"\t"<< "NGmEn3" <<"\t"<< "NGmEn4" <<"\t"<<
    "NGmEn5" <<"\t"<< "NGmEn6" <<"\t"<<
    "NElEn1" <<"\t"<< "NElEn2" <<"\t"<< "NElEn3" <<"\t"<< "NElEn4" <<"\t"<<
    "NElEn5" <<"\t"<< "NElEn6" <<"\t"<<
    "NMuEn1" <<"\t"<< "NMuEn2" <<"\t"<< "NMuEn3" <<"\t"<< "NMuEn4" <<"\t"<<
    "NMuEn5" <<"\t"<< "NMuEn6" <<"\t"<<
    "NNtEn1" <<"\t"<< "NNtEn2" <<"\t"<< "NNtEn3" <<"\t"<< "NNtEn4" <<"\t"<<
    "NNtEn5" <<"\t"<< "NNtEn6" <<"\t"<<
    "NPrEn1" <<"\t"<< "NPrEn2" <<"\t"<< "NPrEn3" <<"\t"<< "NPrEn4" <<"\t"<<
    "NPrEn5" <<"\t"<< "NPrEn6" <<"\t"<<
    "NEClM1" <<"\t"<< "NEClM2" <<"\t"<< "NEClM3" <<"\t"<< "NEClM4" <<"\t"<<
    "NEClM5" <<"\t"<< "NEClM6" <<"\t"<<
    "NEClE1" <<"\t"<< "NEClE2" <<"\t"<< "NEClE3" <<"\t"<< "NEClE4" <<"\t"<<
    "NEClE5" <<"\t"<< "NEClE6" <<"\t"<<
    "NECdT1" <<"\t"<< "NECdT2" <<"\t"<< "NECdT3" <<"\t"<< "NECdT4" <<"\t"<<
    "NECdT5" <<"\t"<< "NECdT6" <<"\t"<<
    "NECTh1" <<"\t"<< "NECTh2" <<"\t"<< "NECTh3" <<"\t"<< "NECTh4" <<"\t"<<
    "NECTh5" <<"\t"<< "NECTh6" <<"\t"<<
    "NECSz1" <<"\t"<< "NECSz2" <<"\t"<< "NECSz3" <<"\t"<< "NECSz4" <<"\t"<<
    "NECSz5" <<"\t"<< "NECSz6" <<"\t"<<
    "NECED1" <<"\t"<< "NECED2" <<"\t"<< "NECED3" <<"\t"<< "NECED4" <<"\t"<<
    "NECED5" <<"\t"<< "NECED6" <<"\t"<<
    "NECEM1" <<"\t"<< "NECEM2" <<"\t"<< "NECEM3" <<"\t"<< "NECEM4" <<"\t"<<
    "NECEM5" <<"\t"<< "NECEM6" <<"\t"<<
    "NECET1" <<"\t"<< "NECET2" <<"\t"<< "NECET3" <<"\t"<< "NECET4" <<"\t"<<
    "NECET5" <<"\t"<< "NECET6" <<"\t"<<
    "NEETS1" <<"\t"<< "NEETS2" <<"\t"<< "NEETS3" <<"\t"<< "NEETS4" <<"\t"<<
    "NEETS5" <<"\t"<< "NEETS6" <<"\t"<<
    "NEGEn1" <<"\t"<< "NEGEn2" <<"\t"<< "NEGEn3" <<"\t"<< "NEGEn4" <<"\t"<<
    "NEGEn5" <<"\t"<< "NEGEn6" <<"\t"<<
    "NGCEn1" <<"\t"<< "NGCEn2" <<"\t"<< "NGCEn3" <<"\t"<< "NGCEn4" <<"\t"<<
    "NGCEn5" <<"\t"<< "NGCEn6" <<"\t"<<
    endl ;
file4 << "# Radial Summary: DetMxD, DetMSg, PCRZhMin, PCRZhMx: " <<"\t"<<
        mxdst <<"\t"<< sigrmx <<"\t"<< zhmin <<"\t"<< zhmax <<"\t"<<
        " HMin, HMax: "<<"\t"<< hmin <<"\t"<< hmax << endl;
//        "Height intervals:  0, 12, 16, 18, 20, 22, 28, 40, 150" << endl;
file4 << "# \t \t \t \t \t \t \t"
    << 0 <<"\t"<< ir1    <<"\t"<< ir2    <<"\t"<< ir3    <<"\t"<< ir4    <<"\t"<< ir5 <<"\t"
    << 0 <<"\t"<< ir1    <<"\t"<< ir2    <<"\t"<< ir3    <<"\t"<< ir4    <<"\t"<< ir5 <<"\t"
    << 0 <<"\t"<< ir1    <<"\t"<< ir2    <<"\t"<< ir3    <<"\t"<< ir4    <<"\t"<< ir5 <<"\t"
    << 0 <<"\t"<< ir1    <<"\t"<< ir2    <<"\t"<< ir3    <<"\t"<< ir4    <<"\t"<< ir5 <<"\t"
    << 0 <<"\t"<< ir1    <<"\t"<< ir2    <<"\t"<< ir3    <<"\t"<< ir4    <<"\t"<< ir5 <<"\t"
    << 0 <<"\t"<< ir1    <<"\t"<< ir2    <<"\t"<< ir3    <<"\t"<< ir4    <<"\t"<< ir5 <<"\t"
    << 0 <<"\t"<< ir1    <<"\t"<< ir2    <<"\t"<< ir3    <<"\t"<< ir4    <<"\t"<< ir5 <<"\t"
    << 0 <<"\t"<< ir1    <<"\t"<< ir2    <<"\t"<< ir3    <<"\t"<< ir4    <<"\t"<< ir5 <<"\t"
    << 0 <<"\t"<< ir1    <<"\t"<< ir2    <<"\t"<< ir3    <<"\t"<< ir4    <<"\t"<< ir5 <<"\t"
    << 0 <<"\t"<< ir1    <<"\t"<< ir2    <<"\t"<< ir3    <<"\t"<< ir4    <<"\t"<< ir5 <<"\t"
    << 0 <<"\t"<< ir1    <<"\t"<< ir2    <<"\t"<< ir3    <<"\t"<< ir4    <<"\t"<< ir5 <<"\t"
    << 0 <<"\t"<< ir1    <<"\t"<< ir2    <<"\t"<< ir3    <<"\t"<< ir4    <<"\t"<< ir5 <<"\t"
    << 0 <<"\t"<< ir1    <<"\t"<< ir2    <<"\t"<< ir3    <<"\t"<< ir4    <<"\t"<< ir5 <<"\t"
    << endl;
file4 <<   // - Radial distributions summary for all clusters
    "#PrimCR"<<"\t"<< "LEnPCR" <<"\t"<<
    "MnHPCR" <<"\t"<<  //  PCR Mean height
    "FlxInt" <<"\t"<<  //  PCR Flux integral
    "NShAna" <<"\t"<< "NtSePs" <<"\t"<< "NtClts" <<"\t"<<
    "NElR1"  <<"\t"<< "NElR2"  <<"\t"<< "NElR3"  <<"\t"<< "NElR4"  <<"\t"<<
    "NElR5"  <<"\t"<< "NElR6"  <<"\t"<<
    "NMuR1"  <<"\t"<< "NMuR2"  <<"\t"<< "NMuR3"  <<"\t"<< "NMuR4"  <<"\t"<<
    "NMuR5"  <<"\t"<< "NMuR6"  <<"\t"<<
    "NEClR1" <<"\t"<< "NEClR2" <<"\t"<< "NEClR3" <<"\t"<< "NEClR4" <<"\t"<<
    "NEClR5" <<"\t"<< "NEClR6" <<"\t"<<
    "NMClR1" <<"\t"<< "NMClR2" <<"\t"<< "NMClR3" <<"\t"<< "NMClR4" <<"\t"<<
    "NMClR5" <<"\t"<< "NMClR6" <<"\t"<<
    "NXClR1" <<"\t"<< "NXClR2" <<"\t"<< "NXClR3" <<"\t"<< "NXClR4" <<"\t"<<
    "NXClR5" <<"\t"<< "NXClR6" <<"\t"<<
    "ElMER1" <<"\t"<< "ElMER2" <<"\t"<< "ElMER3" <<"\t"<< "ElMER4" <<"\t"<<
    "ElMER5" <<"\t"<< "ElMER6" <<"\t"<<
    "ElMTR1" <<"\t"<< "ElMTR2" <<"\t"<< "ElMTR3" <<"\t"<< "ElMTR4" <<"\t"<<
    "ElMTR5" <<"\t"<< "ElMTR6" <<"\t"<<
    "MuMER1" <<"\t"<< "MuMER2" <<"\t"<< "MuMER3" <<"\t"<< "MuMER4" <<"\t"<<
    "MuMER5" <<"\t"<< "MuMER6" <<"\t"<<
    "MuMTR1" <<"\t"<< "MuMTR2" <<"\t"<< "MuMTR3" <<"\t"<< "MuMTR4" <<"\t"<<
    "MuMTR5" <<"\t"<< "MuMTR6" <<"\t"<<
    "ECMER1" <<"\t"<< "ECMER2" <<"\t"<< "ECMER3" <<"\t"<< "ECMER4" <<"\t"<<
    "ECMER5" <<"\t"<< "ECMER6" <<"\t"<<
    "ECMTR1" <<"\t"<< "ECMTR2" <<"\t"<< "ECMTR3" <<"\t"<< "ECMTR4" <<"\t"<<
    "ECMTR5" <<"\t"<< "ECMTR6" <<"\t"<<
    "EMdTR1" <<"\t"<< "EMdTR2" <<"\t"<< "EMdTR3" <<"\t"<< "EMdTR4" <<"\t"<<
    "EMdTR5" <<"\t"<< "EMdTR6" <<"\t"<<
    "EMERR1" <<"\t"<< "EMERR2" <<"\t"<< "EMERR3" <<"\t"<< "EMERR4" <<"\t"<<
    "EMERR5" <<"\t"<< "EMERR6" <<"\t"<<
    endl;
    
file5 << "# Mean Heights Summary: DetMxD, DetMSg, PCRZhMin, PCRZhMx: " <<"\t"<<
            mxdst <<"\t"<< sigrmx <<"\t"<< zhmin <<"\t"<< zhmax <<"\t"<<
            " HMin, HMax: "<<"\t"<< hmin <<"\t"<< hmax << endl;
    //        "Height intervals:  0, 12, 16, 18, 20, 22, 28, 40, 150" << endl;
file5 << "# \t \t \t \t \t \t \t "
    << 0 <<"\t"<< igmen1 <<"\t"<< igmen2 <<"\t"<< igmen3 <<"\t"<< igmen4 <<"\t"<< igmen5 <<"\t"
    << 0 <<"\t"<< ielen1 <<"\t"<< ielen2 <<"\t"<< ielen3 <<"\t"<< ielen4 <<"\t"<< ielen5 <<"\t"
    << 0 <<"\t"<< imuen1 <<"\t"<< imuen2 <<"\t"<< imuen3 <<"\t"<< imuen4 <<"\t"<< imuen5 <<"\t"
    << 0 <<"\t"<< inten1 <<"\t"<< inten2 <<"\t"<< inten3 <<"\t"<< inten4 <<"\t"<< inten5 <<"\t"
    << 0 <<"\t"<< ipren1 <<"\t"<< ipren2 <<"\t"<< ipren3 <<"\t"<< ipren4 <<"\t"<< ipren5 <<"\t"
    << ieclm1 <<"\t"<< ieclm2 <<"\t"<< ieclm3 <<"\t"<< ieclm4 <<"\t"<< ieclm5 <<"\t"<< ieclm5+1 <<"\t"
    << 0 <<"\t"<< iecen1 <<"\t"<< iecen2 <<"\t"<< iecen3 <<"\t"<< iecen4 <<"\t"<< iecen5 <<"\t"
    << 0 <<"\t"<< iecdt1 <<"\t"<< iecdt2 <<"\t"<< iecdt3 <<"\t"<< iecdt4 <<"\t"<< iecdt5 <<"\t"
    << 0 <<"\t"<< iecth1 <<"\t"<< iecth2 <<"\t"<< iecth3 <<"\t"<< iecth4 <<"\t"<< iecth5 <<"\t"
    << 0 <<"\t"<< iecsz1 <<"\t"<< iecsz2 <<"\t"<< iecsz3 <<"\t"<< iecsz4 <<"\t"<< iecsz5 <<"\t"
    << 0 <<"\t"<< ieced1 <<"\t"<< ieced2 <<"\t"<< ieced3 <<"\t"<< ieced4 <<"\t"<< ieced5 <<"\t"
    << 0 <<"\t"<< iecem1 <<"\t"<< iecem2 <<"\t"<< iecem3 <<"\t"<< iecem4 <<"\t"<< iecem5 <<"\t"
    << 0 <<"\t"<< iecet1 <<"\t"<< iecet2 <<"\t"<< iecet3 <<"\t"<< iecet4 <<"\t"<< iecet5 <<"\t"
    << 0 <<"\t"<< ieest1 <<"\t"<< ieest2 <<"\t"<< ieest3 <<"\t"<< ieest4 <<"\t"<< ieest5 <<"\t"
    << 0 <<"\t"<< iegen1 <<"\t"<< iegen2 <<"\t"<< iegen3 <<"\t"<< iegen4 <<"\t"<< iegen5 <<"\t"
    << 0 <<"\t"<< igcen1 <<"\t"<< igcen2 <<"\t"<< igcen3 <<"\t"<< igcen4 <<"\t"<< igcen5 <<"\t"
        << endl;
file5 <<   // - Mean heigth summary for all clusters
    "#PrimCR"<<"\t"<< "LEnPCR" <<"\t"<<
    "MnHPCR" <<"\t"<<  //  PCR Mean height
    "FlxInt" <<"\t"<<  //  PCR Flux integral
    "NShAna" <<"\t"<< "NtSePs" <<"\t"<< "NtClts" <<"\t"<<
    "NGmEn1" <<"\t"<< "NGmEn2" <<"\t"<< "NGmEn3" <<"\t"<< "NGmEn4" <<"\t"<<
    "NGmEn5" <<"\t"<< "NGmEn6" <<"\t"<<
    "NElEn1" <<"\t"<< "NElEn2" <<"\t"<< "NElEn3" <<"\t"<< "NElEn4" <<"\t"<<
    "NElEn5" <<"\t"<< "NElEn6" <<"\t"<<
    "NMuEn1" <<"\t"<< "NMuEn2" <<"\t"<< "NMuEn3" <<"\t"<< "NMuEn4"  <<"\t"<<
    "NMuEn5" <<"\t"<< "NMuEn6" <<"\t"<<
    "NNtEn1" <<"\t"<< "NNtEn2" <<"\t"<< "NNtEn3" <<"\t"<< "NNtEn4"  <<"\t"<<
    "NNtEn5" <<"\t"<< "NNtEn6" <<"\t"<<
    "NPrEn1" <<"\t"<< "NPrEn2" <<"\t"<< "NPrEn3" <<"\t"<< "NPrEn4"  <<"\t"<<
    "NPrEn5" <<"\t"<< "NPrEn6" <<"\t"<<
    "NEClM1" <<"\t"<< "NEClM2" <<"\t"<< "NEClM3" <<"\t"<< "NEClM4" <<"\t"<<
    "NEClM5" <<"\t"<< "NEClM6" <<"\t"<<
    "NEClE1" <<"\t"<< "NEClE2" <<"\t"<< "NEClE3" <<"\t"<< "NEClE4" <<"\t"<<
    "NEClE5" <<"\t"<< "NEClE6" <<"\t"<<
    "NECdT1" <<"\t"<< "NECdT2" <<"\t"<< "NECdT3" <<"\t"<< "NECdT4" <<"\t"<<
    "NECdT5" <<"\t"<< "NECdT6" <<"\t"<<
    "NECTh1" <<"\t"<< "NECTh2" <<"\t"<< "NECTh3" <<"\t"<< "NECTh4" <<"\t"<<
    "NECTh5" <<"\t"<< "NECTh6" <<"\t"<<
    "NECSz1" <<"\t"<< "NECSz2" <<"\t"<< "NECSz3" <<"\t"<< "NECSz4" <<"\t"<<
    "NECSz5" <<"\t"<< "NECSz6" <<"\t"<<
    "NECED1" <<"\t"<< "NECED2" <<"\t"<< "NECED3" <<"\t"<< "NECED4" <<"\t"<<
    "NECED5" <<"\t"<< "NECED6" <<"\t"<<
    "NECEM1" <<"\t"<< "NECEM2" <<"\t"<< "NECEM3" <<"\t"<< "NECEM4" <<"\t"<<
    "NECEM5" <<"\t"<< "NECEM6" <<"\t"<<
    "NECET1" <<"\t"<< "NECET2" <<"\t"<< "NECET3" <<"\t"<< "NECET4" <<"\t"<<
    "NECET5" <<"\t"<< "NECET6" <<"\t"<<
    "NEETS1" <<"\t"<< "NEETS2" <<"\t"<< "NEETS3" <<"\t"<< "NEETS4" <<"\t"<<
    "NEETS5" <<"\t"<< "NEETS6" <<"\t"<<
    "NEGEn1" <<"\t"<< "NEGEn2" <<"\t"<< "NEGEn3" <<"\t"<< "NEGEn4" <<"\t"<<
    "NEGEn5" <<"\t"<< "NEGEn6" <<"\t"<<
    "NGCEn1" <<"\t"<< "NGCEn2" <<"\t"<< "NGCEn3" <<"\t"<< "NGCEn4" <<"\t"<<
    "NGCEn5" <<"\t"<< "NGCEn6" <<"\t"<<
    "MHGmE1" <<"\t"<< "MHGmE2" <<"\t"<< "MHGmE3" <<"\t"<< "MHGmE4" <<"\t"<<
    "MHGmE5" <<"\t"<< "MHGmE6" <<"\t"<<
    "MHElE1" <<"\t"<< "MHElE2" <<"\t"<< "MHElE3" <<"\t"<< "MHElE4" <<"\t"<<
    "MHElE5" <<"\t"<< "MHElE6" <<"\t"<<
    "MHMuE1" <<"\t"<< "MHMuE2" <<"\t"<< "MHMuE3" <<"\t"<< "MHMuE4" <<"\t"<<
    "MHMuE5" <<"\t"<< "MHMuE6" <<"\t"<<
    "MHNtE1" <<"\t"<< "MHNtE2" <<"\t"<< "MHNtE3" <<"\t"<< "MHNtE4" <<"\t"<<
    "MHNtE5" <<"\t"<< "MHNtE6" <<"\t"<<
    "MHPrE1" <<"\t"<< "MHPrE2" <<"\t"<< "MHPrE3" <<"\t"<< "MHPrE4" <<"\t"<<
    "MHPrE5" <<"\t"<< "MHPrE6" <<"\t"<<
    "MHECM1" <<"\t"<< "MHECM2" <<"\t"<< "MHECM3" <<"\t"<< "MHECM4" <<"\t"<<
    "MHECM5" <<"\t"<< "MHECM6" <<"\t"<<
    "MHECE1" <<"\t"<< "MHECE2" <<"\t"<< "MHECE3" <<"\t"<< "MHECE4" <<"\t"<<
    "MHECE5" <<"\t"<< "MHECE6" <<"\t"<<
    "MRECE1" <<"\t"<< "MRECE2" <<"\t"<< "MRECE3" <<"\t"<< "MRECE4" <<"\t"<<
    "MRECE5" <<"\t"<< "MRECE6" <<"\t"<<
    "MTECE1" <<"\t"<< "MTECE2" <<"\t"<< "MTECE3" <<"\t"<< "MTECE4" <<"\t"<<
    "MTECE5" <<"\t"<< "MTECE6" <<"\t"<<
    "MHEdT1" <<"\t"<< "MHEdT2" <<"\t"<< "MHEdT3" <<"\t"<< "MHEdT4" <<"\t"<<
    "MHEdT5" <<"\t"<< "MHEdT6" <<"\t"<<
    "MHETh1" <<"\t"<< "MHETh2" <<"\t"<< "MHETh3" <<"\t"<< "MHETh4" <<"\t"<<
    "MHETh5" <<"\t"<< "MHETh6" <<"\t"<<
    "MHESz1" <<"\t"<< "MHESz2" <<"\t"<< "MHESz3" <<"\t"<< "MHESz4" <<"\t"<<
    "MHESz5" <<"\t"<< "MHESz6" <<"\t"<<
    "MHEED1" <<"\t"<< "MHEED2" <<"\t"<< "MHEED3" <<"\t"<< "MHEED4" <<"\t"<<
    "MHEED5" <<"\t"<< "MHEED6" <<"\t"<<
    "MHEEM1" <<"\t"<< "MHEEM2" <<"\t"<< "MHEEM3" <<"\t"<< "MHEEM4" <<"\t"<<
    "MHEEM5" <<"\t"<< "MHEEM6" <<"\t"<<
    "MHEET1" <<"\t"<< "MHEET2" <<"\t"<< "MHEET3" <<"\t"<< "MHEET4" <<"\t"<<
    "MHEET5" <<"\t"<< "MHEET6" <<"\t"<<
    "MHETS1" <<"\t"<< "MHETS2" <<"\t"<< "MHETS3" <<"\t"<< "MHETS4" <<"\t"<<
    "MHETS5" <<"\t"<< "MHETS6" <<"\t"<<
    "MHEGE1" <<"\t"<< "MHEGE2" <<"\t"<< "MHEGE3" <<"\t"<< "MHEGE4" <<"\t"<<
    "MHEGE5" <<"\t"<< "MHEGE6" <<"\t"<<
    "MHGCE1" <<"\t"<< "MHGCE2" <<"\t"<< "MHGCE3" <<"\t"<< "MHGCE4" <<"\t"<<
    "MHGCE5" <<"\t"<< "MHGCE6" <<"\t"<<
    endl;
}
// -------------------------------------------------------------------------------------
Int_t ndr = 6,     ndt = 6,     nda = 6;         // nmen=8, nece=8;
Int_t vnemcr[ndr], vnmucr// We initialize some variables[ndr], vnmxcr[ndr],
      vnemct[ndt], vnmuct[ndt], vnmxct[ndt],
      vnemca[nda], vnmuca[nda], vnmxca[nda];
Int_t mclem[ndr][ndt][nda], mmucl[ndr][ndt][nda],  mmxcl[ndr][ndt][nda];
// We initialize some variables
for (Int_t ir=0; ir<ndr; ir++)
    for (Int_t it=0; it<ndt; it++)
        for (Int_t ia=0; ia< nda; ia++)
            {mclem[ir][it][ia]=0;    mmucl[ir][it][ia]=0;    mmxcl[ir][it][ia]=0;}
for (Int_t ie=0; ie<ndr; ie++){vnemcr[ie]=0; vnmucr[ie]=0; vnmxcr[ie]=0;} // ClRradius
for (Int_t it=0; it<ndt; it++){vnemct[it]=0; vnmuct[it]=0; vnmxct[it]=0;} // ClTimeWidth
for (Int_t ia=0; ia<nda; ia++){vnemca[ia]=0; vnmuca[ia]=0; vnmxca[ia]=0;} // ClAngWidth
//
// We analize nshana showers
// *************************************************************************************
//
Int_t nssam, isend;
nssam = nshana/nsamp;       // Num Shower in sample
isami = (isami-1) * nssam;  // First Shower in sample
isend = isami + nssam;      // Last Shower in sample
    
cout << "- NtShows. First Last: "  <<"\t"<< nshana <<"\t"<< isami+1 <<"\t"<< isend << endl;
cout << endl;

//for (Long64_t ishow=0; ishow < nshana; ishow++) {
for (Long64_t ishow=isami; ishow < isend; ishow++) {
    Long64_t itree = LoadTree(ishow);
    
    //cout << endl << " ishow "<< "\t" << ishow << endl;
    enemod = ishow%eperiod;
    
    if (itree < 0) break;
    nb = fChain->GetEntry(itree);
    nbytes += nb;
   
    
    // if (Cut(ientry) < 0) continue;
    eb  = b_shower_Energy->GetEntry(itree);        ebytes += eb;
    fb  = b_shower_FirstHeight->GetEntry(itree);   fbytes += fb;
    hb  = b_shower_Theta->GetEntry(itree);         hbytes += hb;
    tb  = b_shower_Phi->GetEntry(itree);           tbytes += tb;
    
    // Primary cosmic energy ---------------------------------
    pcrene = shower_Energy; // GeV
    if(pcrene < emin){emin = pcrene;}
    if(pcrene > emax){emax = pcrene;}
    lpcren  = log10(shower_Energy); //   log of primary energy
    lpcrem = lpcrem + (lpcren-lpcrem)/(ishow+1); // mean log(energy) of primary CR

    // First interaction height selection --------------------
    finhgt = shower_FirstHeight/100000; // h in km
    if((finhgt < hmin || finhgt > hmax )) continue;
    itsana++;               if (itsana == 1) {mnhgt = finhgt;}
    mnhgt = mnhgt + (finhgt-mnhgt)/itsana; //mean height of primary CR
    
    Float_t arr0[100000] = { [0 ... 99999] = -10. };
    Float_t arr1[100000] = { [0 ... 99999] = -10. };
    Float_t arr2[100000] = { [0 ... 99999] = -10. };
    Float_t arr3[100000] = { [0 ... 99999] = -10. };
    Float_t arr4[100000] = { [0 ... 99999] = -10. };
    Float_t arr5[100000] = { [0 ... 99999] = -10. };
    Float_t arr6[100000] = { [0 ... 99999] = -10. };
    Float_t itag[100000] = { [0 ... 99999] =  -1. };

    // Look for the fastest particle in the shower
    t0   = 1000000;
    t0r1 = 1000000; t0r2 = 1000000; t0r3 = 1000000; t0r4 = 1000000; t0r5 = 1000000; t0r6 = 1000000;
    //
    for(Int_t ip=0; ip < particle__; ip++){
        time      = particle__Time[ip];
        pid       = particle__ParticleID[ip];
        if(time < t0){ t0 = time; pidf = pid;}  // fastest particle
        // Fastest particle in each ring
        rx        = particle__x[ip]/100;    ry = particle__y[ip]/100;   // meters
        rr        = sqrt(rx*rx + ry*ry);
        //
        if (rr < ir1)             {if(time < t0r1){t0r1 = time; pidfr1 = pid;}}
        else if (rr>ir1 && rr<ir2){if(time < t0r2){t0r2 = time; pidfr2 = pid;}}
        else if (rr>ir2 && rr<ir3){if(time < t0r3){t0r3 = time; pidfr3 = pid;}}
        else if (rr>ir3 && rr<ir4){if(time < t0r4){t0r4 = time; pidfr4 = pid;}}
        else if (rr>ir4 && rr<ir5){if(time < t0r5){t0r5 = time; pidfr5 = pid;}}
        else {if(time < t0r6){t0r6 = time; pidfr6 = pid;}}
    }
    // radial t0 time
    t0r1 = t0r1 - t0;
    t0r2 = t0r2 - t0;
    t0r3 = t0r3 - t0;
    t0r4 = t0r4 - t0;
    t0r5 = t0r5 - t0;
    t0r6 = t0r6 - t0;
    //
    // First loop in particles ------------------------------------------------
    //particle__ = 1000;  //!
    nsecs = particle__;  // nb. secondary part. in shower
    // We initialize some counters
    isep = -1;  // index nb. of secondary particles
    ncluss = 0; // no. of clusters in the shower
    
// *************************************************************************************
    for(Int_t ip=0; ip < nsecs; ip++){   // First loop, just for counting and storing variables
        pid       = particle__ParticleID[ip];
        time      = particle__Time[ip] - t0;
        rx        = particle__x[ip]/100;   // meters
        ry        = particle__y[ip]/100;   // meters
        rr        = sqrt(rx*rx + ry*ry);
        px        = particle__Px[ip];
        py        = particle__Py[ip];
        pz        = particle__Pz[ip];
        psq       = px*px + py*py + pz*pz;
        pmod      = sqrt(psq);
        zha       = acos(pz/pmod); // * rad2g;
        aza       = atan2(py,px);  // * rad2g;
        
        // cout << " ***  " << ishow  << "\t" << enemod  << "\t" << ipartp  << "\t" << endl;
        
        //      Select type of particle:
        if(pid==1){ // GAMMAS
            if (pmod < gpcut) continue;  // We count only gammas with p > gpcut
            pic = icgam;    // particle code for filling the cluster id: clid
            ngams++;  ngamt++;
            rmgms = rmgms + (rr-rmgms)/ngams;   rmgmt = rmgmt + (rr-rmgmt)/ngamt;
            ene = pmod;
            if (ene < igmen1){ngmen1++;ngmes1++; hgmen1=hgmen1+finhgt;}
            else if (ene>igmen1 && ene<igmen2){ngmen2++;ngmes2++; hgmen2=hgmen2+finhgt;}
            else if (ene>igmen2 && ene<igmen3){ngmen3++;ngmes3++; hgmen3=hgmen3+finhgt;}
            else if (ene>igmen3 && ene<igmen4){ngmen4++;ngmes4++; hgmen4=hgmen4+finhgt;}
            else if (ene>igmen4 && ene<igmen5){ngmen5++;ngmes5++; hgmen5=hgmen5+finhgt;}
            else    {ngmen6++;ngmes6++; hgmen6=hgmen6+finhgt;}
        }
        else if( (pid==2) || (pid==3) ){ // ELECTRONS
            pic = icele;
            nels++; nelt++;
            rmels = rmels + (rr-rmels)/nels;    rmelt = rmelt + (rr-rmelt)/nelt;
            ene = sqrt(pmod*pmod + mele*mele);
            elmens = elmens + ene;    elment = elment + ene;
            // Energy selection
            if (ene<ielen1){nelen1++;neles1++; helen1=helen1+finhgt;}
            else if (ene>ielen1 && ene<ielen2){nelen2++;neles2++; helen2=helen2+finhgt;}
            else if (ene>ielen2 && ene<ielen3){nelen3++;neles3++; helen3=helen3+finhgt;}
            else if (ene>ielen3 && ene<ielen4){nelen4++;neles4++; helen4=helen4+finhgt;}
            else if (ene>ielen4 && ene<ielen5){nelen5++;neles5++; helen5=helen5+finhgt;}
            else    {nelen6++;neles6++; helen6=helen6+finhgt;}
            //   Electron radial distribution  (most probable distance)
            if (rr<ir1){nelr1++;nelrs1++;
                elmer1=elmer1+ene; emers1=emers1+ene;
                elmtr1=elmtr1+time-t0r1; emtrs1=emtrs1+time-t0r1;}
            else if (rr>ir1 && rr<ir2){nelr2++;nelrs2++;
                elmer2=elmer2+ene; emers2=emers2+ene;
                elmtr2=elmtr2+time-t0r2; emtrs2=emtrs2+time-t0r2;}
            else if (rr>ir2 && rr<ir3){nelr3++;nelrs3++;
                elmer3=elmer3+ene; emers3=emers3+ene;
                elmtr3=elmtr3+time-t0r3; emtrs3=emtrs3+time-t0r3;}
            else if (rr>ir3 && rr<ir4){nelr4++;nelrs4++;
                elmer4=elmer4+ene; emers4=emers4+ene;
                elmtr4=elmtr4+time-t0r4; emtrs4=emtrs4+time-t0r4;}
            else if (rr>ir4 && rr<ir5){nelr5++;nelrs5++;
                elmer5=elmer5+ene; emers5=emers5+ene;
                elmtr5=elmtr5+time-t0r5; emtrs5=emtrs5+time-t0r5;}
            else {nelr6++;nelrs6++;
                elmer6=elmer6+ene; emers6=emers6+ene;
                elmtr6=elmtr6+time-t0r6; emtrs6=emtrs6+time-t0r6;}
        }
        else if( (pid==5) || (pid==6) ){ // MUONS
            pic = icmu;
            nmus++; nmut++;
            //   Muon radial distribution
            rmmus = rmmus + (rr-rmmus)/nmus;    rmmut = rmmut + (rr-rmmut)/nmut;
            ene = sqrt(mmu*mmu + pmod*pmod);
            mument = mument + ene;
            mumens = mumens + ene;
            // Radial selection
            if (rr<ir1){nmur1++;nmurs1++;
                mumer1=mumer1+ene; mmers1= mmers1+ene;
                mumtr1=mumtr1+time-t0r1; mmtrs1= mmtrs1+time-t0r1;}
            else if (rr>ir1 && rr<ir2){nmur2++;nmurs2++;
                mumer2=mumer2+ene; mmers2= mmers2+ene;
                mumtr2=mumtr2+time-t0r2; mmtrs2= mmtrs2+time-t0r2;}
            else if (rr>ir2 && rr<ir3){nmur3++;nmurs3++;
                mumer3=mumer3+ene; mmers3= mmers3+ene;
                mumtr3=mumtr3+time-t0r3; mmtrs3= mmtrs3+time-t0r3;}
            else if (rr>ir3 && rr<ir4){nmur4++;nmurs4++;
                mumer4=mumer4+ene; mmers4= mmers4+ene;
                mumtr4=mumtr4+time-t0r4; mmtrs4= mmtrs4+time-t0r4;}
            else if (rr>ir4 && rr<ir5){nmur5++;nmurs5++;
                mumer5=mumer5+ene; mmers5= mmers5+ene;
                mumtr5=mumtr5+time-t0r5; mmtrs5= mmtrs5+time-t0r5;}
            else {nmur6++;nmurs6++;
                mumer6=mumer6+ene; mmers6= mmers6+ene;
                mumtr6=mumtr6+time-t0r6; mmtrs6= mmtrs6+time-t0r6;}
            // Energy selection
            if (ene < imuen1){nmuen1++;nmues1++; hmuen1=hmuen1+finhgt;}
            else if (ene>imuen1 && ene<imuen2){nmuen2++;nmues2++; hmuen2=hmuen2+finhgt;}
            else if (ene>imuen2 && ene<imuen3){nmuen3++;nmues3++; hmuen3=hmuen3+finhgt;}
            else if (ene>imuen3 && ene<imuen4){nmuen4++;nmues4++; hmuen4=hmuen4+finhgt;}
            else if (ene>imuen4 && ene<imuen5){nmuen5++;nmues5++; hmuen5=hmuen5+finhgt;}
            else    {nmuen6++;nmues6++; hmuen6=hmuen6+finhgt;}
        }
        else if( pid==13 ){ // NEUTRONS
            pic = icnt;   // nucleon
            nnt++;  nns++;
            rmnts = rmnts + (rr-rmnts)/nns;     rmntt = rmntt + (rr-rmntt)/nnt;
            ene  = sqrt(mn*mn + pmod*pmod);
            // Energy selection
            if (ene < inten1){nnten1++;nntes1++; hnten1=hnten1+finhgt;}
            else if (ene>inten1 && ene<inten2){nnten2++;nntes2++; hnten2=hnten2+finhgt;}
            else if (ene>inten2 && ene<inten3){nnten3++;nntes3++; hnten3=hnten3+finhgt;}
            else if (ene>inten3 && ene<inten4){nnten4++;nntes4++; hnten4=hnten4+finhgt;}
            else if (ene>inten4 && ene<inten5){nnten5++;nntes5++; hnten5=hnten5+finhgt;}
            else    {nnten6++;nntes6++; hnten6=hnten6+finhgt;}
        }
        else if( pid==14){ // PROTONS
            pic = icpt;   // nucleon
            npt++;  nps++;
            rmpts = rmpts + (rr-rmpts)/nps;     rmptt = rmptt + (rr-rmptt)/npt;
            ene  = sqrt(mp*mp + pmod*pmod);
            // Energy selection
            if (ene < ipren1){npren1++;npres1++; hpren1=hpren1+finhgt;}
            else if (ene>ipren1 && ene<ipren2){npren2++;npres2++; hpren2=hpren2+finhgt;}
            else if (ene>ipren2 && ene<ipren3){npren3++;npres3++; hpren3=hpren3+finhgt;}
            else if (ene>ipren3 && ene<ipren4){npren4++;npres4++; hpren4=hpren4+finhgt;}
            else if (ene>ipren4 && ene<ipren5){npren5++;npres5++; hpren5=hpren5+finhgt;}
            else    {npren6++;npres6++; hpren6=hpren6+finhgt;}
        }
        else{             // other particles
            nopt++; nops++;
            rmots = rmots + (rr-rmots)/nops;     rmott = rmott + (rr-rmott)/nopt;
            ene  = sqrt(mother*mother + pmod*pmod);
        }   // endif in pid
        //
        isep++;     // index of secondary particles
        arr0[isep] = pic;    // picode
        arr1[isep] = rx;     // x - coordinate
        arr2[isep] = ry;     // y - coordinate
        arr3[isep] = time;   // corrected arrival time
        arr4[isep] = pmod;   // momentum
        arr5[isep] = zha;    // zenith angle
        arr6[isep] = aza;    // azimuth angle
        //
        if ( ip == 1 & enemod == 1){ // prints up to mxparp particles
        // Particles summary
        // cout << "ishow, enemod, pid:   "<< "\t" << ishow << "\t" << enemod << "\t" << pid << endl;
        file0 << mpcr <<"\t"<< setprecision(4)<< lpcren   <<"\t"<<
            finhgt <<"\t"<<
            setprecision(0) << ishow+1 <<"\t"<< pid  <<"\t"<<
            setprecision(4) << rr      <<"\t"<< ene <<"\t"<<
            aza    <<"\t"<< zha <<"\t"<< endl;
    }    // end of first scan of particles
        
    ntsep  = ntsep + nsecs;   // total nb. secondary part.
    // Radial distribution of mean arrival time and energies
    emers1 = emers1/(nelrs1+0.001); emers2 = emers2/(nelrs2+0.001);
    emers3 = emers3/(nelrs3+0.001); emers4 = emers4/(nelrs4+0.001);
    emers5 = emers5/(nelrs5+0.001); emers6 = emers6/(nelrs6+0.001);
    emtrs1 = emtrs1/(nelrs1+0.001); emtrs2 = emtrs2/(nelrs2+0.001);
    emtrs3 = emtrs3/(nelrs3+0.001); emtrs4 = emtrs4/(nelrs4+0.001);
    emtrs5 = emtrs5/(nelrs5+0.001); emtrs6 = emtrs6/(nelrs6+0.001);
    mmers1 = mmers1/(nmurs1+0.001); mmers2 = mmers2/(nmurs2+0.001);
    mmers3 = mmers3/(nmurs3+0.001); mmers4 = mmers4/(nmurs4+0.001);
    mmers5 = mmers5/(nmurs5+0.001); mmers6 = mmers6/(nmurs6+0.001);
    mmtrs1 = mmtrs1/(nmurs1+0.001); mmtrs2 = mmtrs2/(nmurs2+0.001);
    mmtrs3 = mmtrs3/(nmurs3+0.001); mmtrs4 = mmtrs4/(nmurs4+0.001);
    mmtrs5 = mmtrs5/(nmurs5+0.001); mmtrs6 = mmtrs6/(nmurs6+0.001);
    //
    // Time difference between electrons and muouns
    emdts1 = emtrs1 - mmtrs1;      emdts2 = emtrs2 - mmtrs2;
    emdts3 = emtrs3 - mmtrs3;      emdts4 = emtrs4 - mmtrs4;
    emdts5 = emtrs5 - mmtrs5;      emdts6 = emtrs6 - mmtrs6;
    // Energy ratio between electrons and muouns
    if(mmers1!= 0){emrts1 = emers1/mmers1;}
    if(mmers2!= 0){emrts2 = emers2/mmers2;}
    if(mmers3!= 0){emrts3 = emers3/mmers3;}
    if(mmers4!= 0){emrts4 = emers4/mmers4;}
    if(mmers5!= 0){emrts5 = emers5/mmers5;}
    if(mmers6!= 0){emrts6 = emers6/mmers6;}
    //
    elmens = elmens/(nels + 0.001);
    mumens = mumens/(nmus + 0.001);  // rate mean elec. energy, mean muon. energy
    //
    emrats = elmens/mumens;          // electron muon energy ratio
    //
    }
    //
    mrgam=0;  mrele=0;  mrmuon=0; mrneut=0; mrprot=0; mroth=0;
// ===================================================================================
// ==================================================  ifp loop:  Looking for clusters
    for(Int_t ifp=0; ifp<nsecs; ifp++){
        if(itag[ifp]!=-1) continue;
        // First particle properties
        pic1  = arr0[ifp];
        if(pic1<0) continue;
        x1    = arr1[ifp];
        y1    = arr2[ifp];
        t1    = arr3[ifp];
        pm1   = arr4[ifp];
        zen1  = arr5[ifp];
        azh1  = arr6[ifp];
        //
        xfpc  = x1;   xlpc  = x1;    yfpc = y1;    ylpc = y1;    tfpc= t1; tlpc   = t1;
        zhfp  = zen1; zhlp  = zen1;  azfp = azh1;  azlp = azh1;
        pfpc  = pm1;  plpc  = pm1;
        picfp = pic1; piclp = pic1;
        //
        if(pic1 == icgam){
            if(pm1 < gpcut) continue; // We ignore gammas with p < gpcut
            ngamc = 1;
            tfg = t1;   tlg = t1;
            gmten = pm1;
            pxgmc = pm1 * sin(zen1) * cos(azh1);
            pygmc = pm1 * sin(zen1) * sin(azh1);
            pzgmc = pm1 * cos(zen1);
        }
        if(pic1 == icele){      // electron
            nelec = 1;
            tfe = t1;   tle = t1;
            // We calculate new variables of the electron cluster
            elten = sqrt(mele*mele + pm1*pm1);     // energy of the electron cluster
            pxelc = pm1 * sin(zen1) * cos(azh1);
            pyelc = pm1 * sin(zen1) * sin(azh1);
            pzelc = pm1 * cos(zen1);
        }
        if(pic1 == icmu){      // muon
            nmuc  = 1;
            tfm = t1;   tlm = t1;
        }
        //
        itag[ifp]++;    cmult = 1;      clid = pic1;
        //
        // =====================================
        // =====================================  isp loop: Adding particles to clusters
        for(Int_t isp = ifp+1; isp < nsecs; isp++){
        // Looking for new particles in the cluster
            if(itag[isp]!=-1) continue;
            pic2  = arr0[isp];
            if(pic2<0) continue;
            x2    = arr1[isp];  y2   = arr2[isp];   t2    = arr3[isp];
            pm2   = arr4[isp];  zen2 = arr5[isp];   azh2  = arr6[isp];
            //
            if(pic2 == icgam & pm2 < gpcut){    // we ignore gammas with p < gpcut
                //cout << endl << " --- isp: gamma ignored. pm= " << pm2 << endl << endl;
                tag++;  itag[isp] = tag; tag = -1;   continue;}
            //
            xcmno = x1;     ycmno = y1;     tcmno = t1;   // old values
            dxc   = x2 - xcmno;         dyc  = y2 - ycmno;  dtc   = t2 - tcmno;
            drcsq = dxc*dxc + dyc*dyc;  drc  = sqrt(drcsq);
            ssxo  = ssx;    ssyo = ssy;     ssto = sst;
            //
            if(drc <= mxdst){      // New candidate particle
                // Recursive formulae
                cmulp1  = cmult + 1;
                xclmn = xcmno + dxc/cmulp1;
                yclmn = ycmno + dyc/cmulp1;
                tclmn = tcmno + dtc/cmulp1;
                ssx  = ssxo + dxc*dxc * (cmulp1/cmult) ;     sigxsq = ssx/cmult;
                ssy  = ssyo + dyc*dyc * (cmulp1/cmult) ;     sigysq = ssy/cmult;
                sst  = ssto + dtc*dtc * (cmulp1/cmult) ;     sigtsq = sst/cmult;
                sigx = sqrt(sigxsq);  sigy = sqrt(sigysq);   sigtc = sqrt(sigtsq);
                ssxo = ssx;           ssyo = ssy;            ssto = sst;
                sigrc = sqrt((sigxsq + sigysq)/2) ;       // cluster stdev

                if(sigrc > sigrmx) continue;     // We add a new particle to the cluster
                //cout << "* nclust, cmult, drc, sigrc sigrmx "  << "\t" << nclust  << "\t" << cmult  << "\t" << drc  << "\t" << sigrc << "\t" << sigrmx << endl;
                cmult++;          clid  = clid + pic2;  tag++;
                xcmno = xclmn;    ycmno = yclmn;        tcmno = tclmn;
                //
                if(t2<=tfpc){picfp = pic2; xfpc = x2; yfpc = y2; tfpc = t2; pfpc = pm2;
                    zhfp = zen2; azfp = azh2;}
                if(t2>tlpc) {piclp = pic2; xlpc = x2; ylpc = y2; tlpc = t2; plpc = pm2;
                    zhlp = zen2; azlp = azh2;}
                // Normal components of first and last particles
                nxfp = sin(zhfp)*cos(azfp); nyfp = sin(zhfp)*sin(azfp); nzfp = cos(zhfp);
                nxlp = sin(zhlp)*cos(azlp); nylp = sin(zhlp)*sin(azlp); nzlp = cos(zhlp);
                clang = acos(nxfp*nxlp+nyfp*nylp+nzfp*nzlp)*rad2g; //cl. openning angle
                //clsz  = sqrt((xfpc-xlpc)*(xfpc-xlpc) + (yfpc-ylpc)*(yfpc-ylpc));
                cltsz = sigrc;
                cltdt = tlpc - tfpc + 0.01;                     // cluster time width
                rclmn = sqrt(xclmn*xclmn + yclmn*yclmn);   // distance to shower center
                //
                if(pic2 == icgam){  // gamma
                    if(ngamc==0){tfg = t2; tlg = t2;}
                    if(t2<=tfg){ tfg = t2;}
                    if(t2> tlg){ tlg = t2;}
                    ngamc++;
                    gmten = gmten + pm2;  // total gamma energy in the cluster
                    pxgmc = pxgmc + pm2 * sin(zen1) * cos(azh1);
                    pygmc = pygmc + pm2 * sin(zen1) * sin(azh1);
                    pzgmc = pzgmc + pm2 * cos(zen1);
                    // egcth  = atan2(sqrt(pxegc*pxegc+pyegc*pyegc), pzegc) * 90/(TMath::Pi()/2); // "Thermalicity"
                }
                if(pic2 == icele){  // electron
                    //cout << "- new electron ifp " << "\t" << ifp << endl;
                    if(nelec==0){tfe = t2; tle = t2;}
                    if(t2<=tfe){ tfe = t2;}
                    if(t2> tle){ tle = t2;}
                    nelec++;
                    ene = sqrt(mele*mele + pm2*pm2);
                    elten = elten + ene;
                    pxelc = pxelc + pm2 * sin(zen1) * cos(azh1);
                    pyelc = pyelc + pm2 * sin(zen1) * sin(azh1);
                    pzelc = pzelc + pm2 * cos(zen1);
                    //eclth = atan2(sqrt(pxelc*pxelc+pyelc*pyelc), pzelc) * 90/(TMath::Pi()/2); // "Thermalicity"
                }
                if(pic2 == icmu){   // muon
                    if(nmuc==0){tfm = t2; tlm = t2;}
                    if(t2<=tfm){ tfm = t2;}
                    if(t2> tlm){ tlm = t2;}
                    nmuc++;
                }
                //cout << "\t \t \t \t \t --- cmult" << "\t" << cmult << endl;
            }  //  endif distance
        
            itag[isp] = tag;   // tags particle
            tag       = -1;    //
        } //  End isp loop looking for more particles to the cluster
        //
        
        //cout << " ** ishow, clid, cmult, ngamc, nelec, nmuc:" << "\t" << ishow << "\t" << clid << "\t" << cmult << "\t" << ngamc << "\t" << nelec << "\t" << nmuc << endl;
        
        ssx = 0;    ssy = 0;    sst = 0;
        //
        if (cmult < 2 && nelec==0 ) { // ECluster + single electron filter
            cmult=0;    clid =0;    sclid = 0;
            ngamc=0;    nelec=0;    nmuc = 0;  gmten=0; elten=0; eclen=0; egcen=0; gmcen=0;
            continue;}
        // ==========================================     New cluster found
        //
        // cout << " \t \t New cluster found: ifp clid"  << "\t" << ifp << "\t" << clid << endl;
        //
        if (clid!=1000){     //  We ignore single electrons
            nclust++;   // Total no. of clusters
            ncluss++;   // No. of clusters in the shower
            // We classify a few general cluster variables
            if (cltdt < iecdt1){jt=0;}    // Time width of the cluster
            else if (cltdt> iecdt1 && cltdt<iecdt2){jt=1;}
            else if (cltdt> iecdt2 && cltdt<iecdt3){jt=2;}
            else if (cltdt> iecdt3 && cltdt<iecdt4){jt=3;}
            else if (cltdt> iecdt4 && cltdt<iecdt5){jt=4;}
            else {jt=5;}
            if (clang < iecan1){ja=0;}    // Openning angle of the cluster
            else if (clang> iecan1 && clang<iecan2){ja=1;}
            else if (clang> iecan2 && clang<iecan3){ja=2;}
            else if (clang> iecan3 && clang<iecan4){ja=3;}
            else if (clang> iecan4 && clang<iecan5){ja=4;}
            else {ja=5;}
            if (rclmn < ir1){jr=0;}  // Distance to Center of the Shower
            else if (rclmn>ir1 && rclmn<ir2){jr=1;}
            else if (rclmn>ir2 && rclmn<ir3){jr=2;}
            else if (rclmn>ir3 && rclmn<ir4){jr=3;}
            else if (rclmn>ir4 && rclmn<ir5){jr=4;}
            else {jr=5;}
        }
        //cout << " ***** ishow, clid, cmult, ngamc, nelec, nmuc:" << "\t" << "\t" << ishow << "\t" << clid << "\t" << cmult << "\t" << ngamc << "\t" << nelec << "\t" << nmuc << endl;
        if (cmult>1 && clid< icele){   // pure gamma clusters
            //cout << "\t --- Gammas: " << "\t" << clid << endl;
            sclid = 5;
            ngmcls++;
            ngmclt++;
            gmcsz  = cltsz;
            gmcmlt = cmult;
            gmcen  = gmten;
            gmcdt  = cltdt;
            rgclsq = rclmn*rclmn;   rgcl   = sqrt(rgclsq);
            rmgcls = rmgcls + (rgcl-rmgcls)/ngmcls;
            rmgclt = rmgclt + (rgcl-rmgclt)/ngmclt;
            gmcth  = atan2(sqrt(pxgmc*pxgmc+pygmc*pygmc), pzgmc) * 90/(TMath::Pi()/2);
            gmced  = gmcen / (gmcsz*gmcsz);     // energy density
            gmcem  = gmcen / ngamc;
            gmcet  = gmcen / gmcdt;
            gmcest = gmced / gmcdt;    // energy space-time density
            // "Thermalicity"
            // gamma cluster energy
            if (gmcen < igcen1){ngcen1++; ngces1++; hgcen1=hgcen1+finhgt;}  // GmCl. Energy
            else if (gmcen>igcen1 && gmcen<igcen2){ngcen2++; neces2++; hgcen2=hgcen2+finhgt;}
            else if (gmcen>igcen2 && gmcen<igcen3){ngcen3++; neces3++; hgcen3=hgcen3+finhgt;}
            else if (gmcen>igcen3 && gmcen<igcen4){ngcen4++; neces4++; hgcen4=hgcen4+finhgt;}
            else if (gmcen>igcen4 && gmcen<igcen5){ngcen5++; neces5++; hgcen5=hgcen5+finhgt;}
            else {ngcen6++; neces6++; hgcen6=hgcen6+finhgt;}
        }
        else if(nelec > 0 && clid < icmu ){  // EM Cluster (including gammas) + single electrons
            //
            //cout << "* nelec, clid: " << "\t" << nelec << "\t" << clid << endl;
            // ECluster multiplicity
            if (nelec == ieclm1){neclm1++; necms1++; heclm1=heclm1+finhgt;} // Cl Mult.
            else if(nelec == ieclm2){neclm2++; necms2++; heclm2=heclm2+finhgt;}
            else if(nelec == ieclm3){neclm3++; necms3++; heclm3=heclm3+finhgt;}
            else if(nelec == ieclm4){neclm4++; necms4++; heclm4=heclm4+finhgt;}
            else if(nelec == ieclm5){neclm5++; necms5++; heclm5=heclm5+finhgt;}
            else {neclm6++; necms6++; heclm6=heclm6+finhgt;}
            //cout << "** neclm1, m2, m3, m4, m5, m6 " << "\t" << neclm1 << "\t" << neclm2 << "\t" << neclm3 << "\t" << neclm4 << "\t" << neclm5 << "\t" << neclm6 << endl;
            eclen = elten;  //  only electron energy
            if (eclen < iecen1){necen1++; neces1++; hecen1=hecen1+finhgt; recen1=recen1+recl; tecen1=tecen1+ecldt;}  // Cl. Energy
            else if (eclen>iecen1 && eclen<iecen2){necen2++; neces2++; hecen2=hecen2+finhgt; recen2=recen2+recl; tecen2=tecen2+ecldt; }
            else if (eclen>iecen2 && eclen<iecen3){necen3++; neces3++; hecen3=hecen3+finhgt; recen3=recen3+recl; tecen3=tecen3+ecldt;}
            else if (eclen>iecen3 && eclen<iecen4){necen4++; neces4++; hecen4=hecen4+finhgt; recen4=recen4+recl; tecen4=tecen4+ecldt;}
            else if (eclen>iecen4 && eclen<iecen5){necen5++; neces5++; hecen5=hecen5+finhgt; recen5=recen5+recl; tecen5=tecen5+ecldt;}
            else {necen6++; neces6++; hecen6=hecen6+finhgt;  recen6=recen6+recl; tecen6=tecen6+ecldt;}
            //
            if (nelec>1) { // Multi electron clusters
                sclid  = 1;    // short IC of the cluster
                eclmlt = nelec;
                eclsz  = cltsz;
                nelcls ++;    // nb. of EM clusters in the shower
                nelclt ++;    // Total nb. of EM clusters
                //cout << "\t multielectron cluster: " << nelclt << endl;
                //cout << endl;
                ecldt   = cltdt;
                // Mean distance of the clusters in the shower
                reclsq = rclmn*rclmn;  recl    = sqrt(reclsq);
                rmecls = rmecls + (recl-rmecls)/nelcls;
                rmeclt = rmeclt + (recl-rmeclt)/nelclt;
                eclth  = atan2(sqrt(pxelc*pxelc+pyelc*pyelc), pzelc) * 90/(TMath::Pi()/2); // "Thermalicity"
                ecled   = eclen / (eclsz*eclsz);     // energy density
                eclem   = eclen / nelec;
                eclet   = eclen / ecldt;
                eclest  = ecled / ecldt;    // energy space-time density
                //cout << "--- eclest " << "\t" << eclest << endl;
                // ECluster time width
                if (ecldt < iecdt1){jt=0; necdt1++; nedts1++; hecdt1=hecdt1+finhgt;}
                else if (ecldt>= iecdt1 && ecldt<iecdt2){jt=1; necdt2++; nedts2++; hecdt2=hecdt2+finhgt;}
                else if (ecldt>= iecdt2 && ecldt<iecdt3){jt=2; necdt3++; nedts3++; hecdt3=hecdt3+finhgt;}
                else if (ecldt>= iecdt3 && ecldt<iecdt4){jt=3; necdt4++; nedts4++; hecdt4=hecdt4+finhgt;}
                else if (ecldt>= iecdt4 && ecldt<iecdt5){jt=4; necdt5++; nedts5++; hecdt5=hecdt5+finhgt;}
                else {jt=5; necdt6++; nedts6++; hecdt6=hecdt6+finhgt;}
                // ECluster Thermalicity
                if (eclth < iecth1){necth1++; neths1++; hecth1=hecth1+finhgt;}
                else if(eclth>=iecth1 && eclth<iecth2){necth2++; neths2++; hecth2=hecth2+finhgt;}
                else if(eclth>=iecth2 && eclth<iecth3){necth3++; neths3++; hecth3=hecth3+finhgt;}
                else if(eclth>=iecth3 && eclth<iecth4){necth4++; neths4++; hecth4=hecth4+finhgt;}
                else if(eclth>=iecth4 && eclth<iecth5){necth5++; neths5++; hecth5=hecth5+finhgt;}
                else {necth6++; neths6++; hecth6=hecth6+finhgt;}
                // Cluster size
                if (eclsz < iecsz1){necsz1++; neszs1++; hecsz1=hecsz1+finhgt;}
                else if (eclsz> iecsz1 && cltsz<iecsz2){necsz2++; neszs2++; hecsz2=hecsz2+finhgt;}
                else if (eclsz> iecsz2 && cltsz<iecsz3){necsz3++; neszs3++; hecsz3=hecsz3+finhgt;}
                else if (eclsz> iecsz3 && cltsz<iecsz4){necsz4++; neszs4++; hecsz4=hecsz4+finhgt;}
                else if (eclsz> iecsz4 && cltsz<iecsz5){necsz5++; neszs5++; hecsz5=hecsz5+finhgt;}
                else {necsz6++; neszs6++; hecsz6=hecsz6+finhgt;}
                // ECluster Energy density
                if (ecled < ieced1){neced1++; needs1++; heced1=heced1+finhgt;}
                else if (ecled>= ieced1 && ecled<ieced2){neced2++; needs2++; heced2=heced2+finhgt;}
                else if (ecled>= ieced2 && ecled<ieced3){neced3++; needs3++; heced3=heced3+finhgt;}
                else if (ecled>= ieced3 && ecled<ieced4){neced4++; needs4++; heced4=heced4+finhgt;}
                else if (ecled>= ieced4 && ecled<ieced5){neced5++; needs5++; heced5=heced5+finhgt;}
                else {neced6++; needs6++; heced6=heced6+finhgt;}
                // ECluster Energy/Multiplicy
                if (eclem < iecem1){necet1++; neems1++; hecem1=hecem1+finhgt;}
                else if (eclem>= iecem1 && eclem<iecem2){necem2++; neems2++; hecem2=hecem2+finhgt;}
                else if (eclem>= iecem2 && eclem<iecem3){necem3++; neems3++; hecem3=hecem3+finhgt;}
                else if (eclem>= iecem3 && eclem<iecem4){necem4++; neems4++; hecem4=hecem4+finhgt;}
                else if (eclem>= iecem4 && eclem<iecem5){necem5++; neems5++; hecem5=hecem5+finhgt;}
                else {necem6++; neems6++; hecem6=hecem6+finhgt;}
                // ECluster Energy/dTime
                if (eclet < iecet1){necet1++; neets1++; hecet1=hecet1+finhgt;}
                else if (eclet>= iecet1 && eclet<iecet2){necet2++; neets2++; hecet2=hecet2+finhgt;}
                else if (eclet>= iecet2 && eclet<iecet3){necet3++; neets3++; hecet3=hecet3+finhgt;}
                else if (eclet>= iecet3 && eclet<iecet4){necet4++; neets4++; hecet4=hecet4+finhgt;}
                else if (eclet>= iecet4 && eclet<iecet5){necet5++; neets5++; hecet5=hecet5+finhgt;}
                else {necet6++; neets6++; hecet6=hecet6+finhgt;}
                // ECluster EnergyDensity/dTime
                if (eclest < ieest1){neest1++; nests1++; heets1=heets1+finhgt;}
                else if (eclest>= ieest1 && eclet<ieest2){neest2++; nests2++; heets2=heets2+finhgt;}
                else if (eclest>= ieest2 && eclet<ieest3){neest3++; nests3++; heets3=heets3+finhgt;}
                else if (eclest>= ieest3 && eclet<ieest4){neest4++; nests4++; heets4=heets4+finhgt;}
                else if (eclest>= ieest4 && eclet<ieest5){neest5++; nests5++; heets5=heets5+finhgt;}
                else {neest6++; nests6++; heets6=heets6+finhgt;}
                //
                if (rclmn < ir1){neclr1++;necrs1++;   // Distance to Shower center
                    clmer1=clmer1+pmod; clmtr1=clmtr1+time-t0r1;}
                else if (rclmn>ir1 && rclmn<ir2){neclr2++;necrs2++;
                    clmer2=clmer2+pmod; clmtr2=clmtr2+time-t0r2;}
                else if (rclmn>ir2 && rclmn<ir3){neclr3++;necrs3++;
                    clmer3=clmer3+pmod; clmtr3=clmtr3+time-t0r3;}
                else if (rclmn>ir3 && rclmn<ir4){neclr4++;necrs4++;
                    clmer4=clmer4+pmod; clmtr4=clmtr4+time-t0r4;}
                else if (rclmn>ir4 && rclmn<ir5){neclr5++;necrs5++;
                    clmer5=clmer5+pmod; clmtr5=clmtr5+time-t0r5;}
                else {neclr6++;necrs6++; clmer6=clmer6+pmod; clmtr6=clmtr6+time-t0r6;}
                //
                mclem[jr][jt][ja]++;
                vnemcr[jr]++;  vnemct[jt]++;  vnemca[ja]++;}
                //
            else if (clid > 1000){ // Electron & Gamma Cluster
                sclid = 4;    // short IC of the cluster
                negcls ++;    // nb. of EM clusters in the shower
                negclt ++;    // Total nb. of EM clusters
                egcmlt = cmult;
                egcsz  = cltsz;
                egcen  = gmten + elten ;
                egcdt  = cltdt;
                pxegc  = pxelc + pxgmc;  pyegc = pyelc + pygmc;  pzegc = pzelc + pzgmc;
                egcth  = atan2(sqrt(pxegc*pxegc+pyegc*pyegc), pzegc) * 90/(TMath::Pi()/2);
                //
                regcsq = rclmn*rclmn;  regc    = sqrt(regcsq);
                rmegcs = rmegcs + (regc-rmegcs)/negcls;
                rmegct = rmegct + (regc-rmegct)/negclt;
                egced  = egcen / (egcsz*egcsz);     // energy density
                egcem  = egcen / (nelec+ngamc);
                egcet  = egcen / egcdt;
                egcest = egced / egcdt;    // energy space-time density
                // EGCluster energy
                if (egcen < iegen1){negen1++; neges1++; hegen1=hegen1+finhgt;}  //
                else if (egcen>iegen1 && egcen<iegen2){negen2++; neges2++; hegen2=hegen2+finhgt;}
                else if (egcen>iegen2 && egcen<iegen3){negen3++; neges3++; hegen3=hegen3+finhgt;}
                else if (egcen>iegen3 && egcen<iegen4){negen4++; neges4++; hegen4=hegen4+finhgt;}
                else if (egcen>iegen4 && egcen<iegen5){negen5++; neges5++; hegen5=hegen5+finhgt;}
                else {negen6++; neges6++; hegen6=hegen6+finhgt;}
            }
        }
        else if (nmuc > 1 && nelec == 0){      //  Mu cluster with >1 muons & no electrons
            sclid = 2;  nmucls ++;  nmuclt ++;
            // Mean distance of the clusters in the shower
            rmclsq = rclmn*rclmn;  rmcl    = sqrt(rmclsq);
            rmmcls = rmmcls + (rmcl-rmmcls)/nmucls;
            rmmclt = rmmclt + (rmcl-rmmclt)/nmuclt;
            //
            mmucl[jr][jt][ja]++;
            vnmucr[jr]++;  vnmuct[jt]++;  vnmuca[ja]++; }

        else if (nelec>0 && nmuc>0) {           // mixed cluster with e's & mu's
            sclid = 3;  nmxcls++;   nmxclt++;
            // Mean distance of the clusters in the shower
            rxclsq  = rclmn*rclmn;  rxcl    = sqrt(rxclsq);
            rmxcls = rmxcls + (rxcl-rmxcls)/nmxcls;
            rmxclt = rmxclt + (rxcl-rmxclt)/nmxclt;
            //
            mmxcl[jr][jt][ja]++;
            vnmxcr[jr] ++;  vnmxct[jt] ++;  vnmxca[ja] ++;
        }  // end mixed cluster
        else{notcls ++; notclt ++;   // other clusters
            sclid = 6;
            roclsq  = rclmn*rclmn;
            rocl    = sqrt(roclsq);
            rmocls  = rmocls + (rocl-rmocls)/notcls;
            rmoclt  = rmoclt + (rocl-rmoclt)/notclt;
        }
        //
        if (cmult>1 && ishowp < mxshop && iclusp < mxclup){ // prints mxclup clusters
            iclusp++;
            // cout << " - new cluster: ishow, ifp, clid, sclid, nmuc: " <<"\t"<< ishow <<"\t"<< ifp <<"\t"<< clid <<"\t"<< sclid <<"\t"<< nmuc << endl;
        // Cluster summary
        file1 << mpcr <<"\t"<< setprecision(4)<< lpcren   <<"\t"<<
            finhgt <<"\t"<< setprecision(0)<< ishow+1 <<"\t"<<
            nclust <<"\t"<< clid  <<"\t"<< sclid <<"\t"<<
            cmult  <<"\t"<< ngamc <<"\t"<< nelec <<"\t"<< nmuc  <<"\t"<<
            setprecision(4) << //<< fixed << setprecision(3)
            xclmn  <<"\t"<< yclmn <<"\t"<< rclmn <<"\t"<< sigrc <<"\t"<<
            tclmn  <<"\t"<< cltdt <<"\t"<< sigtc <<"\t"<<
            picfp  <<"\t"<< xfpc  <<"\t"<< yfpc  <<"\t"<< tfpc  <<"\t"<<
            zhfp   <<"\t"<< azfp  <<"\t"<< pfpc  <<"\t"<<
            piclp  <<"\t"<< xlpc  <<"\t"<< ylpc  <<"\t"<< tlpc  <<"\t"<<
            zhlp   <<"\t"<< azlp  <<"\t"<< plpc  <<"\t"<<
            eclmlt <<"\t"<< eclen <<"\t"<< ecldt <<"\t"<< eclth  <<"\t"<< eclsz <<"\t"<<
            ecled  <<"\t"<< eclem <<"\t"<< eclet <<"\t"<< eclest <<"\t"<<
            egcmlt <<"\t"<< egcen <<"\t"<< egcdt <<"\t"<< egcth  <<"\t"<< egcsz <<"\t"<<
            egced  <<"\t"<< egcem <<"\t"<< egcet <<"\t"<< egcest <<"\t"<<
            gmcmlt <<"\t"<< gmcen <<"\t"<< gmcdt <<"\t"<< gmcth  <<"\t"<< gmcsz <<"\t"<<
            gmced  <<"\t"<< gmcem <<"\t"<< gmcet <<"\t"<< gmcest <<"\t"<<
            endl;
        } // end if print
        cmult=0;    clid=0;     sclid=0;  sigrc=0;    sigtc=0.; ene=0;
        ngamc=0;    nelec=0;    nmuc=0;   cltsz=0;
        eclmlt=0;   eclen=0;    ecldt=0;  eclth=0;
        eclsz=0;    ecled=0;    eclem=0;  eclet=0;    eclest=0;
        pxelc =0;   pyelc=0;    pzelc=0;
        egcmlt=0;   egcen=0;    egcdt=0;  egcth=0;
        egcsz=0;    egced=0;    egcem=0;  egcet=0;    egcest=0;
        gmcmlt=0;   gmcen=0;    gmcdt=0;  gmcth=0;
        gmcsz=0;    gmced=0;    gmcem=0;  gmcet=0;    gmcest=0;
    }  // ends ifp loop

    for(Int_t i=0; i < 100000; i++){
       itag[i]=-1;      arr0[i]=-10.;   arr1[i]=-10.;   arr2[i]=-10.;   arr3[i]=-10.;
       arr4[i]=-10.;    arr5[i]=-10.;   arr6[i]=-10.;}
    //*
    if (ishowp < mxshop && ncluss > 0){      // prints mxshop showers
        ishowp ++;
        iclusp = 0;
        lemin  = log10(pcrene) - denes/2;
        lemax  = log10(pcrene) + denes/2;
        emin   = pow(10,lemin);
        emax   = pow(10,lemax);
        flints = (f0/gm1) * (pow(emin,-gm1) - pow(emax,-gm1));
        lflins = log10(flints);

        // file 2: Shower summary
        file2  << mpcr <<   "\t"
        << setprecision(4)
        << lpcren  <<"\t"<< finhgt <<"\t"<< lflins  <<"\t"
        << setprecision(0)
        << ishow+1<<"\t"
        << nsecs  <<"\t"<< ncluss  <<"\t"<< ngams   <<"\t"<< nels    <<"\t"
        << nmus   <<"\t"<< nns     <<"\t"<< nps     <<"\t"<< nops    <<"\t"
        << nelcls <<"\t"<< nmucls  <<"\t"<< nmxcls  <<"\t"<< negcls  <<"\t"
        << ngmcls <<"\t"<< notcls  <<"\t"
        << setprecision(4)
        << elmens <<"\t"<< mumens  <<"\t"
        << rmgms  <<"\t"<< rmels  <<"\t"<< rmmus   <<"\t" << rmnts   <<"\t" << rmpts  <<"\t"
        << rmots  <<"\t"
        << rmecls <<"\t"<< rmmcls <<"\t"<< rmxcls  <<"\t"
        << rmegcs <<"\t"<< rmgcls <<"\t"<< rmocls  <<"\t"
        << setprecision(0)
        << ngmes1 <<"\t"<< ngmes2 <<"\t"<< ngmes3  <<"\t"<< ngmes4 <<"\t"
        << ngmes5 <<"\t"<< ngmes6 <<"\t"
        << neles1 <<"\t"<< neles2 <<"\t"<< neles3  <<"\t"<< neles4 <<"\t"
        << neles5 <<"\t"<< neles6 <<"\t"
        << nmues1 <<"\t"<< nmues2 <<"\t"<< nmues3  <<"\t"<< nmues4 <<"\t"
        << nmues5 <<"\t"<< nmues6 <<"\t"
        << necms1 <<"\t"<< necms2 <<"\t"<< necms3  <<"\t"<< necms4 <<"\t"
        << necms5 <<"\t"<< necms6 <<"\t"
        << neces1 <<"\t"<< neces2 <<"\t"<< neces3  <<"\t"<< neces4 <<"\t"
        << neces5 <<"\t"<< neces6 <<"\t"
        << nedts1 <<"\t"<< nedts2 <<"\t"<< nedts3  <<"\t"<< nedts4 <<"\t"
        << nedts5 <<"\t"<< nedts6 <<"\t"
        << neths1 <<"\t"<< neths2 <<"\t"<< neths3  <<"\t"<< neths4 <<"\t"
        << neths5 <<"\t"<< neths6 <<"\t"
        << neszs1 <<"\t"<< neszs2 <<"\t"<< neszs3  <<"\t"<< neszs4 <<"\t"
        << neszs5 <<"\t"<< neszs6 <<"\t"
        << neges1 <<"\t"<< neges2 <<"\t"<< neges3  <<"\t"<< neges4 <<"\t"
        << neges5 <<"\t"<< neges6 <<"\t"
        << ngces1 <<"\t"<< ngces2 <<"\t"<< ngces3  <<"\t"<< ngces4 <<"\t"
        << ngces5 <<"\t"<< ngces6 <<"\t"
        << nelrs1 <<"\t"<< nelrs2 <<"\t"<< nelrs3  <<"\t"<< nelrs4 <<"\t"
        << nelrs5 <<"\t"<< nelrs6 <<"\t"
        << nmurs1 <<"\t"<< nmurs2 <<"\t"<< nmurs3  <<"\t"<< nmurs4 <<"\t"
        << nmurs5 <<"\t"<< nmurs6  <<"\t"
        << necrs1 <<"\t"<< necrs2 <<"\t"<< necrs3  <<"\t"<< necrs4 <<"\t"
        << necrs5 <<"\t"<< necrs6 <<"\t"
        << endl;}
    //*/
    // cout << " ncluss " << ncluss << endl;
    // cout << " --- ishow, nels, nmus, nelcls  :"<<"\t"<< ishow <<"\t"<< nels <<"\t"<< nmus <<"\t"<< nelcls <<"\t" << endl;
    //
        ngams = 0;  nels = 0;   nmus = 0;   nns  = 0;   nps  = 0;   nops = 0;
        ncluss = 0; nelcls= 0;  ngmcls=0;   nmucls=0;   nmxcls=0;   notcls = 0;  negcls = 0;
        rmgms = 0.;rmels = 0.; rmmus = 0.; rmnts= 0.;   rmpts = 0.; rmots = 0.;
        elmens = 0; mumens=0;
        rmecls=0;   rmmcls=0;   rmxcls=0;   rmegcs=0;   rmgcls=0;   rmocls=0;
        t0r1  = 0;  t0r2 = 0;   t0r3 = 0;   t0r4 = 0;   t0r5 = 0;   t0r6 = 0;
        pidfr1=-1;  pidfr2=-1;  pidfr3=-1;  pidfr4=-1;  pidfr5=-1;  pidfr6=-1;
        //
        ngmes1=0; ngmes2=0; ngmes3=0; ngmes4=0; ngmes5=0; ngmes6=0;
        neles1=0; neles2=0; neles3=0; neles4=0; neles5=0; neles6=0;
        neces1=0; neces2=0; neces3=0; neces4=0; neces5=0; neces6=0;
        necms1=0; necms2=0; necms3=0; necms4=0; necms5=0; necms6=0;
        nedts1=0; nedts2=0; nedts3=0; nedts4=0; nedts5=0; nedts6=0;
        neths1=0; neths2=0; neths3=0; neths4=0; neths5=0; neths6=0;
        neszs1=0; neszs2=0; neszs3=0; neszs4=0; neszs5=0; neszs6=0;
        needs1=0; needs2=0; needs3=0; needs4=0; needs5=0; needs6=0;
        neems1=0; neems2=0; neems3=0; neems4=0; neems5=0; neems6=0;
        neets1=0; neets2=0; neets3=0; neets4=0; neets5=0; neets6=0;
        neges1=0; neges2=0; neges3=0; neges4=0; neges5=0; neges6=0;
        ngces1=0; ngces2=0; ngces3=0; ngces4=0; ngces5=0; ngces6=0;
        nmues1=0; nmues2=0; nmues3=0; nmues4=0; nmues5=0; nmues6=0;
        // radial slices
        nelrs1=0; nelrs2=0; nelrs3=0; nelrs4=0; nelrs5=0; nelrs6=0;
        necrs1=0; necrs2=0; necrs3=0; necrs4=0; necrs5=0; necrs6=0;
        nmurs1=0; nmurs2=0; nmurs3=0; nmurs4=0; nmurs5=0; nmurs6=0;
        emers1=0; emers2=0; emers3=0; emers4=0; emers5=0; emers6=0;  // electron mean energy
        emtrs1=0; emtrs2=0; emtrs3=0; emtrs4=0; emtrs5=0; emtrs6=0;  // electron mean time
        mmers1=0; mmers2=0; mmers3=0; mmers4=0; mmers5=0; mmers6=0;  // muon mean energy
        mmtrs1=0; mmtrs2=0; mmtrs3=0; mmtrs4=0; mmtrs5=0; mmtrs6=0;  // muon mean time
        emdts1=0; emdts2=0; emdts3=0; emdts4=0; emdts5=0; emdts6=0;  // electron-muon dtime
} // end loop in  all showers
// ------------------------------------------------------------------------------------
// mean values
elment = elment/(nelt + 0.001);
mument = mument/(nmut + 0.001);
// Radial mean values
elmer1 = elmer1/(nelr1+0.01);  elmer2 = elmer2/(nelr2+0.01);
elmer3 = elmer3/(nelr3+0.01);  elmer4 = elmer4/(nelr4+0.01);
elmer5 = elmer5/(nelr5+0.01);  elmer6 = elmer6/(nelr6+0.01);
elmtr1 = elmtr1/(nelr1+0.01);  elmtr2 = elmtr2/(nelr2+0.01);
elmtr3 = elmtr3/(nelr3+0.01);  elmtr4 = elmtr4/(nelr4+0.01);
elmtr5 = elmtr5/(nelr5+0.01);  elmtr6 = elmtr6/(nelr6+0.01);
clmer1 = clmer1/(neclr1+0.01); clmer2 = clmer2/(neclr2+0.01);
clmer3 = clmer3/(neclr3+0.01); clmer4 = clmer4/(neclr4+0.01);
clmer5 = clmer5/(neclr5+0.01); clmer6 = clmer6/(neclr6+0.01);
clmtr1 = clmtr1/(neclr1+0.01); clmtr2 = clmtr2/(neclr2+0.01);
clmtr3 = clmtr3/(neclr3+0.01); clmtr4 = clmtr4/(neclr4+0.01);
clmtr5 = clmtr5/(neclr5+0.01); clmtr6 = clmtr6/(neclr6+0.01);
mumer1 = mumer1/(nmur1+0.01);  mumer2 = mumer2/(nmur2+0.01);
mumer3 = mumer3/(nmur3+0.01);  mumer4 = mumer4/(nmur4+0.01);
mumer5 = mumer5/(nmur5+0.01);  mumer6 = mumer6/(nmur6+0.01);
mumtr1 = mumtr1/(nmur1+0.01);  mumtr2 = mumtr2/(nmur2+0.01);
mumtr3 = mumtr3/(nmur3+0.01);  mumtr4 = mumtr4/(nmur4+0.01);
mumtr5 = mumtr5/(nmur5+0.01);  mumtr6 = mumtr6/(nmur6+0.01);
    
//
emdtr1 = elmer1 - mumer1;     emdtr2 = elmer2 - mumer2;
emdtr3 = elmer3 - mumer3;     emdtr4 = elmer4 - mumer4;
emdtr5 = elmer5 - mumer5;     emdtr6 = elmer6 - mumer6;
//
if(mumer1!= 0){emrat1 = elmer1/mumer1;}  //  mhmer1 = hmuen1/nmuen1;}
if(mumer2!= 0){emrat2 = elmer2/mumer2;}  // mhmer2 = hmuen2/nmuen2;}
if(mumer3!= 0){emrat3 = elmer3/mumer3;}  // mhmer3 = hmuen3/nmuen3;}
if(mumer4!= 0){emrat4 = elmer4/mumer4;}  // mhmer4 = hmuen4/nmuen4;}
if(mumer5!= 0){emrat5 = elmer5/mumer5;}  // mhmer5 = hmuen5/nmuen5;}
if(mumer6!= 0){emrat6 = elmer6/mumer6;}  // mhmer6 = hmuen6/nmuen6;}
// ***
cout << endl;
cout << "= Electrons radial distribution: " <<"\t"<< nelr1 <<"\t"<< nelr2 <<"\t"<< nelr3 <<"\t"<< nelr4<<"\t"<< nelr5 <<"\t"<< nelr6 << endl;
cout << "= Electron clusters energy:       " <<"\t"<< necen1 <<"\t"<<necen2 <<"\t"<<necen3 <<"\t"<<necen4<<"\t"<< necen5<<"\t"<< necen6 << endl;
cout << "= Electron clusters multiplicity: " <<"\t"<< neclm1 <<"\t" << neclm2 <<"\t"<< neclm3 <<"\t"<< neclm4 <<"\t"<< neclm5 <<"\t"<< neclm6 << endl;
cout << "= Electron clusters thermalicity: " <<"\t"<< necth1 <<"\t"<< necth2 <<"\t"<< necth3 <<"\t"<< necth4 <<"\t"<< necth5 <<"\t"<< necth6 << endl;
// ----------------------------------
if(emin==emax){
    emin = emin*(1-1/nshows); emax = emax*(1+1/(nshows)); }  // parche
flint = (f0/gm1) * (pow(emin,-gm1) - pow(emax,-gm1));
    
lflint = log10(flint);
//cout << "flint, emin, emax: " <<"\t"<< flint <<"\t"<< emin <<"\t"<< emax << endl;
//cout << endl;
cout << "* -------------------------" <<endl;
cout << "- Energy limits, FluxInt/m2·s·sr, LogFluxInt:" <<"\t"<<
       emin << " \t " << emax << " \t " << flint << " \t " << lflint <<endl;
cout << "- Height interval:                              " <<"\t"<<
           hmin << " \t " << hmax << endl;
cout << "* -------------------------" <<endl;
// ---------------------------------------------------------------------------
// file 3: Total summary
file3 << mpcr  <<"\t"<< setprecision(4)  << lpcrem   <<"\t"     << mnhgt   <<"\t"
    << setprecision(2) << flint*itsana/nshana          <<"\t"     << setprecision(0)
    << itsana <<"\t"<< ntsep  <<"\t"<< nclust  <<"\t"
    << ngamt  <<"\t"<< nelt   <<"\t"<< nmut    <<"\t"
    << nnt    <<"\t"<< npt    <<"\t"<< nopt    <<"\t"
    << nelclt <<"\t"<< nmuclt <<"\t"<< nmxclt  <<"\t"<< negclt  <<"\t"
    << ngmclt <<"\t"<< notclt <<"\t"
    << setprecision(4)
    << rmgmt  <<"\t"<< rmelt  <<"\t"<< rmmut  <<"\t"<< rmntt  <<"\t" << rmptt  <<"\t"
    << rmott  <<"\t"
    << rmeclt <<"\t"<< rmmclt <<"\t"<< rmxclt <<"\t"
    << rmegct <<"\t"<< rmgclt <<"\t"<< rmoclt <<"\t"
    << setprecision(0)
    << ngmen1 <<"\t"<< ngmen2 <<"\t"<< ngmen3 <<"\t"<< ngmen4 <<"\t"
    << ngmen5 <<"\t"<< ngmen6 <<"\t"
    << nelen1 <<"\t"<< nelen2 <<"\t"<< nelen3 <<"\t"<< nelen4 <<"\t"
    << nelen5 <<"\t"<< nelen6 <<"\t"
    << nmuen1 <<"\t"<< nmuen2 <<"\t"<< nmuen3 <<"\t"<< nmuen4 <<"\t"
    << nmuen5 <<"\t"<< nmuen6 <<"\t"
    << nnten1 <<"\t"<< nnten2 <<"\t"<< nnten3 <<"\t"<< nnten4 <<"\t"
    << nnten5 <<"\t"<< nnten6 <<"\t"
    << npren1 <<"\t"<< npren2 <<"\t"<< npren3 <<"\t"<< npren4 <<"\t"
    << npren5 <<"\t"<< npren6 <<"\t"
    << neclm1 <<"\t"<< neclm2 <<"\t"<< neclm3 <<"\t"<< neclm4 <<"\t"
    << neclm5 <<"\t"<< neclm6 <<"\t"
    << necen1 <<"\t"<< necen2 <<"\t"<< necen3 <<"\t"<< necen4 <<"\t"
    << necen5 <<"\t"<< necen6 <<"\t"
    << necdt1 <<"\t"<< necdt2 <<"\t"<< necdt3 <<"\t"<< necdt4 <<"\t"
    << necdt5 <<"\t"<< necdt6 <<"\t"
    << necth1 <<"\t"<< necth2 <<"\t"<< necth3 <<"\t"<< necth4 <<"\t"
    << necth5 <<"\t"<< necth6 <<"\t"
    << necsz1 <<"\t"<< necsz2 <<"\t"<< necsz3 <<"\t"<< necsz4 <<"\t"
    << necsz5 <<"\t"<< necsz6 <<"\t"
    << neced1 <<"\t"<< neced2 <<"\t"<< neced3 <<"\t"<< neced4 <<"\t"
    << neced5 <<"\t"<< neced6 <<"\t"
    << necem1 <<"\t"<< necem2 <<"\t"<< necem3 <<"\t"<< necem4 <<"\t"
    << necem5 <<"\t"<< necem6 <<"\t"
    << necet1 <<"\t"<< necet2 <<"\t"<< necet3 <<"\t"<< necet4 <<"\t"
    << necet5 <<"\t"<< necet6 <<"\t"
    << neest1 <<"\t"<< neest2 <<"\t"<< neest3 <<"\t"<< neest4 <<"\t"
    << neest5 <<"\t"<< neest6 <<"\t"
    << negen1 <<"\t"<< negen2 <<"\t"<< negen3 <<"\t" << negen4 <<"\t"
    << negen5 <<"\t"<< negen6 <<"\t"
    << ngcen1 <<"\t"<< ngcen2 <<"\t"<< ngcen3 <<"\t" << ngcen4 <<"\t"
    << ngcen5 <<"\t"<< ngcen6 <<"\t"
    << endl;
// file 4: Radial summary
file4 << mpcr    <<"\t"<< setprecision(4)   << lpcrem    <<"\t"<< mnhgt     <<"\t"
    << setprecision(2)   << flint*itsana/nshana            <<"\t"<< setprecision(0)
    << itsana <<"\t"<< ntsep  <<"\t"<< nclust <<"\t"
    << nelr1  <<"\t"<< nelr2  <<"\t"<< nelr3  <<"\t"<< nelr4   <<"\t"
    << nelr5  <<"\t"<< nelr6  <<"\t"
    << nmur1  <<"\t"<< nmur2  <<"\t"<< nmur3  <<"\t"<< nmur4   <<"\t"
    << nmur5  <<"\t"<< nmur6  <<"\t"
    << neclr1 <<"\t"<< neclr2 <<"\t"<< neclr3 <<"\t"<< neclr4  <<"\t"
    << neclr5 <<"\t"<< neclr6 <<"\t"
    << vnmucr[0] <<"\t"<< vnmucr[1] <<"\t"<< vnmucr[2] <<"\t"<< vnmucr[3] <<"\t"
    << vnmucr[4] <<"\t"<< vnmucr[5] <<"\t"
    << vnmxcr[0] <<"\t"<< vnmxcr[1] <<"\t"<< vnmxcr[2] <<"\t"<< vnmxcr[3] <<"\t"
    << vnmxcr[4] <<"\t"<< vnmxcr[5] <<"\t"
    << setprecision(4)
    << elmer1 <<"\t"<< elmer2 <<"\t"<< elmer3 <<"\t"<< elmer4 <<"\t"
    << elmer5 <<"\t"<< elmer6 <<"\t"
    << elmtr1 <<"\t"<< elmtr2 <<"\t"<< elmtr3 <<"\t"<< elmtr4 <<"\t"
    << elmtr5 <<"\t"<< elmtr6 <<"\t"
    << mumer1 <<"\t"<< mumer2 <<"\t"<< mumer3 <<"\t"<< mumer4 <<"\t"
    << mumer5 <<"\t"<< mumer6 <<"\t"
    << mumtr1 <<"\t"<< mumtr2 <<"\t"<< mumtr3 <<"\t"<< mumtr4 <<"\t"
    << mumtr5 <<"\t"<< mumtr6 <<"\t"
    << clmer1 <<"\t"<< clmer2 <<"\t"<< clmer3 <<"\t"<< clmer4 <<"\t"
    << clmer5 <<"\t"<< clmer6 <<"\t"
    << clmtr1 <<"\t"<< clmtr2 <<"\t"<< clmtr3 <<"\t"<< clmtr4 <<"\t"
    << clmtr5 <<"\t"<< clmtr6 <<"\t"
    << emdtr1 <<"\t"<< emdtr2 <<"\t"<< emdtr3  <<"\t"<< emdtr4 <<"\t"
    << emdtr5 <<"\t"<< emdtr6 <<"\t"
    << emrat1 <<"\t"<< emrat2 <<"\t"<< emrat3  <<"\t"<< emrat4 <<"\t"
    << emrat5 <<"\t"<< emrat6 <<"\t"
    << endl;
// << endl;
cout << "- Radial slices/m:   " <<"\t"<<
    "-0-"<<"\t"<<"-100-"<<"\t"<<"-200-"<<"\t"<<"-500-"<<"\t"<<"-1000-"<<"\t"<<"-4000- " << endl;
cout << "emc R Distribution: " <<"\t"<<
    vnemcr[0] <<"\t"<< vnemcr[1] <<"\t"<< vnemcr[2] <<"\t"<< vnemcr[3] <<"\t"<<
    vnemcr[4] <<"\t"<< vnemcr[5] << endl;
cout << "muc R Distribution: " <<"\t"<<
    vnmucr[0] <<"\t"<< vnmucr[1] <<"\t"<< vnmucr[2] <<"\t"<< vnmucr[3] <<"\t"<<
    vnmucr[4] <<"\t"<< vnmucr[5] << endl;
cout << "mxc R Distribution: " <<"\t"<<
    vnmxcr[0] <<"\t"<< vnmxcr[1] <<"\t"<< vnmxcr[2] <<"\t"<< vnmxcr[3] <<"\t"<<
    vnmxcr[4] <<"\t"<< vnmxcr[5] << endl;
    cout << endl;
    
cout << "- dTime intervals/ns:" <<"\t"<<
    "-0- "<<"\t"<<" -0.5- "<<"\t"<<" -2- "<<"\t"<<" -5- "<<"\t"<<" -20- "<<"\t"<<" -50- "<< endl;
cout << "emc dT distribution: " <<  "\t"<<
   vnemct[0] <<"\t"<< vnemct[1] <<"\t"<< vnemct[2] <<"\t"<<
   vnemct[3] <<"\t"<< vnemct[4] <<"\t"<< vnemct[5] << endl;
   cout << "muc dT distribution: " <<"\t"<<
   vnmuct[0] <<"\t"<< vnmuct[1] <<"\t"<< vnmuct[2] <<"\t"<<
   vnmuct[3] <<"\t"<< vnmuct[4] <<"\t"<< vnmuct[5] << endl;
   cout << "mxc dT distribution: " <<"\t"<<
   vnmxct[0] <<"\t"<< vnmxct[1] <<"\t"<< vnmxct[2] <<"\t"<<
   vnmxct[3] <<"\t"<< vnmxct[4] <<"\t"<< vnmxct[5] << endl;
    cout << endl;

cout << "--- N. Events by EClDens/dT: "  <<"\t"   <<"\t"<<
         neest1 <<"\t"<< neest2 <<"\t"<< neest3   <<"\t"<<
         neest4 <<"\t"<< neest5 <<"\t"<< neest6 << endl;
//    cout << "-- gamp, ngmen1, hgmen1: " <<"\t" << pmod <<"\t" << ngmen1 <<"\t" << hgmen1 << endl;
mhgme1 = hgmen1/ngmen1; mhgme2 = hgmen2/ngmen2; mhgme3 = hgmen3/ngmen3;
mhgme4 = hgmen4/ngmen4; mhgme5 = hgmen5/ngmen5; mhgme6 = hgmen6/ngmen6;
mhele1 = helen1/nelen1; mhele2 = helen2/nelen2; mhele3 = helen3/nelen3;
mhele4 = helen4/nelen4; mhele5 = helen5/nelen5; mhele6 = helen6/nelen6;
mhmue1 = hmuen1/nmuen1; mhmue2 = hmuen2/nmuen2; mhmue3 = hmuen3/nmuen3;
mhmue4 = hmuen4/nmuen4; mhmue5 = hmuen5/nmuen5; mhmue6 = hmuen6/nmuen6;
mhnte1 = hnten1/nnten1; mhnte2 = hnten2/nnten2; mhnte3 = hnten3/nnten3;
mhnte4 = hnten4/nnten4; mhnte5 = hnten5/nnten5; mhnte6 = hnten6/nnten6;
mhpre1 = hpren1/npren1; mhpre2 = hpren2/npren2; mhpre3 = hpren3/npren3;
mhpre4 = hpren4/npren4; mhpre5 = hpren5/npren5; mhpre6 = hpren6/npren6;
mhecm1 = heclm1/neclm1; mhecm2 = heclm2/neclm2; mhecm3 = heclm3/neclm3;
mhecm4 = heclm4/neclm4; mhecm5 = heclm5/neclm5; mhecm6 = heclm6/neclm6;
mhece1 = hecen1/necen1; mhece2 = hecen2/necen2; mhece3 = hecen3/necen3;
mhece4 = hecen4/necen4; mhece5 = hecen5/necen5; mhece6 = hecen6/necen6;
mrece1 = recen1/necen1; mrece2 = recen2/necen2; mrece3 = recen3/necen3;
mrece4 = recen4/necen4; mrece5 = recen5/necen5; mrece6 = recen6/necen6;
mtece1 = tecen1/necen1; mtece2 = tecen2/necen2; mtece3 = tecen3/necen3;
mtece4 = tecen4/necen4; mtece5 = tecen5/necen5; mtece6 = tecen6/necen6;
mhedt1 = hecdt1/necdt1; mhedt2 = hecdt2/necdt2; mhedt3 = hecdt3/necdt3;
mhedt4 = hecdt4/necdt4; mhedt5 = hecdt5/necdt5; mhedt6 = hecdt6/necdt6;
mheth1 = hecth1/necth1; mheth2 = hecth2/necth2; mheth3 = hecth3/necth3;
mheth4 = hecth4/necth4; mheth5 = hecth5/necth5; mheth6 = hecth6/necth6;
mhcsz1 = hecsz1/necsz1; mhcsz2 = hecsz2/necsz2; mhcsz3 = hecsz3/necsz3;
mhcsz4 = hecsz4/necsz4; mhcsz5 = hecsz5/necsz5; mhcsz6 = hecsz6/necsz6;
mheed1 = heced1/neced1; mheed2 = heced2/neced2; mheed3 = heced3/neced3;
mheed4 = heced4/neced4; mheed5 = heced5/neced5; mheed6 = heced6/neced6;
mheem1 = hecem1/necem1; mheem2 = hecem2/necem2; mheem3 = hecem3/necem3;
mheem4 = hecem4/necem4; mheem5 = hecem5/necem5; mheem6 = hecem6/necem6;
mheet1 = hecet1/necet1; mheet2 = hecet2/necet2; mheet3 = hecet3/necet3;
mheet4 = hecet4/necet4; mheet5 = hecet5/necet5; mheet6 = hecet6/necet6;
mhets1 = heets1/neest1; mhets2 = heets2/neest2; mhets3 = heets3/neest3;
mhets4 = heets4/neest4; mhets5 = heets5/neest5; mhets6 = heets6/neest6;
mhege1 = hegen1/negen1; mhege2 = hegen2/negen2; mhege3 = hegen3/negen3;
mhege4 = hegen4/negen4; mhege5 = hegen5/negen5; mhege6 = hegen6/negen6;
mhgce1 = hgcen1/ngcen1; mhgce2 = hgcen2/ngcen2; mhgce3 = hgcen3/ngcen3;
mhgce4 = hgcen4/ngcen4; mhgce5 = hgcen5/ngcen5; mhgce6 = hgcen6/ngcen6;
//
file5   << mpcr    <<"\t"<< setprecision(4) << lpcrem <<"\t"<< mnhgt <<"\t"
        << setprecision(2) << flint*itsana/nshana       <<"\t"<< setprecision(0)
        << itsana <<"\t"<< ntsep  <<"\t"<< nclust <<"\t"
        << ngmen1 <<"\t"<< ngmen2 <<"\t"<< ngmen3 <<"\t"
        << ngmen4 <<"\t"<< ngmen5 <<"\t"<< ngmen6 <<"\t"
        << nelen1 <<"\t"<< nelen2 <<"\t"<< nelen3 <<"\t"
        << nelen4 <<"\t"<< nelen5 <<"\t"<< nelen6 <<"\t"
        << nmuen1 <<"\t"<< nmuen2 <<"\t"<< nmuen3 <<"\t"
        << nmuen4 <<"\t"<< nmuen5 <<"\t"<< nmuen6 <<"\t"
        << nnten1 <<"\t"<< nnten2 <<"\t"<< nnten3 <<"\t"
        << nnten4 <<"\t"<< nnten5 <<"\t"<< nnten6 <<"\t"
        << npren1 <<"\t"<< npren2 <<"\t"<< npren3 <<"\t"
        << npren4 <<"\t"<< npren5 <<"\t"<< npren6 <<"\t"
        << neclm1 <<"\t"<< neclm2 <<"\t"<< neclm3 <<"\t"
        << neclm4 <<"\t"<< neclm5 <<"\t"<< neclm6 <<"\t"
        << necen1 <<"\t"<< necen2 <<"\t"<< necen3 <<"\t"
        << necen4 <<"\t"<< necen5 <<"\t"<< necen6 <<"\t"
        << necdt1 <<"\t"<< necdt2 <<"\t"<< necdt3 <<"\t"
        << necdt4 <<"\t"<< necdt5 <<"\t"<< necdt6 <<"\t"
        << necth1 <<"\t"<< necth2 <<"\t"<< necth3 <<"\t"
        << necth4 <<"\t"<< necth5 <<"\t"<< necth6 <<"\t"
        << necsz1 <<"\t"<< necsz2 <<"\t"<< necsz3 <<"\t"
        << necsz4 <<"\t"<< necsz5 <<"\t"<< necsz6 <<"\t"
        << neced1 <<"\t"<< neced2 <<"\t"<< neced3 <<"\t"
        << neced4 <<"\t"<< neced5 <<"\t"<< neced6 <<"\t"
        << necem1 <<"\t"<< necem2 <<"\t"<< necem3 <<"\t"
        << necem4 <<"\t"<< necem5 <<"\t"<< necem6 <<"\t"
        << necet1 <<"\t"<< necet2 <<"\t"<< necet3 <<"\t"
        << necet4 <<"\t"<< necet5 <<"\t"<< necet6 <<"\t"
        << neest1 <<"\t"<< neest2 <<"\t"<< neest3 <<"\t"
        << neest4 <<"\t"<< neest5 <<"\t"<< neest6 <<"\t"
        << negen1 <<"\t"<< negen2 <<"\t"<< negen3 <<"\t"
        << negen4 <<"\t"<< negen5 <<"\t"<< negen6 <<"\t"
        << ngcen1 <<"\t"<< ngcen2 <<"\t"<< ngcen3 <<"\t"
        << ngcen4 <<"\t"<< ngcen5 <<"\t"<< ngcen6 <<"\t"
        << setprecision(4)
        << mhgme1 <<"\t"<< mhgme2 <<"\t"<< mhgme3 <<"\t"
        << mhgme4 <<"\t"<< mhgme5 <<"\t"<< mhgme6 <<"\t"
        << mhele1 <<"\t"<< mhele2 <<"\t"<< mhele3 <<"\t"
        << mhele4 <<"\t"<< mhele5 <<"\t"<< mhele6 <<"\t"
        << mhmue1 <<"\t"<< mhmue2 <<"\t"<< mhmue3 <<"\t"
        << mhmue4 <<"\t"<< mhmue5 <<"\t"<< mhmue6 <<"\t"
        << mhnte1 <<"\t"<< mhnte2 <<"\t"<< mhnte3 <<"\t"
        << mhnte4 <<"\t"<< mhnte5 <<"\t"<< mhnte6 <<"\t"
        << mhpre1 <<"\t"<< mhpre2 <<"\t"<< mhpre3 <<"\t"
        << mhpre4 <<"\t"<< mhpre5 <<"\t"<< mhpre6 <<"\t"
        << mhecm1 <<"\t"<< mhecm2 <<"\t"<< mhecm3 <<"\t"
        << mhecm4 <<"\t"<< mhecm5 <<"\t"<< mhecm6 <<"\t"
        << mhece1 <<"\t"<< mhece2 <<"\t"<< mhece3 <<"\t"
        << mhece4 <<"\t"<< mhece5 <<"\t"<< mhece6 <<"\t"
        << mrece1 <<"\t"<< mrece2 <<"\t"<< mrece3 <<"\t"
        << mrece4 <<"\t"<< mrece5 <<"\t"<< mrece6 <<"\t"
        << mtece1 <<"\t"<< mtece2 <<"\t"<< mtece3 <<"\t"
        << mtece4 <<"\t"<< mtece5 <<"\t"<< mtece6 <<"\t"
        << mhedt1 <<"\t"<< mhedt2 <<"\t"<< mhedt3 <<"\t"
        << mhedt4 <<"\t"<< mhedt5 <<"\t"<< mhedt6 <<"\t"
        << mheth1 <<"\t"<< mheth2 <<"\t"<< mheth3 <<"\t"
        << mheth4 <<"\t"<< mheth5 <<"\t"<< mheth6 <<"\t"
        << mhcsz1 <<"\t"<< mhcsz2 <<"\t"<< mhcsz3 <<"\t"
        << mhcsz4 <<"\t"<< mhcsz5 <<"\t"<< mhcsz6 <<"\t"
        << mheed1 <<"\t"<< mheed2 <<"\t"<< mheed3 <<"\t"
        << mheed4 <<"\t"<< mheed5 <<"\t"<< mheed6 <<"\t"
        << mheem1 <<"\t"<< mheem2 <<"\t"<< mheem3 <<"\t"
        << mheem4 <<"\t"<< mheem5 <<"\t"<< mheem6 <<"\t"
        << mheet1 <<"\t"<< mheet2 <<"\t"<< mheet3 <<"\t"
        << mheet4 <<"\t"<< mheet5 <<"\t"<< mheet6 <<"\t"
        << mhets1 <<"\t"<< mhets2 <<"\t"<< mhets3 <<"\t"
        << mhets4 <<"\t"<< mhets5 <<"\t"<< mhets6 <<"\t"
        << mhege1 <<"\t"<< mhege2 <<"\t"<< mhege3 <<"\t"
        << mhege4 <<"\t"<< mhege5 <<"\t"<< mhege6 <<"\t"
        << mhgce1 <<"\t"<< mhgce2 <<"\t"<< mhgce3 <<"\t"
        << mhgce4 <<"\t"<< mhgce5 <<"\t"<< mhgce6
        << endl;
    
cout << "*** MHight Dif.  Muons - Electrons: "  <<"\t"<<
            mhmue1 - mhece1 <<"\t"<< mhmue2 - mhece2 <<"\t"<< mhmue3 - mhece3 <<"\t"<<
            mhmue4 - mhece4 <<"\t"<< mhmue5 - mhece5 <<"\t"<< mhmue6 - mhece6 << endl;
cout << "*** MHight Ratio Electrons/Muons "  <<"\t"<<
            mhece1/mhmue1 <<"\t"<< mhece2/mhmue2 <<"\t"<< mhece3/mhmue3 <<"\t"<<
            mhece4/mhmue4 <<"\t"<< mhece5/mhmue5 <<"\t"<< mhece6/mhmue6 << endl;
cout << endl;
cout << "- dAngle intervals/º:" <<"\t"<<
    "-0- "<<"\t"<<" -2- "<<"\t"<<" -4- "<<"\t"<<" -8- "<<"\t"<<" -20- "<<"\t"<<" -40- "<< endl;
cout << "emc dAng distribution: " <<"\t"<<
   vnemca[0] <<"\t"<< vnemca[1] <<"\t"<< vnemca[2] <<"\t"<<
   vnemca[3] <<"\t"<< vnemca[4] <<"\t"<< vnemca[5] << endl;
cout << "muc dAng distribution: " <<"\t"<<
   vnmuca[0] <<"\t"<< vnmuca[1] <<"\t"<< vnmuca[2] <<"\t"<<
   vnmuca[3] <<"\t"<< vnmuca[4] <<"\t"<< vnmuca[5] << endl;
cout << "mxc dAng distribution: " <<"\t"<<
   vnmxca[0] <<"\t"<< vnmxca[1] <<"\t"<< vnmxca[2] <<"\t"<<
   vnmxca[3] <<"\t"<< vnmxca[4] <<"\t"<< vnmxca[5] << endl;
cout << endl;
//cout << "---neceds: " <<"\t"<<
    //neced1 <<"\t"<< neced2 <<"\t"<< neced3 <<"\t"<<
    //neced4 <<"\t"<< neced5 <<"\t"<< neced6 << endl;
//cout << "---necems: " <<"\t"<<
    //necem1 <<"\t"<< necem2 <<"\t"<< necem3 <<"\t"<<
    //necem4 <<"\t"<< necem5 <<"\t"<< necem6 << endl;
//cout << "---necets: " <<"\t"<<
    //necet1 <<"\t"<< necet2 <<"\t"<< necet3 <<"\t"<<
    //necet4 <<"\t"<< necet5 <<"\t"<< necet6 << endl;
//cout << endl;
    
// Close all files
printf("******  Proceso completado  ******");
    file0.close();file1.close(); file2.close(); file3.close(); file4.close(); file5.close(); cout << endl;

    if(ntsep==0){cout << endl; cout << "*** ntsep = 0!!! " << endl;
        cout<< endl;}
    else{
    cout << "* Num Showers analizados: " <<"\t" <<"\t" << itsana <<"\t ( " << 100*itsana/nshana << " %)" << endl;
    cout << "* Parametros Iniciales. MxDist, MxSig:"<<"\t"<< mxdst <<"\t"<< sigrmx << endl;
    //cout << endl;
    cout << "* NTSec, NTSec/Shw, NTSec/km2·s·sr:" <<"\t"<< ntsep <<"\t"<< ntsep/itsana <<"\t"<< 1.E6 * ntsep * flint/itsana << endl;
    cout << "* Contador (gemnpo):       " <<"\t"<<"\t"<< ngamt <<"\t"<< nelt <<"\t"<< nmut <<  "\t" << nnt <<  "\t" << npt <<  "\t" << nopt << endl;
    cout << "* Distribucion/% (gemnpo):     " <<"\t"<<"\t"<< 100*ngamt/ntsep <<"\t"<< 100*nelt/ntsep <<"\t"<< 100*nmut/ntsep <<  "\t" << 100*nnt/ntsep <<  "\t" << 100*npt/ntsep <<  "\t" << 100*nopt/ntsep << endl;
    cout << "* NParticulas (gemnpo): " <<"\t"<<"\t"<< ngamt <<"\t"<< nelt <<"\t"<< nmut <<"\t"<< nnt <<  "\t"<< npt <<  "\t" << nopt  << endl;
    cout << "* NParticulas/1kShow (gemnpo):" <<"\t\t"<< 1000* ngamt/itsana <<"\t"<< 1000*nelt/itsana <<"\t"<< 1000*nmut/itsana <<"\t"<< 1000*nnt/itsana <<"\t"<< 1000*npt/itsana <<"\t"<< 1000*nopt/itsana << endl;
        cout << " --- " << endl;
        cout << "* Num. Clusters, NClt/1kShower:  " <<"\t"<< nclust <<"\t"<< 1000*nclust/itsana << endl;
        if(nclust > 0){
            cout << "* NClstId  : ECl, MuC, MxC, EGC, GmC, OtC: " <<"\t"<< nelclt <<"\t"<< nmuclt  <<"\t"<< nmxclt <<"\t" << negclt <<"\t" << ngmclt  <<"\t"<< notclt << endl;
            cout << "* NClstId/%: ECl, MuC, MxC, EGC, GmC, OtC: " <<"\t"<< 100*nelclt/nclust <<"\t"<< 100*nmuclt/nclust <<"\t"<< 100*nmxclt/nclust <<"\t"<< 100*negclt/nclust <<"\t"<< 100*ngmclt/nclust <<"\t"<< 100*notclt/nclust << endl;
    }
//
}
cout << endl;
//TH2D *h1 = new TH2D("h1", "cluster distribution", 100, 0., 1000., 100, 0., 1000.);
    
// Estimamos densidades de particulas y sigmas correspondientes
}

/*
 void bigpm() {
    TFile *f = new TFile("cernstaff.root");
    TTree *T = (TTree*)f->Get("T");
    TCanvas *c1 = new TCanvas("c1");
    T->Draw("Cost:Age>>hist","","goff");
    TH2F *h = new TH2F("h","Cost vs Age",60,10,70,20,0,20000);
    h->Draw();
    c1->Update();
    T->SetMarkerStyle(3);
    T->SetMarkerColor(3); T.Draw("Cost:Age","Grade==3","same");
    T->SetMarkerColor(4); T.Draw("Cost:Age","Grade==4","same");
    T->SetMarkerColor(5); T.Draw("Cost:Age","Grade==5","same");
    T->SetMarkerColor(6); T.Draw("Cost:Age","Grade==6","same");
    T->SetMarkerColor(1); T.Draw("Cost:Age","Grade==10","same");
    T->SetMarkerColor(2); T.Draw("Cost:Age","Grade==12","same");
 }
 */

