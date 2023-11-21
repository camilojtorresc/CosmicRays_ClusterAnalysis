import math
import numpy as np

# Primary cosmic ray properties
mpcr = 1       # A mass of primary cosmic ray: 1, 4, 12, 56
f0 = 1.8E4     # Intercepts H, He, C, Fe: 1.8E4, 5.0E3, 1.0E3, 1.5E3
spindex = 2.7  # Spectral index
zhmin = 0     # Zenith angle min 0, 25, 37, 46, 53
zhmax = 25    # Zenith angle max
dene = 0.25   # Width of the LogEnergy interval generated

# H:  8 Height Intervals: 0, 12, 16, 18, 20, 22, 28, 40, 150 km
# He: 8 Height Intervals: 0, 16, 20, 22, 24, 26, 30, 36, 150 km
# C:  8 Height Intervals: 0, 24, 30, 36, 40, 44, 48, 64, 150 km
# Fe: 8 Height Intervals: 0, 32, 40, 44, 48, 52, 56, 64, 150 km
hmin = 0.
hmax = 150.

nshana = 100000  # Max. N.Showers to analyze
mxshop = 50     # nshana/50; # 100; # Mx. N.Showers to print
mxclup = 50     # 10;  # Mx. N.Clusters/Shw to print
mxparp = 500    # 10; # Mx. N.Particles/Energy to print
eperiod = 0
enemod = 0      # shower period and module to print particles
itsana = 0
ishowp = 0
iclusp = 0      # iCounters

# Data Files
# Simulation with 4 samples
# nh: n alturas. an: num. de intervalo angular
nsamp = 1       # Number of samples per shower.
isami = 1       # N. Initial Sample
gm1 = spindex - 1  # pcr energy distribution decrease factor

# Secondary masses
mele = 0.000511   # GeV
mmu = 0.10566     # GeV
mgam = 0.0        # GeV
mp = 0.9383       # GeV
mn = 0.9396       # GeV
mother = 0.14     # GeV. we assume they are mainly pions

# Cluster size analyzed
dsize = 1.8       # m^2. Size of the detector

# Printing options
iprhd = 1         # Index for printing headers. 0=No Print.  1=Print

# Output file names
fout0 = "xmdatp_p_E30b_2m2.txt"  # Particle output
fout1 = "xmdatc_p_E30b_2m2.txt"  # Clusters output
fout2 = "xmdats_p_E30b_2m2.txt"  # Showers output
fout3 = "xmdatm_p_E30b_2m2.txt"  # Main Summary
fout4 = "xmdatr_p_E30b_2m2.txt"  # Cluster Radial distribution
fout5 = "xmdath_p_E30b_2m2.txt"  # Mean Height estimation


class Mclust:
    def __init__(self):
        self.fChain = None
        # ... Any other member variables initialization can go here

    def loop(self):

        if self.fChain == 0:
            return
        
        nshows = len(self.fChain)  # No. of showers in file. Assuming fChain is a list or similar.

        if iprhd == 1:
            print("\n\n***   Print Header\n\n\n")
        
        print("***************************** ")
        print("- nShowers in file:  ", nshows)
        if nshana > nshows:
            nshana = nshows
        print("- nShowers analyzed: ", nshana)
        
        denes = dene / nshana  # ~ Energy Interval Width associated to each Shower
        
        print("***************************** ")
        
        eperiod = nshana / mxparp

        nbytes, nb = 0, 0
        fbytes, fb = 0, 0
        ebytes, eb = 0, 0
        hbytes, hb = 0, 0
        tbytes, tb = 0, 0

        # identificadores corsika 1 gammas,  2 3      e- e+, 5 6      mu- mu+, 13 14    n p

        rad2g = 180 / math.pi

        mxdst = 2 * math.sqrt(dsize / math.pi)  # diameter of a dsize area
        sigrmx = mxdst / math.sqrt(2)

        tag = 0
        gpcut = 0.1    # gammas minimum momentum cut

        zero, iconta, icontb, icontc = 0, 0, 0, 0  # auxiliary counters for prints
        ncont, n, nsecs, isep, ntsep = 0, 0, 0, 0, 0

        icgam, icele, icmu, icpt, icnt = 1, 1000, 100000, 10000000, 100000000  # codes

        pid = None
        pidf, pidfr1, pidfr2, pidfr3, pidfr4, pidfr5, pidfr6 = None, -1, -1, -1, -1, -1, -1
        pic, pic1, pic2, picfp, piclp = 0, 0, 0, 0, 0  # id code of particles, first, last...
        clid, sclid = 0, 0  # ClusterId, ShortClusterId

        ngamt, nelt, nmut, nnt, npt, nopt = 0, 0, 0, 0, 0, 0  # N. of particles
        ngams, nels, nmus, nns, nps, nops = 0, 0, 0, 0, 0, 0  # N. of particles in shower
        ngamc, nelec, nmuc = 0, 0, 0  # N. of particles in cluster

        ncluss, nclust, cmult, cmulp1 = 0, 0, 0, 0
        nelcls, nmucls, nmxcls, negcls, ngmcls, notcls = 0, 0, 0, 0, 0, 0  # cluster count in showers
        nelclt, nmuclt, nmxclt, negclt, ngmclt, notclt = 0, 0, 0, 0, 0, 0  # total cluster count
        pcrene, emin, emax, lemin, lemax, flint, flints, lflint, lflins = 1.E20, 1.E20, 0., 0., 0., 0., 0., 0., 0.
        ene, gmten, elten = 0., 0., 0.
        elment, mument, elmens, mumens, emrats = 0., 0., 0., 0., 0.
        rx, ry, rr, dxc, dyc, dtc, edc, x1, x2, y1, y2, drc, rclmn = [0.] * 13
        xfpc, xlpc, yfpc, ylpc, drcsq, cltdt, clang, cltst, cltsz = [0.] * 9
        eclmlt, eclen, ecldt, eclth, eclsz, ecled, eclem, eclet, eclest = [0.] * 9
        egcmlt, egcen, egcdt, egcth, egcsz, egced, egcem, egcet, egcest = [0.] * 9
        gmcmlt, gmcen, gmcdt, gmcth, gmcsz, gmced, gmcem, gmcet, gmcest = [0.] * 9
        dt, dtsq, dxsq, dysq, sigrc = [0.] * 5
        pxelc, pyelc, pzelc, pxegc, pyegc, pzegc, pxgmc, pygmc, pzgmc = [0.] * 9
        rmgms, rmgmt, rmels, rmelt, rmmus, rmmut, rmnts, rmntt, rmpts, rmptt, rmots, rmott = [0.] * 12
        reclsq, recl, rmecls, rmeclt, rmclsq, rmcl, rmmcls, rmmclt, rxclsq, rxcl, rmxcls, rmxclt = [0.] * 12
        regcsq, regc, rmegcs, rmegct, rgclsq, rgcl, rmgcls, rmgclt, roclsq, rocl, rmocls, rmoclt = [0.] * 12
        zha, zen1, zen2, aza, azh1, azh2 = [0.] * 6
        t0, time, t1, t2, tfg, tlg, tfe, tle, tfm, tlm, tfpc, tlpc = [0.] * 12
        zhfp, zhlp, azfp, azlp, zhfe, zhle, zhfm, zhlm, azfe, azle, azfm, azlm = [0.] * 12
        px, py, pz, px1, py1, pz1, px2, py2, pz2, pmod, pm1, pm2, psq, pfpc, plpc = [0.] * 15
        xclmn, yclmn, tclmn, xcmno, ycmno, tcmno = [0.] * 6
        sigx, sigxsq, sigy, sigysq, sigtc, sigtsq = [0.] * 6
        ssx, ssy, sst, ssxo, ssyo, ssto = [0.] * 6
        nxfp, nyfp, nzfp, nxlp, nylp, nzlp = [0.] * 6
        finhgt, mnhgt = 0., 0.
        lpcren, lpcrem = 0., 0.
        jr, jt, ja = [0] * 3
        ir1, ir2, ir3, ir4, ir5 = 100, 200, 500, 1000, 2500
        # MuEnergy
        imuen = [2, 4, 6, 10, 20]
        # NeutEnergy
        inten = [2, 3, 5, 8, 15]
        # ProtEnergy
        ipren = [2, 3, 5, 8, 15]
        # ECl.Mult.
        ieclm = [1, 2, 3, 4, 5]
        # ECl Energy
        iecen = [0.1, 0.2, 0.3, 0.6, 1.0]
        # ECl dTime
        iecdt = [0.3, 0.6, 1.0, 1.5, 2.5]
        # ECl.Therm.
        iecth = [5, 10, 15, 20, 25]
        # ECl AnglInt.
        iecan = [2, 4, 8, 20, 40]
        # ECl Size
        iecsz = [0.2, 0.4, 0.6, 0.8, 1.0]
        # ECl EnDens
        ieced = [0.5, 1, 2, 5, 20]
        # ECl En/Mult
        iecem = [0.1, 0.2, 0.3, 0.4, 0.8]
        # ECl En/dTime
        iecet = [0.4, 0.8, 1.5, 3, 6]
        # ECl EnDens/dT
        ieest = [0.2, 0.4, 0.8, 2, 5]
        # EGCl. Energy
        iegen = [0.3, 0.4, 0.6, 1, 1.5]
        # GmCl Energy
        igcen = [0.3, 0.4, 0.6, 1, 1.5]
        # Initialize counters with zeroes
        zeroes = [0] * 6
        # Counters
        counters = {
    "GamEnergy": zeroes.copy(),
    "GamEnergyShw": zeroes.copy(),
    "EleEnergy": zeroes.copy(),
    "EleEnergyShw": zeroes.copy(),
    "MuEnergy": zeroes.copy(),
    "MuEnergyShw": zeroes.copy(),
    "NeutEnergy": zeroes.copy(),
    "NeutEnergyShw": zeroes.copy(),
    "ProtEnergy": zeroes.copy(),
    "ProtEnergyShw": zeroes.copy(),
    "EClMult": zeroes.copy(),
    "EClMultShw": zeroes.copy(),
    "EClEnergy": zeroes.copy(),
    "EClEnergyShw": zeroes.copy(),
    "ECldT": zeroes.copy(),
    "ECldTShw": zeroes.copy(),
    "EClTherm": zeroes.copy(),
    "EClThermShw": zeroes.copy(),
    "EClSize": zeroes.copy(),
    "EClSizeShw": zeroes.copy(),
    "EClEnDens": zeroes.copy(),
    "EClEnDensShw": zeroes.copy(),
    "EClEnMult": zeroes.copy(),
    "EClEnMultShw": zeroes.copy(),
    "EClEndT": zeroes.copy(),
    "EClEndTShw": zeroes.copy(),
    "EClEnSzdT": zeroes.copy(),
    "EClEnSzdTShw": zeroes.copy(),
    "EGClEnergy": zeroes.copy(),
    "EGClEnergyShw": zeroes.copy(),
    "GmClEnergy": zeroes.copy(),
    "GmClEnergyShw": zeroes.copy(),
}

# For example, to access the third value of "GamEnergy":
#value = counters["GamEnergy"][2]
#print(value)  # This will print 0

# Create a list of prefixes
prefixes = [
    "hgmen", "helen", "heclm", "hecen", "hecdt", "hecth", 
    "hecsz", "heced", "hecem", "hecet", "heets", "hegen", 
    "hgcen", "hmuen", "hnten", "hpren", "mhgme", "mhele", 
    "mhecm", "mhece"
]

# Create a dictionary to hold each list of values
data = {}

# Initialize each list with six zeros
for prefix in prefixes:
    data[prefix] = [0.0 for _ in range(6)]

# Now, data holds all the values. For example, to access 'hgmen3' in C++,
# you'd use data["hgmen"][2] in Python (remember Python indices start from 0).
# Create a list of prefixes
prefixes = [
    "mhedt", "mheth", "mhcsz", "mheed", "mheem", "mheet", 
    "mhets", "mhege", "mhgce", "mhmue", "mhnte", "mhpre"
]

# Create a dictionary to hold each list of values
data = {}

# Initialize each list with six zeros
for prefix in prefixes:
    data[prefix] = [0.0 for _ in range(6)]

# Now, data holds all the values. For example, to access 'mhedt3' in C++,
# you'd use data["mhedt"][2] in Python (remember Python indices start from 0).
#Mean Distance (radius)
# Creating different groups of variable prefixes and their associated counts
groups = {
    "recen": 6, "mrece": 6, "tecen": 6, "mtece": 6,
    "rrecm": 6, "rrece": 6, "rredt": 6, "rreth": 6,
    "rrcsz": 6, "rreed": 6, "rrmue": 6
}

# Special cases where the variable doesn't have numeric suffixes
specials = ["mrgam", "mrele", "mrmuon", "mrneut", "mrprot", "mroth"]

# Creating a dictionary to hold the variables
data = {}

# Initializing variables based on the given groups
for key, count in groups.items():
    data[key] = [0.0 for _ in range(count)]

# Initializing special cases
for special in specials:
    data[special] = 0.0

# Now, for example, to access 'recen3' in C++, you'd use data["recen"][2] in Python.
#Radial analysis
# Variable prefixes grouped by their type
float_groups = {
    "t0r": 6, "elmer": 6, "emers": 6, "ecmer": 6, "elmtr": 6, "emtrs": 6, "clmer": 6,
    "clmtr": 6, "mumer": 6, "mmers": 6, "mumtr": 6, "mmtrs": 6, "emdtr": 6, "emdts": 6,
    "emrat": 6, "emrts": 6
}

int_groups = {
    "nelr": 6, "nelrs": 6, "nmur": 6, "nmurs": 6, "neclr": 6, "necrs": 6
}

# Creating a dictionary to hold the variables
data = {}

# Initializing float variables
for key, count in float_groups.items():
    data[key] = [0.0 for _ in range(count)]

# Initializing int variables
for key, count in int_groups.items():
    data[key] = [0 for _ in range(count)]

# Now, for example, to access 't0r3' in C++, you'd use data["t0r"][2] in Python.
# Assuming fout0, fout1, ... are defined file paths or names

with open(fout0, 'a') as file0:  # Particle Summary
    pass  # Replace with any operations you want to perform on this file

with open(fout1, 'a') as file1:  # Clusters Summary
    pass  # Replace with any operations you want to perform on this file

with open(fout2, 'a') as file2:  # Shower Summary
    pass  # Replace with any operations you want to perform on this file

with open(fout3, 'a') as file3:  # Main Summary
    pass  # Replace with any operations you want to perform on this file

with open(fout4, 'a') as file4:  # Radial distributions Summary
    pass  # Replace with any operations you want to perform on this file

with open(fout5, 'a') as file5:  # Mean height estimation
    pass  # Replace with any operations you want to perform on this file
# Assuming you have defined the necessary variables (e.g., iprhd, mxdst, sigrmx, etc.) previously in Python.

if iprhd == 1:
    # Writing to file0 (Particle Summary)
    with open(fout0, 'a') as file0:
        file0.write("# Particle Summary: DetMxD, DetMSg, PCRZhMin, PCRZhMx: \t{}\t{}\t{}\t{}\tHMin, HMax: \t{}\t{}\n".format(mxdst, sigrmx, zhmin, zhmax, hmin, hmax))
        # Uncomment below if you need the height intervals line
        # file0.write("Height intervals:  0, 12, 16, 18, 20, 22, 28, 40, 150\n")
        # file0.write("---------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
        file0.write("#PrCR\tLPCREn\tHghPCR\tIShow\tPId\tRPart\tenert\tAzAPar\tZhAPar\n")
    
    # Writing to file1 (Clusters Summary)
    with open(fout1, 'a') as file1:
        file1.write("# Clusters Summary: DetMxD, DetMSg, PCRZhMin, PCRZhMx: \t{}\t{}\t{}\t{}\tHMin, HMax: \t{}\t{}\n".format(mxdst, sigrmx, zhmin, zhmax, hmin, hmax))
        # Uncomment below if you need the height intervals line
        # file1.write("Height intervals:  0, 12, 16, 18, 20, 22, 28, 40, 150\n")
        # file1.write("---------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
        headers = [
            "#PrCR", "LPCREn", "HghPCR", "IShow", "iClust", 
            "ClsIdC", "sClIdC", "NpartC", "NGamC", "NEleC", "NMuCl", 
            "XmClst", "YmClst", "RmClst", "SigRCl", 
            "TmClst", "dTClst", "sTClst",
            "FstPID", "FstPCX", "FstPCY", "FstPCT", 
            "FstPZh", "FstPAz", "FstPPm",
            "LstPID", "LstPCX", "LstPCY", "LstPCT", 
            "LstPZh", "LstPAz", "LstPPm",
            "EClMlt", "ECltEn", "ECltdT", "ECltTh", 
            "ECltSz", "ECltED", "ECltEM", "ECltET",
            "EClETS", 
            "EGmMlt", "EGmCEn", "EGmCdT", "EGmCTh", 
            "EGmCSz", "EGmCED", "EGmCEM", "EGmCET",
            "EGmETS",
            "GmCMlt", "GmClEn", "GmCldT", "GmClTh", 
            "GmClSz", "GmClED", "GmClEM", "GmClET",
            "GmCETS"
        ]
        file1.write("\t".join(headers) + "\n")
file2.write(f"# Showers Summary: DetMxD, DetMSg, PCRZhMin, PCRZhMx: \t{mxdst}\t{sigrmx}\t{zhmin}\t{zhmax}\t HMin, HMax: \t{hmin}\t{hmax}\n")
# "Height intervals:  0, 12, 16, 18, 20, 22, 28, 40, 150" - Está comentado en el código original, así que no lo incluiremos

line = "#\t" * 31  # Hay 31 tabulaciones antes de las variables
line += f"{0}\t{igmen1}\t{igmen2}\t{igmen3}\t{igmen4}\t{igmen5}\t"
line += f"{0}\t{ielen1}\t{ielen2}\t{ielen3}\t{ielen4}\t{ielen5}\t"
line += f"{0}\t{imuen1}\t{imuen2}\t{imuen3}\t{imuen4}\t{imuen5}\t"
line += f"{ieclm1}\t{ieclm2}\t{ieclm3}\t{ieclm4}\t{ieclm5}\t{ieclm5+1}\t"
line += f"{0}\t{iecen1}\t{iecen2}\t{iecen3}\t{iecen4}\t{iecen5}\t"
line += f"{0}\t{iecdt1}\t{iecdt2}\t{iecdt3}\t{iecdt4}\t{iecdt5}\t"
line += f"{0}\t{iecth1}\t{iecth2}\t{iecth3}\t{iecth4}\t{iecth5}\t"
line += f"{0}\t{iecsz1}\t{iecsz2}\t{iecsz3}\t{iecsz4}\t{iecsz5}\t"
line += f"{0}\t{iegen1}\t{iegen2}\t{iegen3}\t{iegen4}\t{iegen5}\t"
line += f"{0}\t{igcen1}\t{igcen2}\t{igcen3}\t{igcen4}\t{igcen5}\t"
line += f"{0}\t{ir1}\t{ir2}\t{ir3}\t{ir4}\t{ir5}\t"
line += f"{0}\t{ir1}\t{ir2}\t{ir3}\t{ir4}\t{ir5}\t"
line += f"{0}\t{ir1}\t{ir2}\t{ir3}\t{ir4}\t{ir5}\n"

file2.write(line)

header = (
    "#PrimCR\tLPCREn\tHghPCR\tFlinSh\tIShow\tNSecP\tNClst\t"
    "NGamS\tNEleS\tNMuS\tNNeutS\t"
    "NProtS\tNOthrS\t"
    "NElClS\tNMuClS\tNMxClS\tNEGClS\t"
    "NGmClS\tNOtClS\t"
    "ElMEnS\tMuMEnS\t"
    "GmMnR\tElMnR\tMuMnR\tNtMnR\tPtMnR\t"
    "OtMnR\t"
    "ElCMnR\tMuCMnR\tMxCMnR\t"
    "EGCMnR\tGmCMnR\tOtCMnR\t"
    "NGmEn1\tNGmEn2\tNGmEn3\tNGmEn4\t"
    "NGmEn5\tNGmEn6\t"
    "NElEn1\tNElEn2\tNElEn3\tNElEn4\t"
    "NElEn5\tNElEn6\t"
    "NMuEn1\tNMuEn2\tNMuEn3\tNMuEn4\t"
    "NMuEn5\tNMuEn6\t"
    "NEClM1\tNEClM2\tNEClM3\tNEClM4\t"
    "NEClM5\tNEClM6\t"
    "NEClE1\tNEClE2\tNEClE3\tNEClE4\t"
    "NEClE5\tNEClE6\t"
    "NECdT1\tNECdT2\tNECdT3\tNECdT4\t"
    "NECdT5\tNECdT6\t"
    "NECTh1\tNECTh2\tNECTh3\tNECTh4\t"
    "NECTh5\tNECTh6\t"
    "NECSz1\tNECSz2\tNECSz3\tNECSz4\t"
    "NECSz5\tNECSz6\t"
    "NEGEn1\tNEGEn2\tNEGEn3\tNEGEn4\t"
    "NEGEn5\tNEGEn6\t"
    "NGCEn1\tNGCEn2\tNGCEn3\tNGCEn4\t"
    "NGCEn5\tNGCEn6\t"
    "NEleR1\tNEleR2\tNEleR3\tNEleR4\t"
    "NEleR5\tNEleR6\t"
    "NMunR1\tNMunR2\tNMunR3\tNMunR4\t"
    "NMunR5\tNMunR6\t"
    "NEClR1\tNEClR2\tNEClR3\tNEClR4\t"
    "NEClR5\tNEClR6\n"
)

file2.write(header)


