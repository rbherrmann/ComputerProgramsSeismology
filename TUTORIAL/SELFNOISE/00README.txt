Estimation of instrumental noise through a corss-correlation of
outputs of two instruments recording side-by-side

     Holcomb, G. L. (1989).
     A Direct Method for Calculating Instrument Noise Levels
     in Side-by-Side Seismometer Evaluations,
     USGS Open-File Report 89-214


Data processing steps

1. cd FPSD
   This directory has the code required for the coherency analysis
   and also codes for PSD noise estimateion

   This codes here are an implementation of the USGS Open File Report on
   noise measurements

   There are several subdirectories

   src - source code and Makefile
   bin - location of executables
   TESTPSD - location of a test data set

   

   a) cd src
      edit the Makefile to give the name of your compilers.
      The Makefile is set up for gfortran and gcc

      I also include the plot program xyplt2 which  requires the
      graphics libraries of Computer Programs in FORTRAN
      You do not require this since your can plot the output 
      using GMT
   
      The psd and coher programs are pure FORTRAN programs

      make all

      This will create two programs: psd for noise analysis and
      coher for the comparison of output of two sensors

    b) cd TESTPSD

       This directory contains the codes for performing a 
       noise analysis
    
       A data set is provided for the GSN station ALQ

       The files ending with .cmp are the data files that I was
       given

       The files ending with .asac are Sac files in ASCII format. 
       You can use the Computer Programs in Seismology program
       asctosac to covert from this format to the SAC binary format
       for your computer, e.g.,

       asctosac anmbz_08.asac anmbz_08.sac

       The scripts DOITALQBIN and DOITALQASC are essentially the same
       except that the first uses the binary sac file and the other
       uses the corresponding ascii file.

       The program psd outputs the noise spectrum in acceleration

       The program is controlled from the command line:

       psd [-h] [-D -V -A] -R response_file sac_files

       where response_file is an ASCII file with the first column
       indicating frequency and the second indicating the response.
       The second column is in units of counts/meter for the -D flag,
       counts/m/sec for the -V flag and coutns/m/s/s for the -A flag


       The output of this program is redicted to a file, here called
       junk.alq which consists of two columns. The first column is
       period and the second is thepower spectral density of the noise
       in decibels or mathematically

       PSD 10*log10(m^2/s^4 .Hz)

       The program xyplt2 creates a CALPLOT graphics file plotXXXXX,
       where XXXXX is a unique number. The PSD from the program 
       is plotted and is compared to the low and high noise models,
       nlnm.acc and nhnm.acc respectively

       The CALPLOT graphcis file is renamed PLOTALQ and converted to
       an Encapsulated PostScript file (EPS) using the
       Computer Program in Seismology routine plotnps
       


   
    
2. cd TEST.ANMO
       This directory has the scripts for performing the coherency analysis for the
       two LHZ sensors at ANMO


	* CHANNEL(NSCL)IUANMO LHZ00
	* NETWORK      IU
	* STATION      ANMO
	* COMPONENT    LHZ
	* LOCATION     00
	* RATE (HZ)    1.0
	* INSTRMNTTYPE Geotech KS-54000 Borehole Seismometer

	* CHANNEL(NSCL)IUANMO LHZ10
	* NETWORK      IU
	* STATION      ANMO
	* COMPONENT    LHZ
	* LOCATION     10
	* RATE (HZ)    1.0
	* INSTRMNTTYPE Guralp CMG3-T Seismometer (borehole)

	This example starts from a SEED volume obtained from IRIS

	The shell script DOALL 
	a) unpacks the SEED volume to get the waveforms in a SAC format and
           the response in a RESPONSE file
	b) uses the resposne file to create an amplitude and phase table
	   for the acceleration sensitivity

        YOU WILL NEED rdseed, evalresp, gsac 

     I got the pdf_E2010.220_S2010.220_cBHZ_l10_nIU_sANMO.png indepdendent estimate form
     IRIS for this day so that I can check my code,  The KS borehole is much 
     quieter than the Guralp

3. cd KMI.TEST
    This uses KMI data from January 2010 and compares the STS-1 to
    STS-2 sensors. The STS-2 high gain is noisier  but more sensitive





