Image processing and analyzing code for chromatin structures in electron microscopy images.

Files:
 
denoiseTiff.py -removes unspecific noise from the EM tiff stack
denoiseSettings.py- configuration file for denoiseTiff.py
AnalyzeStructure.py- analises EMISH signal- alignes structure in principle axis coordinates, calculates mass centers and
segments domains
AnalysisParameters.py- configuration file for AnalyzeStructure.py
get_gradient (.cpp and .pyd)- auxialiry files 
SegmentDomain.py, Misc.py-auxialiry files

Requirements: Python2.7 with data science packages (tested using Anaconda with Python 2.7.14 under Windows10)
If the working directory is not in python library path, copy the file get_gradient.pyd into the library path.
For other systems (e.g. Linux) get_gradient.cpp has to be compiled. Execution time denoiseSettings.py~ 30 sec. AnalyzeStructure.py ~ 10 mins. for a single example 

Usage:
1. Set correct paths (including input file name /three dimensional tiff/ and output folder) and options in denoiseSettings.py
2. Run python denoiseTiff.py to remove unspecific background
3. Set correct paths (including input file name with segmented structure /obtained from the previous step/ and output folder) and options in AnalysisParameters.py
4. Run python AnalyzeStructure.py to center and align it in princliple axis coordinate system, 
to determine local density centers and to segment it into domains.

The description of adjustable parameters given in files AnalysisParameters.py and denoiseSettings.py

Developed by Blazej Ruszczycki, 2019, e-mail:b.ruszczycki@nencki.gov.pl
