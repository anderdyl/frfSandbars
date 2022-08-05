# frfSandbars

Place-holder for zenodo doi: 
[![DOI](https://zenodo.org/badge/381731651.svg)](https://zenodo.org/badge/latestdoi/381731651)

### Overview

Codes in this repo process topo-bathy surveys from the Field Research Facility to assess sandbar dynamics observed since 1980.
The intent is to develop a cross-shore perspective on morphology evolution that considers the complete profile evolution (as opposed to 1-dimensional variables such as sandbar cross-shore loction or sandbar depth).
Complex Empirical Orthogonal Functions (CEOFs) are used to isolate offshore propagrations (dominant mode) and onshore migrations (2nd mode).




![image](https://github.com/anderdyl/frfSandbars/blob/master/figure1.png)

a) All elevations in a 200-m stretch along the southern portion of Duck were averaged to derive 635 different transects shown in b) the cross-shore with mean (solid) and standard deviation (dotted) and c) bathymetric anomalies with respect to time. d) The non-dimensional fall velocity offshore of Duck with hourly values (gray) and month-long rolling averages (black).

### Data

Surveys are freely available on the Coastal & Hydraulics Laboratory [THREDDS Server](https://chldata.erdc.dren.mil/), specially under /frf_data/geomorphology/elevationTransects/survey/data/. 
More than 1000 surveys exist in this directory, including daily surveys associated with focused experiments at the FRF. The daily frequency surveys have been removed from this analysis, focusing on the fort-nightly to monthly intervals.
All data within a 200 meter alongshore segment of beach to the south of the pier was averaged to a single cross-shore profile to reduce the effect of 3-dimensional morphology.
All transects where subsequently aligned such that the origin of the profile began at the mean high high water (MHHW), effectively isolating to process occurring subaqueously and creating a proces-oriented reference frame.
Loading, extracting transects in the area of interest, and averaging in the alongshore with accomplished with [A1loadExtractData.py](./A1loadExtractData.py).

Wave data is also freely available at the Coastal & Hydraulics Laboratory [THREDDS Server](https://chldata.erdc.dren.mil/).

### Complex EOFs

![image](https://github.com/anderdyl/frfSandbars/blob/master/figure2.png)

Historically observed coherent sandbar behavior captured in a) mode 1 representing offshore migrations and b) mode 2 representing onshore migration. Average migrations with respect to temporal amplitude and propagation speed for c) offshore migrations and d) onshore migrations provided as a function of the degrees of mode progression and the average time elapsed.
Constructive interference of the two modes results in a profile reconstruction, where examples in e) and f) are the summation of the phases denoted by the corresponding yellow and green lines in c) and d).


![image](https://github.com/anderdyl/frfSandbars/blob/master/figure3.png)
