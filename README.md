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
The script [downloadThredds.py](./downloadThredds.py) was used to locally download monthly hindcast files from Wave Information Studies node 63218. See [getDataFRF.py](https://github.com/erdc/getdatatestbed) for functions to automatically download FRF buoy observations from the server rather than the hindcast. 
The above figure showing all data used in this study can be generated with [F1studyMap.py](./F1studyMap.py) using either the provided 'sandbarsSouthernTransect' file or by running [AlloadExtractData.py](./A1loadExtractData.py) to generate a new binary file from a directory of FRF data files (for example, with newer surveys since July 2022).


### Complex EOFs

Complex EOFs are applied to the dataset of cross-shore profiles using [B1complexEOFs.py](./B1complexEOFs.py) to produce the following two dominant modes:


![image](https://github.com/anderdyl/frfSandbars/blob/master/figure2.png)

Historically observed coherent sandbar behavior captured in a) mode 1 representing offshore migrations and b) mode 2 representing onshore migration. Average migrations with respect to temporal amplitude and propagation speed for c) offshore migrations and d) onshore migrations provided as a function of the degrees of mode progression and the average time elapsed.
Constructive interference of the two modes results in a profile reconstruction, where examples in e) and f) are the summation of the phases denoted by the corresponding yellow and green lines in c) and d).

### Phase Space

![image](https://github.com/anderdyl/frfSandbars/blob/master/figure3.png)



![image](https://github.com/anderdyl/frfSandbars/blob/master/conceptualFlowchart.png)

### Files
- [A1loadExtractData.py](./A1loadExtractData.py) Loading, averaging, interpolating transects from a directory of netcdf files containing FRF survey transects.
- [B1complexEOFs.py](./B1complexEOFs.py) Applying CEOFs to the geomorphology data set.
- [C1phaseSpace.py](./C1phaseSpace.py) Identifying profile characteristics for all interferences of the dominant two modes of sandbar migration.
- [F1studyMap.py](./F1studyMap.py) Creating the top figure in this README.
- [F2ceofOnOff.py](./F2ceofOnOff.py) Creating above figure with two different dominant mode of sandbar migration.
- [F3phaseDiagram.py](./F3phaseDiagram.py) Creating annotated phase space diagram provided above.
- [downloadThredds.py](./downloadThredds.py) Local download of CHL Thredds data, with script set for WIS waves.