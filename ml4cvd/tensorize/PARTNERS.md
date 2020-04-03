# Partners ECG
Organizing and Tensorizing MUSE 12-lead ECGs

## Table of Contents
1. [Organizing XMLs and Removing Duplicates](#organizing-xmls-and-removing-duplicates)
2. [Tensorizing XMLs to HDF5](#tensorizing-xmls-to-hdf5)
3. [ECG Data Structure](#ecg-data-structure)
4. [Extracting ECG Metadata](#extracting-ecg-metadata)
5. [Other documentation](#other-documentation)

## Organizing XMLs and Removing Duplicates
`ingest/partners_ecg/organize_xml.py` moves XML files from a single directory into the appropriate yyyy-mm directory.

`ingest/partners_ecg/remove_xml_duplicates.py` finds and removes exact duplicate XML files, as defined by every bit of two files being identical, determined via SHA-256 hashing.  

## Tensorizing XMLs to HDF5
`tensorize_partners` mode in `recipes.py` extracts data from all XML files and saves as [HDF5 files](https://www.hdfgroup.org). Tensorization also removes duplicates that contain nearly the same information, except for minor differences, for example minor version changes in acquisition software. This duplicate detection is done by matching patient-date-time fields.  

This mode is called with the following arguments:  
`--xml_folder` to specify the directory containing ECG XMLs.  
`--tensors` to specify the directory where tensorized HD5 files should be saved.  

All the ECGs belonging to one patient, identified by medical record number (MRN), will be saved to one HD5, indexed by ECG acquisition date and time:  
```
<MRN>.hd5
└--partners_ecg_rest
   |
   |--date_1
   |  └--ECG Data
   |
   └--date_2
      └--ECG Data
```

## ECG Data Structure
Voltage is saved from XMLs as a dictionary of numpy arrays indexed by leads in the set `("I", "II", "V1", "V2", "V3", "V4", "V5", "V6")`, e.g.:

```
voltage = {'I': array([0, -4, -2, ..., 7]),
          {'II': array([2, -9, 0, ..., 5]),
          ...
          {'V6': array([1, -4, -3, ..., 4]),
```

Every other element extracted from the XML is returned as a string, even if the underlying primitive type is a number (e.g. age). Here are some of the more important elements:

```
acquisitiondate
atrialrate
dateofbirth
diagnosis_computer
diagnosis_md
ecgsamplebase
ecgsampleexponent
gender
heightin
location
locationname
overreaderfirstname
overreaderid
overreaderlastname
patientid
paxis
poffset
ponset
printerval
qoffset
qonset
qrscount
qrsduration
qtcfrederica
qtcorrected
qtinterval
race
raxis
taxis
toffset
ventricularrate
weightlbs
```

## Extracting ECG metadata

`explore` mode in `recipes.py` extracts data specified by `--input_tensors` from all HD5 files given to `--tensors` and calculates summary statistics. Additionally, all metadata is saved to a large CSV file:  

This CSV file will be used to construct a performant, queryable database to identify future cohorts for research projects.

## Other documentation
GE documentation is stored in a shared Partners Dropbox folder ([link](https://www.dropbox.com/sh/c5tgm0lory72ge0/AADqKvUicDdyWzHYhtad0lU4a?dl=0)), including 1. physician's guide to the Marquette 12SL ECG analysis program, 2. guide to MuseDB search, and 3. Muse v9 XML developer's guide.
