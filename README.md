# DICOMConvertServer
Automatic service to convert DICOM image series into matlab file format

## Introduction
The DICOM Convert Server (DCS) is a Matlab tool that can be used to convert DICOM files to *.mat format for faster and easier file access with FlowQuant. The DCS is intended to run as a background service on a host computer, typically on a network. DCS has been designed to work with 3D and 4D image volumes from a variety of modalities and vendors. Nevertheless, it may need to adapated to new data formats.

The DCS monitors and incoming directory for new files and converts them in to *.mat files that are saved to an outgoing directory. A temporary intermediate directory is also used as part of the conversion. These directories, along with other settings, are specified in a DICOMServer.dat file that in ASCII text format.

DCS can be paired with a DICOM receiving node to automatically convert files as they are received.

## Configuring the DICOM Convert Server
All configurations of the DICOM Convert Server are through the configuration file (DICOMServer.dat) which is stored in the same directory as the server executable (DICOMConvertServer.exe). Entries follow the format of entryName = entryValue. Unknown entry names are ignored and default values (see table above) will be used. All entry names are lower case and case sensitive. Comment lines are preceded by %.

If no configuration file exists the following default values will be used:

### Default Configuration Settings
Field name - Meaning {Default value}
* indir - Directory in which DICOM files to convert may be found.	{C:\DICOMImport}
* outdir - Directory in which to store the converted files. {C:\DataStash}
* logdir - Directory in which to store log files of the conversion process.	{C:\DataStash\Log}
* tempdir - Directory in which temporary intermediate files are stored. {C:\DataStash\temp}
* cleanup - Option for cleaning up the input directory as part of conversion. 1-enable, 0-disable. It is strongly recommended that this field be left as 1. {1}
* overwrite -  Option to overwrite output files with the same name (1), or give conflicting file names sequential numbers (0). {0}
