# STIR2PETSIRD converter

The purpose of this repo is to provide a starting point for converting list mode PET data read by STIR to PETSIRD. 


## Background
The [Emission Tomography Standardization Initiative (ETSI)](https://etsinitiative.org/)
is working towards establishing a standard for PET Raw Data, called PETSIRD ("PET ETSI Raw Data").
More information is on https://github.com/ETSInitiative/PETSIRD.

## Status
This is still work in progress, but does seem to work on "cylindrical" geometries, including normalisation (no singles nor dead-time info is passed on).

## Example usage
For the [mMR acquisition of the NEMA phantom on Zenodo](https://zenodo.org/records/1304454):

1. fix the EOL conventions in the Siemens .n.hdr
   ```
   sed -i.bak 's/\r\([^\n]\)/\r\n\1/g' < 20170809_NEMA_UCL.n.hdr > 20170809_NEMA_UCL.n.hdr.fixedEOL

2. create a file `norm.par`
   ```
   Bin Normalisation parameters:=
   type:= From ECAT8
    Bin Normalisation From ecat8 :=
    normalisation filename:= 20170809_NEMA_UCL.n.hdr.fixedEOL
    End Bin Normalisation From ecat8:=
   END:=
   ```
3. run the convertor
   ```
   build/STIR_PETSIRD_convertor  mMR_NEMA_norm.petsird 20170809_NEMA_60min_UCL.l.hdr norm.par
   ```
