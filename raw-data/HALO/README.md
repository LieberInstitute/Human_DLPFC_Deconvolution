The HALO files used in the final publication are those from [`Annotation_refinement_20230922`](Annotation_refinement_20230922/). They are used by [`code/03_HALO/01_import_HALO_data.R`](https://github.com/LieberInstitute/Human_DLPFC_Deconvolution/blob/66491598482f157df86f3c9487863179177d4e6b/code/03_HALO/01_import_HALO_data.R). Other files were from earlier analysis iterations.


### Metadata + Notes
`Deconvolution_HALO_Analysis.xlsx`  
Downloaded from [this google sheet](https://docs.google.com/spreadsheets/d/1bld6g-7MN2G18b8hwXEkC0dDOhWlzmFTAmpcPYdgXjI) on 2025-05-04. 
Other info is in [this google doc](https://docs.google.com/document/d/1GXbV134CdPmMcSw9pPeSKeHQ8tBMCf0gLuf5Sq4RF4g)



### Re-exported STAR R5 data
Uploaded my Kristen Maynard on 8/3/32022 to fix 'SLC17A' vs. 'SLC17A7' typo effecting star samples from R5

Copied to working dir as follows:
```
cp Re-export/Exported_Typo_Correction_Files/8325A_star_redo/ObjectData/HALO_B1_R5_Star_Br8325A_Fused.tif_object_Data.csv Deconvolution_HALO_analysis/STAR/R5_HA_Re-export_Analysis/HA_R5_8325A_Star.csv
cp Re-export/Exported_Typo_Correction_Files/8325M_star_redo/ObjectData/HALO_B1_R5_Star_Br8325M_Fused.tif_object_Data.csv Deconvolution_HALO_analysis/STAR/R5_HA_Re-export_Analysis/HA_R5_8325M_Star.csv
cp Re-export/Exported_Typo_Correction_Files/8667A_star_redo/ObjectData/HALO_B1_R5_Star_Br8667A_Fused.tif_object_Data.csv Deconvolution_HALO_analysis/STAR/R5_HA_Re-export_Analysis/HA_R5_8325M_Star.csv

```
