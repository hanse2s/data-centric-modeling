# DataCentricSDM
This is a repository for the manuscript: Hansen, S. E., M. J. Monfils, R. A. Hackett, R. T. Goebel, and A. K. Monfils. 2024. Data-centric species distribution modeling: Impacts of modeler decisions in a case study of invasive European frog-bit. Applications in Plant Sciences 12: e11573. https://doi.org/10.1002/aps3.11573.

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

## If you plan to use or adapt this code, please cite the Open Science Framework project: https://osf.io/y6meq/

The scripts on the OSF project and this repository are nearly identical. The one change is that the OSF project uses a series of setwd() commands for accessibility to new R users, while this is structured like a project. Use whatever format you prefer! No matter which one you use, please cite the OSF project.

I am always open to comments and revisions. Get in touch!

## Repository structure
* code
  * SDM1_DataDownload
  * SMD2_ExplanatoryDataImportAndProcessing
  * SDM3_ReponseDataDownloadImportProcessing
  * SDM4_GetBackgroundPoints_FinalProcessingStep
  * SDM5_FitAndEvaluateModels (start here if you just want to run the model code)
  * SDM6_AnalysisAndVisualization (start here if you just want to recreate paper tables and figures)
  * SDM7_DataExplorationForManuscript
  * SDM8_AlternativeQuestion3
  * SDM9_DelineateLargeScale
 
* data
  * If you choose to use the Michigan European frog-bit dataset for anything else, cite the [GBIF version of it](https://www.gbif.org/dataset/71454d8a-6e9c-49f5-bf37-353f9ad2e2b9), not this one.

 * figures

 * output
   * Many of the raw data files are large and not mine to share, so I'm not including them in the actual output folder. I am including all of the final data used in the models, so you can start with the models in Step 5 without having to run earlier steps.

  * tables
