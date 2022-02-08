Longitudinal examination of perfusion and angiogenesis markers in primary colorectal tumors shows distinct signatures for metronomic and maximum-tolerated dose strategies
=========================================================================================================================================================================

[![License: CC BY 4.0](https://img.shields.io/badge/License%20All-CC%20BY%204.0-lightgrey)](https://creativecommons.org/licenses/by/4.0/) 
[![License: CC0-1.0](https://img.shields.io/badge/License%20Parts-CC0%201.0-lightgrey)](http://creativecommons.org/publicdomain/zero/1.0/)

This repository contains code, raw data, scripts and manuscript files for a study of tumor response to treatment (metronomic and maximum-tolerated dose chemotherapy) in a primary murine model of colorectal cancer.

Repository structure:

- code: Contains scripts that create figures, run statistical analyses (in R)
- data: Files called by scripts in `code`. Contain raw qPCR, DRS, and imaging data.
- manuscript: Contains all necessary files (LaTeX, bibliography, etc. )to create the manuscript. The manuscript is created by knitting `Complete_manuscript.Rmd`
  - The directory `supplementary _materials` contains files used to create the Appendix, which can be compiled to PDF by knitting `SM-doc.Rmd`


## How to use this repository

The best way to use this work is to fork the repo. From there the following approches can be used:

- The manuscript can be compiled.
- Individual scripts can be run, but they require the libraries found in the setup chunk in `Complete_manuscript.Rmd` to work.
- Imaging data can be found in [BioStudies](https://www.ebi.ac.uk/biostudies/studies/S-BIAD323) under accession number S-BIAD323.

# License

This entire repository is licensed under a CC BY 4.0 License, which allows reuse with attribution. However, certain files are released under the CC0 1.0 public domain dedication. The files indicated below are dual licensed under CC BY 4.0 and CC0 1.0:

- `02-Materials-and-Methods.Rmd`
- `03-Results.Rmd`
- `DRS_GAMs.R`
- `Imaging_data_analysis.R`
- `qPCR_GAMs_diagnostics.R`
- `qPCR_GAMs_pwcomps.R`
- `qPCR_GAMs.R`

All other files are under CC BY 4.0.

# Contact

If you have questions or comments please contact Ariel Mundo (aimundo AT uark.edu)
