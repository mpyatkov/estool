# EStool

This script identifies enriched peaks in a foreground set of peak regions (set of genomic regions, in .bed file format, with bed files placed in the foreground directory) when compared to a set of background peaks (.bed files placed in the background directory). Overlap analysis is calculated for each peak set against a database comprised of known sets of peak regions, which are organized into a directory of bioregions (bio_regions directory, with each subdirectory containing sets of related genomic regions of interest). Enrichment is quantified using a Fisher's exact test, which assesses the statistical significance of the overlap between the foreground peaks and database peaks relative to the background peaks. 

---
### Acknowledgements

Created by: Maxim Pyatkov, laboratory of David Waxman, Dept. of Biology, Boston University

Grant support: Development of this script was carried out in the laboratory of David J Waxman, Boston University, supported in part by NIH grants R01-DK121998 (Growth Hormone Regulation of Sex Differences in Liver Metabolism) and R01-ES024421 (Epigenetic Actions of Environmental Chemicals) to DJW.
