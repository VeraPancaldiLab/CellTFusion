# Rebuilds the tutorial data that depend on CellTFusion's defaults (run from the package root).
# raw.counts.tuto, counts.norm.tuto, deconv.tuto and deconv_subgroups.tuto are not rebuilt here:
# - raw.counts.tuto / traitdata.tuto: 192 samples of Mariathasan et al. (2018) (IMvigor210)
# - counts.norm.tuto: ADImpute::NormalizeTPM(raw.counts.tuto, log = TRUE)
# - deconv.tuto: multideconv::compute.deconvolution() on raw.counts.tuto (Quantiseq, CBSX, Epidish, DeconRNASeq, DWLS)
# - deconv_subgroups.tuto: multideconv::compute.deconvolution.analysis(deconv.tuto, corr = 0.7) (multideconv >= 0.2.0)

devtools::load_all()

# Clinical data: overall survival (os, in months, and censOS, 1 = death) added from the same study
# `clinical_file` is the clinical table of Mariathasan et al. (2018) with the samples as row names
clinical_file <- "~/Documents/CellTFusion_paper/input/Mariathasan/Input/dataInput_Mariathasan.RData"
clinical_env <- new.env()
load(clinical_file, envir = clinical_env)
data("traitdata.tuto")
traitdata.tuto <- traitdata.tuto[, c("Best.Confirmed.Overall.Response", "binaryResponse", "Enrollment.IC")]
traitdata.tuto$os <- clinical_env$traitData[rownames(traitdata.tuto), "os"]
traitdata.tuto$censOS <- clinical_env$traitData[rownames(traitdata.tuto), "censOS"]
usethis::use_data(traitdata.tuto, overwrite = TRUE)

# TF activity: CollecTRI regulons, ulm, scaled (defaults of compute.TFs.activity())
# `universe` is the CollecTRI collection (decoupleR::get_collectri(organism = "human", split_complexes = FALSE))
data("counts.norm.tuto")
universe <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)
tfs.tuto <- compute.TFs.activity(counts.norm.tuto, universe = universe, return = FALSE)
usethis::use_data(tfs.tuto, overwrite = TRUE)

# TF network: defaults of compute.WTCNA()
network.tuto <- compute.WTCNA(tfs.tuto, return = FALSE)
usethis::use_data(network.tuto, overwrite = TRUE)
