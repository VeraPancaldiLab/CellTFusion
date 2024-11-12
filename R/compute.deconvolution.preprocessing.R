#' Compute deconvolution preprocessing
#'
#' Give consistent names and patterns following the method_signature_cell structure to the deconvolution features
#'
#' @param deconv A dataframe with the unprocessed deconvolution features
#'
#' @return A matrix of the preprocessed deconvolution features with fixed and consistent names across the different methods and signatures
#'
#' @export
#'
#' @examples
#'
#' deconvolution = compute.deconvolution.preprocessing(raw_deconvolution)
#'
#'
compute.deconvolution.preprocessing = function(deconv){
  cat("Preprocessing deconvolution features...............................................................\n\n")

  #Convert mcp and xcell features to proportions by row-scaling
  for (i in 1:nrow(deconv)) {
    deconv[,grep("MCP", colnames(deconv))][i,] = deconv[,grep("MCP", colnames(deconv))][i,]/sum(deconv[,grep("MCP", colnames(deconv))][i,])
    deconv[,grep("XCell", colnames(deconv))][i,] = deconv[,grep("XCell", colnames(deconv))][i,]/sum(deconv[,grep("XCell", colnames(deconv))][i,])
  }

  ##### Edit cell names for consistency across features
  ##### Macrophages (M0, M1, M2)
  Macrophages = deconv[,grep("acrophage", colnames(deconv))]
  M0 = deconv[,grep("M0", colnames(deconv))]
  M1 = deconv[,grep("M1", colnames(deconv))]
  M2 <- deconv[,grep("M2", colnames(deconv))]
  if(length(grep("LM22", colnames(M2)))>0){M2 <- M2[,-grep("LM22", colnames(M2))]}else{M2 <- M2}
  test = deconv[,grep("LM22", colnames(deconv))]
  test = test[,grep("Macrophages.M2", colnames(test))]
  M2 = cbind(M2, test)
  Macrophages = Macrophages[,-which(colnames(Macrophages)%in%c(colnames(M0), colnames(M1), colnames(M2)))]
  deconv = deconv[,-which(colnames(deconv)%in%c(colnames(Macrophages), colnames(M0), colnames(M1), colnames(M2)))]

  colnames(Macrophages) = stringr::str_replace(colnames(Macrophages), "Macrophages", "Macrophages.cells")
  colnames(Macrophages) = stringr::str_replace(colnames(Macrophages), "Macrophage(?!.)", "Macrophages.cells")
  colnames(M0) = stringr::str_replace(colnames(M0), "Macrophages_M0", "Macrophages.M0")
  colnames(M0) = stringr::str_replace(colnames(M0), "_M0", "_Macrophages.M0")
  colnames(M1) = stringr::str_replace(colnames(M1), "Macrophages_M1", "Macrophages.M1")
  colnames(M1) = stringr::str_replace(colnames(M1), "Macrophages_M1", "Macrophages.M1")
  colnames(M1) = stringr::str_replace(colnames(M1), "Macrophage_M1", "Macrophages.M1")
  colnames(M1) = stringr::str_replace(colnames(M1), "_M1", "_Macrophages.M1")
  colnames(M2) = stringr::str_replace(colnames(M2), "Macrophage_M2", "Macrophages.M2")
  colnames(M2) = stringr::str_replace(colnames(M2), "Macrophages_M2", "Macrophages.M2")
  colnames(M2) = stringr::str_replace(colnames(M2), "_M2", "_Macrophages.M2")

  ##### Monocytes
  Monocytes = deconv[,grep("Mono|mono", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(Monocytes))]
  colnames(Monocytes) = stringr::str_replace(colnames(Monocytes), "Monocytic_lineage", "Monocytes")
  colnames(Monocytes) = stringr::str_replace(colnames(Monocytes), "Monocyte(?!s)", "Monocytes")
  colnames(Monocytes) = stringr::str_replace(colnames(Monocytes), "Mono(?!cytes)", "Monocytes")
  colnames(Monocytes) = stringr::str_replace(colnames(Monocytes), "Mono(?!cytes)", "Monocytes")

  ##### Neutrophils
  Neutrophils <- deconv[,grep("Neu", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(Neutrophils))]

  colnames(Neutrophils) = stringr::str_replace(colnames(Neutrophils), "Neutrophil(?!s)", "Neutrophils")
  colnames(Neutrophils) = stringr::str_replace(colnames(Neutrophils), "Neu(?!trophils)", "Neutrophils")

  ###NK cells
  NK = deconv[,grep("NK", colnames(deconv))]
  NKT = NK[,grep("NKT", colnames(NK))]
  NK.activated <- NK[,grep("activated", colnames(NK), value = TRUE)]
  NK.resting <- NK[,grep("resting", colnames(NK), value = TRUE)]
  NK = NK[,-which(colnames(NK)%in%c(colnames(NK.activated), colnames(NK.resting), colnames(NKT)))]
  deconv = deconv[,-which(colnames(deconv)%in%c(colnames(NK), colnames(NK.activated), colnames(NK.resting), colnames(NKT)))]

  colnames(NK) = stringr::str_replace(colnames(NK), "NK(?!.)", "NK.cells")
  colnames(NK) = stringr::str_replace(colnames(NK), "NK_cells", "NK.cells")
  colnames(NK) = stringr::str_replace(colnames(NK), "NK_cell", "NK.cells")
  colnames(NKT) = stringr::str_replace(colnames(NKT), "NKT_", "NKT.")
  colnames(NK.activated) = stringr::str_replace(colnames(NK.activated), "NK.cells.activated", "NK.activated")
  colnames(NK.activated) = stringr::str_replace(colnames(NK.activated), "NK.cells_activated", "NK.activated")
  colnames(NK.resting) = stringr::str_replace(colnames(NK.resting), "NK.cells.resting", "NK.resting")
  colnames(NK.resting) = stringr::str_replace(colnames(NK.resting), "NK.cells_resting", "NK.resting")

  ###CD4 cells
  CD4 <- deconv[,grep("CD4", colnames(deconv))]
  CD4.memory.activated = CD4[,grep("activated", colnames(CD4))]
  CD4.memory.resting = CD4[,grep("resting", colnames(CD4))]
  CD4.naive = CD4[,grep("naive", colnames(CD4))]
  CD4.non.regulatory = CD4[,grep("regulatory", colnames(CD4))]
  CD4 = CD4[,-which(colnames(CD4)%in%c(colnames(CD4.memory.activated), colnames(CD4.memory.resting), colnames(CD4.naive), colnames(CD4.non.regulatory)))]
  deconv = deconv[,-which(colnames(deconv)%in%c(colnames(CD4), colnames(CD4.memory.activated), colnames(CD4.memory.resting), colnames(CD4.naive), colnames(CD4.non.regulatory)))]

  colnames(CD4) = stringr::str_replace(colnames(CD4), "CD4", "CD4.cells")
  colnames(CD4) = stringr::str_replace(colnames(CD4), "T.cells.CD4.cells", "CD4.cells")
  colnames(CD4.memory.activated) = stringr::str_replace(colnames(CD4.memory.activated), "CD4_memory_activated", "CD4.memory.activated")
  colnames(CD4.memory.resting) = stringr::str_replace(colnames(CD4.memory.resting), "CD4_memory_resting", "CD4.memory.resting")
  colnames(CD4.naive) = stringr::str_replace(colnames(CD4.naive), "CD4_naive", "CD4.naive")
  colnames(CD4.naive) = stringr::str_replace(colnames(CD4.naive), "CD4._naive", "CD4.naive")
  colnames(CD4.naive) = stringr::str_replace(colnames(CD4.naive), "T.cells.CD4.naive", "CD4.naive")
  colnames(CD4.naive) = stringr::str_replace(colnames(CD4.naive), "T_cells_CD4.naive", "CD4.naive")
  colnames(CD4.non.regulatory) = stringr::str_replace(colnames(CD4.non.regulatory), "T_cell_CD4._.non.regulatory.", "T.cells.non.regulatory")

  ####CD8
  CD8 <- deconv[,grep("CD8", colnames(deconv))]

  deconv = deconv[,-which(colnames(deconv)%in%colnames(CD8))]

  colnames(CD8) = stringr::str_replace(colnames(CD8), "T_cells_CD8", "CD8.cells")
  colnames(CD8) = stringr::str_replace(colnames(CD8), "T_cell_CD8", "CD8.cells")
  colnames(CD8) = stringr::str_replace(colnames(CD8), "CD8_T_cells", "CD8.cells")
  colnames(CD8) = stringr::str_replace(colnames(CD8), "T.cells.CD8", "CD8.cells")
  colnames(CD8) = stringr::str_replace(colnames(CD8), "CD8(?!.)", "CD8.cells")
  colnames(CD8) = stringr::str_replace(colnames(CD8), "CD8.cells.", "CD8.cells")

  ##### Regulatory T cells
  Tregs = deconv[,grep("regs", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(Tregs))]

  colnames(Tregs) = stringr::str_replace(colnames(Tregs), "T_cell_regulatory_.Tregs.", "T.cells.regulatory")
  colnames(Tregs) = stringr::str_replace(colnames(Tregs), "T.cells.regulatory..Tregs.", "T.cells.regulatory")
  colnames(Tregs) = stringr::str_replace(colnames(Tregs), "Tregs", "T.cells.regulatory")

  ##### Helper T cells
  Thelper = deconv[,grep("helper", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(Thelper))]

  colnames(Thelper) = stringr::str_replace(colnames(Thelper), "T.cells.follicular.helper", "T.cells.helper")
  colnames(Thelper) = stringr::str_replace(colnames(Thelper), "T_cells_follicular_helper", "T.cells.helper")

  ##### Gamma delta T cells
  Tgamma = deconv[,grep("gamma", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(Tgamma))]

  colnames(Tgamma) = stringr::str_replace(colnames(Tgamma), "T_cells_gamma_delta", "T.cells.gamma.delta")
  colnames(Tgamma) = stringr::str_replace(colnames(Tgamma), "T_cell_gamma_delta", "T.cells.gamma.delta")

  ##### Dendritic cells (activated, resting)
  Dendritic = deconv[,grep("endritic", colnames(deconv))]
  Dendritic.activated = Dendritic[,grep("activated", colnames(Dendritic))]
  Dendritic.resting = Dendritic[,grep("resting", colnames(Dendritic))]
  Dendritic = Dendritic[,-which(colnames(Dendritic)%in%c(colnames(Dendritic.activated), colnames(Dendritic.resting)))]
  deconv = deconv[,-which(colnames(deconv)%in%c(colnames(Dendritic), colnames(Dendritic.activated), colnames(Dendritic.resting)))]

  colnames(Dendritic) = stringr::str_replace(colnames(Dendritic), "Myeloid_dendritic_cells", "Dendritic.cells")
  colnames(Dendritic) = stringr::str_replace(colnames(Dendritic), "Myeloid_dendritic_cell", "Dendritic.cells")
  colnames(Dendritic) = stringr::str_replace(colnames(Dendritic), "Dendritic_cells", "Dendritic.cells")

  colnames(Dendritic.activated) = stringr::str_replace(colnames(Dendritic.activated), "dendritic_cell_activated", "Dendritic.activated.cells")
  colnames(Dendritic.activated) = stringr::str_replace(colnames(Dendritic.activated), "Dendritic.cells.activated", "Dendritic.activated.cells")
  colnames(Dendritic.activated) = stringr::str_replace(colnames(Dendritic.activated), "Dendritic_cells_activated", "Dendritic.activated.cells")
  colnames(Dendritic.resting) = stringr::str_replace(colnames(Dendritic.resting), "Dendritic.cells.resting", "Dendritic.resting.cells")
  colnames(Dendritic.resting) = stringr::str_replace(colnames(Dendritic.resting), "Dendritic_cells_resting", "Dendritic.resting.cells")

  ##### CAF cells
  CAF = deconv[,grep("CAF|Cancer_associated_fibroblast", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(CAF))]

  colnames(CAF) = stringr::str_replace(colnames(CAF), "Cancer_associated_fibroblast", "CAF")
  colnames(CAF) = stringr::str_replace(colnames(CAF), "CAFs", "CAF")

  ##### Cancer cells
  Cancer = deconv[,grep("ancer", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(Cancer))]

  colnames(Cancer) = stringr::str_replace(colnames(Cancer), "cancer", "Cancer")
  colnames(Cancer) = stringr::str_replace(colnames(Cancer), "Cancer.cells", "Cancer")

  malignant = deconv[,grep("alignant", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(malignant))]

  colnames(malignant) = stringr::str_replace(colnames(malignant), "Malignant", "Cancer")
  colnames(malignant) = stringr::str_replace(colnames(malignant), "Cancer_cells", "Cancer")
  colnames(malignant) = stringr::str_replace(colnames(malignant), "Cancer.cells", "Cancer")
  Cancer = cbind(Cancer, malignant)

  ##### Endothelial cells
  Endothelial = deconv[,grep("dothelial", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(Endothelial))]

  colnames(Endothelial) = stringr::str_replace(colnames(Endothelial), "Endothelial_cells", "Endothelial")
  colnames(Endothelial) = stringr::str_replace(colnames(Endothelial), "Endothelial.cells", "Endothelial")
  colnames(Endothelial) = stringr::str_replace(colnames(Endothelial), "Endothelial_cell", "Endothelial")

  ##### Eosinophils cells
  Eosinophils = deconv[,grep("osinophil", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(Eosinophils))]
  colnames(Eosinophils) = stringr::str_replace(colnames(Eosinophils), "Eosinophil(?!.)", "Eosinophils")

  ##### Plasma cells
  Plasma = deconv[,grep("lasma", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(Plasma))]

  colnames(Plasma) = stringr::str_replace(colnames(Plasma), "plasma(?!.)", "Plasma.cells")
  colnames(Plasma) = stringr::str_replace(colnames(Plasma), "Plasma_cells", "Plasma.cells")

  ##### Myocytes cells
  Myocytes = deconv[,grep("yocytes", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(Myocytes))]
  #colnames(Myocytes) = stringr::str_replace(colnames(Myocytes), "Myocytes", "Myocytes")

  ##### Fibroblasts cells :1 column
  Fibroblasts = data.frame(deconv[,grep("ibroblast", colnames(deconv))])
  colnames(Fibroblasts) = colnames(deconv)[grep("ibroblast", colnames(deconv))]
  deconv = deconv[,-which(colnames(deconv)%in%colnames(Fibroblasts))]
  #colnames(Fibroblasts) = stringr::str_replace(colnames(Fibroblasts), "Fibroblasts", "Fibroblasts")

  ##### Mast cells
  Mast = deconv[,grep("Mast", colnames(deconv))]
  Mast.activated = Mast[,grep("activated", colnames(Mast))]
  Mast.resting = Mast[,grep("resting", colnames(Mast))]
  Mast = Mast[,-which(colnames(Mast)%in%c(colnames(Mast.activated), colnames(Mast.resting)))]
  deconv = deconv[,-which(colnames(deconv)%in%c(colnames(Mast), colnames(Mast.activated), colnames(Mast.resting)))]

  colnames(Mast) = stringr::str_replace(colnames(Mast), "Mast_cell(?!.)", "Mast.cells")
  colnames(Mast) = stringr::str_replace(colnames(Mast), "Mast_cells", "Mast.cells")
  colnames(Mast.activated) = stringr::str_replace(colnames(Mast.activated), "Mast.cells.activated", "Mast.activated.cells")
  colnames(Mast.activated) = stringr::str_replace(colnames(Mast.activated), "Mast_cells_activated", "Mast.activated.cells")
  colnames(Mast.resting) = stringr::str_replace(colnames(Mast.resting), "Mast.cells.resting", "Mast.resting.cells")
  colnames(Mast.resting) = stringr::str_replace(colnames(Mast.resting), "Mast_cells_resting", "Mast.resting.cells")

  ##### B cells
  B.naive = deconv[,grep("naive", colnames(deconv))]
  B.memory = deconv[,grep("memory", colnames(deconv))]
  B = deconv[,-which(colnames(deconv)%in%c(colnames(B.naive), colnames(B.memory)))] #"B" will include also the cells haven't been re-named, as it is the last cell
  colnames(B) = stringr::str_replace(colnames(B), "B_cells", "B.cells")
  colnames(B) = stringr::str_replace(colnames(B), "B_cell", "B.cells")
  colnames(B) = stringr::str_replace(colnames(B), "B_lineage", "B.cells")
  colnames(B) = stringr::str_replace(colnames(B), "_B(?!.)", "_B.cells")
  colnames(B.naive) = stringr::str_replace(colnames(B.naive), "B.cells.naive", "B.naive.cells")
  colnames(B.naive) = stringr::str_replace(colnames(B.naive), "B_cells_naive", "B.naive.cells")
  colnames(B.naive) = stringr::str_replace(colnames(B.naive), "B_cell_naive", "B.naive.cells")
  colnames(B.memory) = stringr::str_replace(colnames(B.memory), "B.cells.memory", "B.memory.cells")
  colnames(B.memory) = stringr::str_replace(colnames(B.memory), "B_cells_memory", "B.memory.cells")
  colnames(B.memory) = stringr::str_replace(colnames(B.memory), "B_cell_memory", "B.memory.cells")


  cell_types = cbind(B, B.naive, B.memory, Macrophages, M0, M1, M2, Monocytes, Neutrophils, NK, NK.activated, NK.resting, NKT, CD4, CD4.memory.activated,
                     CD4.memory.resting, CD4.naive, CD4.non.regulatory, CD8, Tregs, Thelper, Tgamma, Dendritic, Dendritic.activated, Dendritic.resting, Cancer,
                     Endothelial, Eosinophils, Plasma, Myocytes, Fibroblasts, Mast, Mast.activated, Mast.resting, CAF)

  cat("Checking consistency in deconvolution cell fractions across patients...............................................................\n\n")

  combinations = c("Quantiseq", "MCP", "XCell", "Epidish_BPRNACan_",  "Epidish_BPRNACanProMet", "Epidish_BPRNACan3DProMet", "Epidish_CBSX.HNSCC.scRNAseq", "Epidish_CBSX.Melanoma.scRNAseq",
                   "Epidish_CBSX.NSCLC.PBMCs.scRNAseq", "Epidish_CCLE.TIL10", "Epidish_TIL10", "Epidish_LM22", "DeconRNASeq_BPRNACan_", "DeconRNASeq_BPRNACanProMet",
                   "DeconRNASeq_BPRNACan3DProMet", "DeconRNASeq_CBSX.HNSCC.scRNAseq", "DeconRNASeq_CBSX.Melanoma.scRNAseq", "DeconRNASeq_CBSX.NSCLC.PBMCs.scRNAseq", "DeconRNASeq_CCLE.TIL10",
                   "DeconRNASeq_TIL10", "DeconRNASeq_LM22", "CBSX_BPRNACan_", "CBSX_BPRNACanProMet", "CBSX_BPRNACan3DProMet", "CBSX_CBSX.HNSCC.scRNAseq",
                   "CBSX_CBSX.Melanoma.scRNAseq", "CBSX_CBSX.NSCLC.PBMCs.scRNAseq", "CBSX_CCLE.TIL10", "CBSX_TIL10", "CBSX_LM22", "DWLS_BPRNACan_",  "DWLS_BPRNACanProMet", "DWLS_BPRNACan3DProMet", "DWLS_CBSX.HNSCC.scRNAseq", "Epidish_CBSX.Melanoma.scRNAseq",
                   "DWLS_CBSX.NSCLC.PBMCs.scRNAseq", "DWLS_CCLE.TIL10", "DWLS_TIL10", "DWLS_LM22")

  for(i in combinations){
    mat = cell_types[,grep(i, colnames(cell_types))]
    print(paste("Total sum across samples of combination", i, "is", round(sum(mat[1, ]), 2)))
  }

  return(cell_types)
}
