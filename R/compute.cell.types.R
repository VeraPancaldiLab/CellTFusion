#' Cell types split from deconvolution
#'
#' @param data A matrix with the deconvolution results
#'
#' @return A list containing:
#' - A sublist with each cell type features as an element recover from the different signatures
#' - Discarded cell types (this will happen if the cell types are not supported. See the READme for more information about this)
#' @export
#'
#' @examples
#'
#'   cells_types = compute.cell.types(deconvolution)
#'   cells = cells_types[[1]]
#'   cells_discarded = cells_types[[2]]
#'
compute.cell.types = function(data){
  ##### B cells
  B = grep("B.cells", colnames(data))
  B = data[, B, drop = FALSE]
  ##### B naive
  B.naive = grep("B.naive", colnames(data))
  B.naive = data[, B.naive, drop = FALSE]
  ##### B memory
  B.memory = grep("B.memory", colnames(data))
  B.memory = data[, B.memory, drop = FALSE]
  ##### Macrophages (M0, M1, M2)
  Macrophages = grep("Macrophages.cells", colnames(data))
  Macrophages = data[, Macrophages, drop = FALSE]
  M0 = grep("Macrophages.M0", colnames(data))
  M0 = data[, M0, drop = FALSE]
  M1 = grep("Macrophages.M1", colnames(data))
  M1 = data[, M1, drop = FALSE]
  M2 = grep("Macrophages.M2", colnames(data))
  M2 = data[, M2, drop = FALSE]
  ##### Monocytes
  Monocytes = grep("Monocytes", colnames(data))
  Monocytes = data[, Monocytes, drop = FALSE]
  ##### Neutrophils
  Neutrophils = grep("Neutrophils", colnames(data))
  Neutrophils = data[, Neutrophils, drop = FALSE]
  ##### NK cells (activated, resting)
  NK = grep("NK.cells", colnames(data))
  NK = data[, NK, drop = FALSE]
  NK.activated = grep("NK.activated", colnames(data))
  NK.activated = data[, NK.activated, drop = FALSE]
  NK.resting = grep("NK.resting", colnames(data))
  NK.resting = data[, NK.resting, drop = FALSE]
  ##### NKT cells
  NKT = grep("NKT.cells", colnames(data))
  NKT = data[, NKT, drop = FALSE]
  ##### CD4 cells (activated, resting)
  CD4 = grep("CD4.cells", colnames(data))
  CD4 = data[, CD4, drop = FALSE]
  #memory = CD4[,grep("memory", colnames(CD4))]
  #helper = CD4[,grep("Th", colnames(CD4))]
  #CD4 = CD4[,-which(colnames(CD4)%in%c(colnames(memory), colnames(helper)))]
  CD4.memory.activated = grep("CD4.memory.activated", colnames(data))
  CD4.memory.activated = data[, CD4.memory.activated, drop = FALSE]
  CD4.memory.resting = grep("CD4.memory.resting", colnames(data))
  CD4.memory.resting = data[, CD4.memory.resting, drop = FALSE]
  CD4.naive = grep("CD4.naive", colnames(data))
  CD4.naive = data[, CD4.naive, drop = FALSE]
  ##### CD8 cells
  CD8 = grep("CD8.cells", colnames(data))
  CD8 = data[, CD8, drop = FALSE]
  #naive = grep("naive", colnames(CD8))
  #naive = CD8[, naive, drop = FALSE]
  #memory = grep("memory", colnames(CD8))
  #memory = CD8[, memory, drop = FALSE]
  #CD8 = CD8[,-which(colnames(CD8)%in%c(colnames(memory), colnames(naive)))]
  ##### Regulatory T cells
  Tregs = grep("T.cells.regulatory", colnames(data))
  Tregs = data[, Tregs, drop = FALSE]
  ##### Non regulatory T cells
  T.non.regs = grep("T.cells.non.regulatory", colnames(data))
  T.non.regs = data[, T.non.regs, drop = FALSE]
  ##### Helper T cells
  Thelper = grep("T.cells.helper", colnames(data))
  Thelper = data[, Thelper, drop = FALSE]
  ##### Gamma delta T cells
  Tgamma = grep("T.cells.gamma.delta", colnames(data))
  Tgamma = data[, Tgamma, drop = FALSE]
  ##### Dendritic cells (activated, resting)
  Dendritic = grep("Dendritic.cells", colnames(data))
  Dendritic = data[, Dendritic, drop = FALSE]
  Dendritic.activated = grep("Dendritic.activated", colnames(data))
  Dendritic.activated = data[, Dendritic.activated, drop = FALSE]
  Dendritic.resting = grep("Dendritic.resting", colnames(data))
  Dendritic.resting = data[, Dendritic.resting, drop = FALSE]
  ##### Cancer cells
  Cancer = grep("Cancer", colnames(data))
  Cancer = data[, Cancer, drop = FALSE]
  ##### Endothelial cells
  Endothelial = grep("Endothelial", colnames(data))
  Endothelial = data[, Endothelial, drop = FALSE]
  ##### Eosinophils cells
  Eosinophils = grep("Eosinophils", colnames(data))
  Eosinophils = data[, Eosinophils, drop = FALSE]
  ##### Plasma cells
  Plasma = grep("Plasma.cells", colnames(data))
  Plasma = data[, Plasma, drop = FALSE]
  ##### Myocytes cells
  Myocytes = grep("Myocytes", colnames(data))
  Myocytes = data[, Myocytes, drop = FALSE]
  ##### Fibroblasts cells
  Fibroblasts = grep("Fibroblasts", colnames(data))
  Fibroblasts = data[, Fibroblasts, drop = FALSE]
  ##### Mast cells
  Mast = grep("Mast.cells", colnames(data))
  Mast = data[, Mast, drop = FALSE]
  Mast.activated = grep("Mast.activated", colnames(data))
  Mast.activated = data[, Mast.activated, drop = FALSE]
  Mast.resting = grep("Mast.resting", colnames(data))
  Mast.resting = data[, Mast.resting, drop = FALSE]
  ##### CAF cells
  CAF = grep("CAF", colnames(data))
  CAF = data[, CAF, drop = FALSE]

  #####Output list
  cell_types = list(B, B.naive, B.memory, Macrophages, M0, M1, M2, Monocytes, Neutrophils, NK, NK.activated, NK.resting, NKT, CD4, CD4.memory.activated, CD4.memory.resting, CD4.naive,
                    CD8, Tregs, T.non.regs, Thelper, Tgamma, Dendritic, Dendritic.activated, Dendritic.resting, Cancer, Endothelial, Eosinophils, Plasma, Myocytes, Fibroblasts, Mast, Mast.activated,
                    Mast.resting, CAF)

  names(cell_types) = c("B.cells", "B.naive", "B.memory", "Macrophages.cells", "Macrophages.M0", "Macrophages.M1", "Macrophages.M2", "Monocytes", "Neutrophils", "NK.cells", "NK.activated", "NK.resting", "NKT.cells", "CD4.cells", "CD4.memory.activated",
                        "CD4.memory.resting", "CD4.naive", "CD8.cells", "T.cells.regulatory", "T.cells.non.regulatory","T.cells.helper", "T.cells.gamma.delta", "Dendritic.cells", "Dendritic.activated", "Dendritic.resting", "Cancer", "Endothelial",
                        "Eosinophils", "Plasma.cells", "Myocytes", "Fibroblasts", "Mast.cells", "Mast.activated", "Mast.resting", "CAF")

  ####Discarded cell types
  cell_types_matrix = cbind(B, B.naive, B.memory, Macrophages, M0, M1, M2, Monocytes, Neutrophils, NK, NK.activated, NK.resting, NKT, CD4, CD4.memory.activated, CD4.memory.resting, CD4.naive,
                            CD8, Tregs, T.non.regs, Thelper, Tgamma, Dendritic, Dendritic.activated, Dendritic.resting, Cancer, Endothelial, Eosinophils, Plasma, Myocytes, Fibroblasts, Mast, Mast.activated,
                            Mast.resting, CAF)

  cell_types_discarded = data[,!(colnames(data)%in%colnames(cell_types_matrix)), drop = F]

  return(list(cell_types, cell_types_discarded))

}
