calculate_dendrogram_cuts = function(cell.group.dendrogram, n_cuts = NULL){

  cuts = list()
  for (i in 1:length(cell.group.dendrogram)) {
    dend_heights <- dendextend::heights_per_k.dendrogram(as.dendrogram(cell.group.dendrogram[[i]])) #Calculate dendrogram heights

    sorted_heights <- sort(dend_heights) # Sort heights
    height_diffs <- diff(sorted_heights) # Calculate differences between consecutive heights

    buffer <- max(median(height_diffs[height_diffs > 0]), 1) # Buffer: take the median of the non-zero differences

    # Add and rest the buffer to the minimum and maximum height respectively to avoid trivial cuts (clusters with 1 feature and cluster with all features)
    min_height <- min(sorted_heights) + buffer
    max_height <- max(sorted_heights) - buffer

    # Ensure the buffer-adjusted min height is valid
    if (min_height >= max_height) {
      stop("Buffered minimum height exceeds maximum height for dendrogram")
    }

    if(is.null(n_cuts)==T){
      number_cuts = floor(max_height) #Truncate based on highest height of dendrogram
    }else{
      number_cuts = n_cuts
    }

    cut_sequence <- seq(min_height, max_height, length.out = number_cuts)  # Generate a sequence of cut heights between the buffered min and max
    #cut_sequence = cut_sequence[-c(1, length(cut_sequence))] #Remove first and last cut to avoid trivial groups (either clusters of 1 feature or a single cluster with all of them)
    cut_sequence <- round(cut_sequence, 2)  # Round the cut heights for cleaner values (2 decimal place)

    cuts[[i]] = cut_sequence
  }

  #Give format to the list to have a sequence of cuts per dendrogram
  if(is.null(n_cuts)==F){
    #Adjust to return the combinations of cuts across all dendrograms
    combined_cuts <- matrix(NA, nrow = n_cuts, ncol = length(cell.group.dendrogram))

    for (i in 1:length(cuts)) {
      combined_cuts[1:length(cuts[[i]]), i] <- cuts[[i]]
    }

    combined_cuts_list <- split(combined_cuts, row(combined_cuts))
    return(combined_cuts_list)
  }else{
    return(cuts)
  }


}
