

calculate_mz <- function(mass, charge) {
  
  mz <- mass + (charge * hydrogen_mass) / charge
  return(mz)
}


frag_masses <- function(sequence, frag_types, min_size = 1, charges_considered = c(1) ) {
  # Create an empty list to store the substrings
  hydrogen_mass <- 1.00728
  masses <- c()
  
  for (i in 1:nchar(sequence)) {
    for (j in i:nchar(sequence)) {
      substr <- substr(sequence, i, j)
      if (nchar(substr) >= min_size && nchar(substr) < nchar(sequence)) {
        
        #Terminal fragments
        if(i == 1 || j == nchar(sequence)){
          ic = ion_caps %>% filter(internal == F) %>% filter(type %in% frag_types)
          for(r in 1:nrow(ic)){
            for(c in charges_considered){
              masses <- c(masses, (mw(substr)+ic[r,]$delta_mass + ( c * hydrogen_mass)) / c )
            }
          }
        }else{#Internal fragments
          ic = ion_caps %>% filter(internal == T) %>% filter(type %in% frag_types)
          if(nrow(ic) >0 ){
            for(r in 1:nrow(ic)){
              for(c in charges_considered){
                masses <- c(masses, (mw(substr)+ic[r,]$delta_mass + ( c * hydrogen_mass)) / c )
              }
            }
          }
        }
      }
    }
  }
  # Return the list of substrings
  

  return(masses)
}


error_intervals <- function(values, error_ppm) {
  # Calculate the error in absolute terms for each value
  errors <- values * error_ppm / 1e6
  
  # Create a list to store the error intervals
  error_intervals <- list()
  
  # Calculate the error intervals for each value
  for (i in 1:length(values)) {
    value <- values[i]
    error <- errors[i]
    interval <- c(value - error, value + error)
    error_intervals[[i]] <- interval
  }
  
  return(error_intervals)
}


generate_peptide_sequence <- function(n) {
  # Create a vector of the 20 amino acids
  amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                   "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  
  # Create a vector of the relative frequencies of each amino acid in proteins
  amino_acid_frequencies <- c(0.0825, 0.0139, 0.0547, 0.0675, 0.0407,
                              0.0684, 0.0228, 0.0586, 0.0585, 0.0944,
                              0.0228, 0.0411, 0.0473, 0.0522, 0.0549,
                              0.0871, 0.0595, 0.0682, 0.0119, 0.0319)
  
  # Generate a random sequence of amino acids
  sequence <- sample(amino_acids, size = n, prob = amino_acid_frequencies, replace = TRUE)
  
  # Join the amino acids into a single string
  sequence <- paste(sequence, collapse = "")
  
  return(sequence)
}

coverage <- function(segments, range) {
  
  segments =   segments[order(sapply(segments, "[", 1))]
  totalInterval <- 0
  
  for (i in 1:length(segments)) {
    if (segments[[i]][1] < range[1]) {
      segments[[i]][1] <- range[1]
    }
    if (segments[[i]][2] > range[2]) {
      segments[[i]][2] <- range[2]
    }
    
    currentIntrval <- segments[[i]][2] - segments[[i]][1]
    
    if (i == 1) {
      if (segments[[i]][1] < range[1]) {
        differenceFromPrevious <- range[1] - segments[[i]][1]
      } else {
        differenceFromPrevious <- 0
      }
    } else {
      if (segments[[i]][1] < segments[[i - 1]][2]) {
        differenceFromPrevious <- segments[[i - 1]][2] - segments[[i]][1]
      } else {
        differenceFromPrevious <- 0
      }
    }
    
    totalInterval <- totalInterval + currentIntrval - differenceFromPrevious
    
    if (i == length(segments)) {
      if (range[2] < segments[[i]][2]) {
        differenceFromPrevious <- segments[[i]][2] - range[2]
        totalInterval <- totalInterval + currentIntrval - differenceFromPrevious
      }
    }
  }
  
  return (totalInterval / (range[2] - range[1]))
}




calculate_overlapping <- function(value_list, error_tolerance) {
  num_values <- length(value_list)
  overlap_no <- c()
  
  for (i in 1:num_values) {
    value <- value_list[i]
    overlapping <- value_list[abs(value_list - value) <= (value * error_tolerance / 1e6)]
    overlap_no = c(overlap_no, length(overlapping))
  }
  return(data.frame(table(overlap_no)))
}