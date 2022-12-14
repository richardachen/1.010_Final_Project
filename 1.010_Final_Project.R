library('pacman')
pacman::p_load(pacman, dplyr, ggplot2)

# PART A
#reading files: geeksforgeeks.org/how-to-import-data-from-a-file-in-r-programming/
path <- "/Users/richy/Downloads/Classes/Fall 2022/1.010/water-treatment.csv"
all_columns <- c("Date", "Input Flow", "Input Zinc", "Input pH", 
                 "Input BOD", "Input COD", "Input Suspended Solids", 
                 "Input Volatile Suspended Solids", "Input Sediment", 
                 "Input Conductivity", "Primary Settler pH", "Primary Settler BOD", 
                 "Primary Settler Suspended Solids", "Primary Settler Volatile Suspended Solids", 
                 "Primary Settler Sediments", "Primary Settler Conductivity", 
                 "Secondary Settler pH", "Secondary Settler BOD", 
                 "Secondary Settler COD", "Secondary Settler Suspended Solids", 
                 "Secondary Settler Volatile Suspended Solids", "Secondary Settler Sediments", 
                 "Secondary Settler Conductivity", "Output pH", 
                 "Output BOD", "Output COD", "Output Suspended Solids", 
                 "Output Volatile Suspended Solids", "Output Sediments", 
                 "Output Conductivity", "Primary Settler BOD Performance", 
                 "Primary Settler Suspended Solids Performance", "Primary Settler Sediment Performance",
                 "Secondary Settler BOD Performance", "Secondary Settler COD Performance", 
                 "Global BOD Performance", "Global COD Performance", 
                 "Global Suspended Solids Performance", "Global Sediments Performance")
water_qualities <- read.csv(path, col.names = all_columns)

#clean up the data slightly
remove_q_marks <- c()

for (i in names(water_qualities)) {
  for (a in water_qualities[[i]]) {
    if (a == "?"){
      remove_q_marks <- append(remove_q_marks, i)
      # https://www.delftstack.com/howto/r/break-a-for-loop-in-r/
      break
    }
  }
}  



#https://www.tutorialspoint.com/how-to-remove-a-column-from-an-r-data-frame
water_qualities <- water_qualities %>% select(-all_of(remove_q_marks))
water_qualities <- water_qualities %>% select(-Date)

# PART B

# binarize all data values
binarize <- function(vector){
  # function to binarize any vector, with the cutoff being the median.
  # input: vector of numbers
  # output: vector of same length as input, with 0's for below median, 
  # 1's for above median
  median_cutoff <- median(vector)
  binarized_form <- c()
  for (v in vector){
    binarized_form <- append(binarized_form, if (v <= median_cutoff){0} else {1})
  }
  return(binarized_form)
}

for (i in names(water_qualities)) {
  water_qualities[[i]] <- binarize(water_qualities[[i]])
}

# there are only 28 different values between input and primary settler conductivity
for (row in 1:nrow(water_qualities)){
  if (water_qualities[row, 2] != water_qualities[row, 5]){
    print('yes')
  }
}


# set up independence calculations
source('/Users/richy/Downloads/Classes/Fall 2022/1.010/pcalg.R')

subsets <- function(li){
  # function to generate a list of all subsets of a list
  # Input: li: a list
  # Output: a list of subsets of li
  if (length(li) == 0){
    return(list(list()))
  }
  else if (length(li) == 1){
    return(list(list(), list(li[[1]])))
  }
  else{
    original <- subsets(li[2:length(li)])
    more <- list()
    for (o in original){
      more <- append(more, list(append(li[1], o)))
    }
    
    return(append(original, more))
  }
}

# https://sparkbyexamples.com/r-programming/convert-vector-to-list-in-r/#:~:text=Convert%20Vector%20to%20List%20in%20R%20or%20create%20a%20list,contain%20elements%20of%20different%20types.

# PART C AND D
# find causation relations
check_connection <- function(idx_v,idx_w,df){
  # function to see if variables v, w are 
  # independent regardless of what variables it is conditioned on
  # Input: idx_v: index of v in names(df), idx_w: index of w in names(df), 
  # df: dataframe of values 
  # Output: tuple (tf,vec)
  # tf: True if v, w are always correlated, False otherwise.
  # vec: if tf is True, c(). if tf is False, vector of indexes 
  # of variables in df that, when conditioned on, make v, w independent.
  
  # https://www.statology.org/remove-element-from-vector-r/
  num_idx <- 1:length(names(df))
  num_idx <- num_idx[! num_idx %in% c(idx_v, idx_w)]
  # careful, possibilities has lists!
  possibilities <- subsets(as.list(num_idx))
  for (p in possibilities) {
    vec_subset <- unlist(p)
    p_val <- gSquareBin(idx_v, idx_w, vec_subset, as.matrix(df))
    if (p_val > 0.01){
      return(c(FALSE, vec_subset))
    }
  }
  return(c(TRUE, c()))
}

all_connections <- function(df){
  # function to carry out IC algorithm. For every pair of values in df, 
  # checks to see if they are connected.
  # Input: data frame with data values
  # Output: vector containing connections
  num_idx <- 1:length(names(df))
  connections <- list()
  for (i in num_idx[1:length(num_idx)-1]){
    for (j in num_idx[(i+1):length(num_idx)]){
      connect <- check_connection(i, j, df)
      if (connect[1] == TRUE){
        connections <- append(connections, list(list(i,j)))
      }
    }
  }
  return(connections)
}

connections <- all_connections(water_qualities)

# PART F

potential <- function(idx_v, idx_w, df){
  # function to see if variable v is a potential cause of variable w
  # Input: idx_v: index of v, idx_w: index of w, df: dataframe of values
  # Output: TRUE if v is a potential cause of w, FALSE if v is not a 
  # potential cause of w
  num_idx <- 1:length(df)
  num_idx <- num_idx[! num_idx %in% c(idx_v, idx_w)]
  for (z in num_idx){
    z_idx <- num_idx[! num_idx %in% c(z)]
    contexts <- subsets(z_idx)
    for (con in contexts){
      context <- unlist(con)
      v_z <- gSquareBin(idx_v, z, context, as.matrix(df))
      w_z <- gSquareBin(idx_w, z, context, as.matrix(df))
      if ((v_z > 0.01)&(w_z < 0.01)){
        return(TRUE)
      }
    }
  }
  return(FALSE)
}

genuine <- function(idx_v, idx_w, df){
  # function to see if variable v is a genuine cause of variable w
  # Input: idx_v: index of v, idx_w: index of w, df: dataframe of values
  # Output: TRUE if v is a genuine cause of w, FALSE if v is not a genuine 
  # cause of w
  num_idx <- 1:length(df)
  num_idx <- num_idx[! num_idx %in% c(idx_v, idx_w)]
  for (z in num_idx){
    z_idx <- num_idx[! num_idx %in% c(z)]
    contexts <- subsets(z_idx)
    if (potential(z, idx_v, df)){
      for (con in contexts){
        context <- unlist(con)
        w_z <- gSquareBin(idx_w, z, context, as.matrix(df))
        w_x_z <- gSquareBin(idx_w, z, c(context, idx_v), as.matrix(df))
        if ((w_z < 0.01)&(w_x_z > 0.01)){
          return(TRUE)
        }
      }
    }
  }
  
  return (FALSE)
}

for (connect in connections){
  pot <- potential(connect[[1]], connect[[2]], water_qualities)
  if (pot){
    print(paste(names(water_qualities)[connect[[1]]], 'is a potential cause of', names(water_qualities)[connect[[2]]], '.'))
  }
  gen <- genuine(connect[[1]], connect[[2]], water_qualities)
  if (gen){
    print(paste(names(water_qualities)[connect[[1]]], 'is a genuine cause of', names(water_qualities)[connect[[2]]], '.'))
  }
}

# PART G
ace <- function(idx_v, idx_w, df){
  # function to calculate the ACE of variable v on variable w.
  # Input: idx_v: the index of v, idx_w: the index of w, df: dataframe containing values
  # Output: the ACE of v on w
  # (v,w) first position: (0,0), second: (1,0), third: (0,1), fourth: (1,1)
  # P(W|V)
  confusion_mat <- list(0,0,0,0)
  
  for (row in 1:nrow(water_qualities)){
    plus <- water_qualities[row, idx_v] + 2*water_qualities[row, idx_w] + 1
    confusion_mat[plus] <- confusion_mat[[plus]] + 1
  }
  
  ace <- confusion_mat[[4]]/(confusion_mat[[2]]+confusion_mat[[4]]) - confusion_mat[[3]]/(confusion_mat[[1]]+confusion_mat[[3]])
  return(ace)
}

for (i in 1:(length(water_qualities)-1)){
  for (j in (i+1):length(water_qualities)){
    print(paste('The ACE of ', names(water_qualities)[i], ' on ', names(water_qualities)[j], ' is ', ace(i,j,water_qualities), '.'))
  }
}


# PART I

all_one <- function(indices, exclude, row){
  # function to check if a row of a dataframe has indices 1 and everything else 0
  # Input: indices: a vector of indices, row: a row of a dataframe
  # Output: TRUE if matching, FALSE if not
  
  for (i in 1:length(row)){
    if (i %in% indices){
      if (row[[i]] == 0){
        return(FALSE)
      }
    }
    else if (!(i %in% exclude)){
      if (row[[i]] == 1){
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

multi_variate <- function(idx_v, idx_w, df){
  # function to calculate the multivariate linear regression between v and w
  # Input: idx_v: the index of v, idx_w: the index of w, df: dataframe containing values
  # Output: the multivariate linear regression between v and w
  
  num_idx <- 1:length(df)
  num_idx <- num_idx[! num_idx %in% c(idx_v, idx_w)]
  contexts <- subsets(num_idx)
  # Cumulative sum for v=1
  v_1 <- 0
  # Cumulative sum for v=0
  v_0 <- 0
  
  # How many rows have 1 in the v spot
  total_v <- 0  
  for (row in 1:nrow(df)){
    if (df[row, idx_v] == 1){
      total_v <- total_v + 1
    }
  }

  # Sum over contexts 
  for (con in contexts){
    context <- unlist(con)
    # https://www.statology.org/r-select-rows-by-condition/
    confusion_mat_v_w <- list(0,0,0,0)
    row_count <- 0
    v_count <- 0
    
    for (row in 1:nrow(df)){
      if (all_one(context, c(idx_v, idx_w), df[row, ])){
        plus <- df[row, idx_v] + 2*df[row, idx_w] + 1
        confusion_mat_v_w[plus] <- confusion_mat_v_w[[plus]] + 1
        
        row_count <- row_count + 1
        if (df[row, idx_v] == 1){
          v_count <- v_count + 1
        }
      }
    }
    
    if ((v_count != 0)&(total_v != 0)){
      # Sum with v=1
      v_1 <- v_1 + (confusion_mat_v_w[[4]]/v_count)*((confusion_mat_v_w[[2]] + confusion_mat_v_w[[4]])/total_v)
    }
    if (((row_count-v_count) != 0)&((nrow(df)-total_v) != 0)){
      # Sum with v=0
      v_0 <- v_0 + (confusion_mat_v_w[[3]]/(row_count-v_count))*((confusion_mat_v_w[[1]] + confusion_mat_v_w[[3]])/(nrow(df)-total_v))
    }
  }
  return(v_1-v_0)
} 

for (i in 1:(length(water_qualities)-1)){
  for (j in (i+1):length(water_qualities)){
    print(paste('The multivariate linear regression of ', names(water_qualities)[i], ' on ', names(water_qualities)[j], ' is ', multi_variate(i,j,water_qualities), '.'))
  }
}

#clean up everything at the end
pacman::p_unload(all)
