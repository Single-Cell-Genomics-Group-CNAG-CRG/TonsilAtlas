
# Objective
# Convert the outputs from Albert to bedgraph (profiles and candidates)

# Parameters
path_to_dir_profile <- here::here("scATAC-seq/Cicero/results/acProfile/profile/")
path_to_dir_candidates <- here::here("scATAC-seq/Cicero/results/acProfile/candidates/")

# Function1
bed2bedgraph <- function(bed_file, path_to_dir = path_to_dir_profile)
  {
  print(bed_file)
  bed <- read.table(paste0(path_to_dir,bed_file))
  bed$acScore [is.na(bed$acScore)] <- 0
  write.table(bed[,c(1,2,3,6)], 
              row.names = F, 
              col.names = F,
              quote = F, 
              file = paste0(path_to_dir,gsub(".txt",".bedGraph",bed_file)))
  }

# Load objects

# Profiles
lapply(list.files(path_to_dir_profile, pattern=".txt"),bed2bedgraph)


# Function2
bed2bedgraph <- function(bed_file, path_to_dir = path_to_dir_candidates)
{
  print(bed_file)
  bed <- read.table(paste0(path_to_dir,bed_file))
  bed$acScore [is.na(bed$acScore)] <- 0
  write.table(bed[,c(1,2,3,6)], 
              row.names = F, 
              col.names = F,
              quote = F, 
              file = paste0(path_to_dir,gsub(".txt",".bedGraph",bed_file)))
}

# Candidates
lapply(list.files(path_to_dir_candidates, pattern=".txt"),bed2bedgraph)