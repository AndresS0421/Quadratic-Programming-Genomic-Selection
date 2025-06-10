# MEASURE THE SIMILARITIES -------------------------------------------------------
# Cosine similarity - Direction of the vectors
cosine_similarity <- function(vec1, vec2) {
  raw_result <- sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2)))
  normalized_result <- (raw_result + 1)/2
  return(normalized_result)
}

# Euclidean distance - Measures the straight-line distance between two points
max_euclidean <- function(data) {
  # Assume data is a matrix or data frame with 2 columns (2D vectors)
  max_dist <- 0
  for(i in 1:(nrow(data) - 1)) {
    for(j in (i + 1):nrow(data)) {
      dist <- sqrt(sum((data[i, ] - data[j, ])^2))
      if(dist > max_dist) max_dist <- dist
    }
  }
  return(max_dist)
}
max_euclidean_value <- max_euclidean(Geno)

euclidean_distance <- function(vec1, vec2) {
  raw_result <- sqrt(sum((vec1 - vec2)^2))
  normalized_result <- 1 - (raw_result / max_euclidean_value)
  return(normalized_result)
}

# Pearson correlation - Measures the linear relationship, using covariance, and sd
pearson_correlation <- function(vec1, vec2) {
  raw_result <- cor(vec1, vec2, method = "pearson")
  normalized_result <- (raw_result + 1)/2
  return(normalized_result)
}

# Manhattan distance - Measures the sum of absolute differences
max_manhattan <- function(data) {
  # Assume data is a matrix or data frame with 2 columns (2D vectors)
  max_dist <- 0
  for(i in 1:(nrow(data) - 1)) {
    for(j in (i + 1):nrow(data)) {
      dist <- sum(abs(data[i, ] - data[j, ]))
      if(dist > max_dist) max_dist <- dist
    }
  }
  return(max_dist)
}
max_manhattan_value <- max_manhattan(Geno)

manhattan_distance <- function(vec1, vec2) {
  raw_result <- sum(abs(vec1 - vec2))
  normalized_result <- 1 - (raw_result / max_manhattan_value)
  return(normalized_result)
}

line_1 <- Geno[1, ]
line_2 <- Geno[2, ]

print(paste("Cosine:", cosine_similarity(line_1, line_2), "Euclidean:", euclidean_distance(line_1, line_2), "Pearson:", pearson_correlation(line_1, line_2), "Manhattan:", manhattan_distance(line_1, line_2)))