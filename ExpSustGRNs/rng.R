importar_red <- function(fuente) {
  setwd(fuente)
  network <- read.table("originalnetwork.txt")
  adjacency_matrix <- as.matrix(network)
  return(adjacency_matrix)
}

cuantos_nodos <- function(adjacency_matrix) {
  num_nodes <- nrow(adjacency_matrix)
  return(num_nodes)
}

cuantas_interacciones <- function(out_degrees) {
  num_interacciones <- sum(out_degrees)
  return(num_interacciones)
}

cuantas_negativas <- function(adjacency_matrix, num_nodes) {
  cuantas_negativas <- 0L
  for (i in 1:num_nodes) {
    for (j in 1:num_nodes) {
      if (adjacency_matrix[i, j] < 0 && i != j) {
        cuantas_negativas <- cuantas_negativas + 1
      }
    }
  }
  for (i in 1:num_nodes) {
    j <- i
    if (adjacency_matrix[i, j] < 0) {
      cuantas_negativas <- cuantas_negativas + 1
    }
  }
  return(cuantas_negativas)
}

matriz_positiva <- function(adjacency_matrix, num_nodes) {
  adjacency_matrix_positive <- matrix(0, nrow = num_nodes, ncol = num_nodes)
  adjacency_matrix_positive <- adjacency_matrix
  for (i in 1:num_nodes) {
    for (j in 1:num_nodes) {
      if (adjacency_matrix_positive[i, j] < 0) {
        adjacency_matrix_positive[i, j] <- 1
      }
    }
  }
  return(adjacency_matrix_positive)
}

obtain_out_degrees <- function(adjacency_matrix_positive) {
  out_degrees <- rowSums(adjacency_matrix_positive)
  return(out_degrees)
}

obtain_in_degrees <- function(adjacency_matrix_positive) {
  in_degrees <- colSums(adjacency_matrix_positive)
  return(in_degrees)
}

eliminate_once <- function(vector, sacar) {
  for (elemento in sacar) {
    pos <- match(elemento, vector)
    if (!is.na(pos)) {
      vector <- vector[-pos]
    }
  }
  return(vector)
}

has_duplicates <- function(vec) {
  return(length(vec) != length(unique(vec)))
}

generate_random_adjacency <- function(adjacency_matrix, in_degrees, out_degrees) {
  num_nodes <- cuantos_nodos(adjacency_matrix)
  num_interacciones <- cuantas_interacciones(out_degrees)
  la_matrix <- matrix(0, nrow = num_nodes, ncol = num_nodes)
  hay_lugar <- lapply(1:num_nodes, function(x) rep(x, in_degrees[x]))
  hay_lugar <- unlist(hay_lugar)
  i <- 1
  contador_repeat <- 1
  contador_unique <- 1
  
  while (i <= num_nodes) {
    if (out_degrees[i] == 0) {
      i <- i + 1
      next
    }
    reiniciar <- FALSE
    unique_lugar <- unique(hay_lugar)
    if (length(hay_lugar) < out_degrees[i]) {
      stop("Híjole, ya no hay espacio, carnal ", i)
    }
    if (length(unique_lugar) < out_degrees[i]) {
      contador_unique <- contador_unique + 1
      la_matrix <- matrix(0, nrow = num_nodes, ncol = num_nodes)
      hay_lugar <- c()
      hay_lugar <- lapply(1:num_nodes, function(x) rep(x, in_degrees[x]))
      hay_lugar <- unlist(hay_lugar)
      contador_repeat <- 1
      reiniciar <- TRUE
      i <- 1
      next
    }
    repeat {
      contador_repeat <- contador_repeat + 1
      targets <- sample(hay_lugar, out_degrees[i])
      if (!has_duplicates(targets)) {
        break
      }
    }
    if (!reiniciar) {
      hay_lugar <- eliminate_once(hay_lugar, targets)
      unique_lugar <- unique(hay_lugar)
      for (target in targets) {
        la_matrix[i, target] <- 1
      }
      i <- i + 1
    }
  }
  
  interacciones_negativas <- cuantas_negativas(adjacency_matrix, num_nodes)
  entradas_interaccion <- list()
  cont_interacciones <- 1
  for (i in 1:num_nodes) {
    for (j in 1:num_nodes) {
      if (la_matrix[i, j] > 0 && i != j) {
        entradas_interaccion[[cont_interacciones]] <- c(i, j)
        cont_interacciones <- cont_interacciones + 1
      }
    }
  }
  for (i in 1:num_nodes) {
    j <- i
    if (la_matrix[i, j] > 0) {
      entradas_interaccion[[cont_interacciones]] <- c(i, j)
      cont_interacciones <- cont_interacciones + 1
    }
  }
  
  if ((cont_interacciones - 1) != num_interacciones) {
    print("ERROR FATAL ¡ARCARSENAL!")
  }
  
  if (interacciones_negativas > 0) {
    di_cuales <- sample(1:num_interacciones, interacciones_negativas, replace = FALSE)
    for (i in 1:interacciones_negativas) {
      en_cual_ponerla <- entradas_interaccion[[di_cuales[i]]]
      la_matrix[en_cual_ponerla[1], en_cual_ponerla[2]] <- -1
    }
  }
  return(la_matrix)
}

exportar_red <- function(la_matrix, r) {
  file_name <- paste("network", r, ".txt", sep = "") #aquí la r es el contador de red, dentro de un ciclo for (r in 1:1000)
  write.table(la_matrix, file = file_name, row.names = FALSE, col.names = FALSE)
}