#' simulate peptide level data
#' @export
#' @param Nprot number of porteins
#' @param probability_of_success
#' @examples
#'
#' res <- sim_data()
#' View(res)
sim_data <- function(Nprot = 1000,
                     probability_of_success = 0.3,
                     fc = c(D = -2, N = 0, U = 2),
                     prop = c(D = 10, N = 80, U = 10),
                     mean_prot = 20,
                     sd = 1.2
                     ) {
  proteins <- stringi::stri_rand_strings(Nprot,6)

  # simulate number of peptides per protein
  nrpeptides <- rgeom(n = Nprot, probability_of_success) + 1

  prot <- data.frame(
    proteinID = proteins,
    nrPeptides = nrpeptides,
    FC = rep(fc, prop / 100 * Nprot),
    DE = rep(names(fc), prop / 100 * Nprot),
    sd = rep(1, Nprot),
    average_prot_abundance = rlnorm(Nprot,log(20),sdlog = log(sd))
  )
  prot <- prot |> tidyr::unite("proteinID", proteinID, DE)
  prot <- prot |> dplyr::mutate(mean_A = 0 + FC, N_A = 4, mean_B = 0, N_B = 4)

  # add row for each protein
  peptide_df <- prot |> uncount(nrPeptides)
  # create peptide ids
  peptide_df$peptideID <- stringi::stri_rand_strings(sum(prot$nrPeptides),8)

  # transform into long format
  peptide_df2 <- peptide_df |> pivot_longer(cols = starts_with(c("mean", "N_")), names_to = "group" , values_to = "mean")
  peptide_df2 <-  peptide_df2 |> tidyr::separate(group, c("what", "group"))
  peptide_df2 <- peptide_df2 |> pivot_wider(names_from = "what", values_from = mean)

  peptide_df2$avg_peptide_abd <-
    with(peptide_df2,
         rlnorm(nrow(peptide_df2),
                meanlog = log(average_prot_abundance - mean),
                sdlog = log(sd)))
  sample_from_normal <- function(mean, sd, n) {
    rnorm(n, mean, sd)
  }

  sampled_data <- with(
    peptide_df2 ,
    mapply(sample_from_normal, avg_peptide_abd, sd, N)
  )
  sampled_data <- as.data.frame(t(sampled_data))
  x <- dplyr::bind_cols(peptide_df2,sampled_data)


  peptideAbudances <- x |>
    tidyr::pivot_longer(
      tidyselect::starts_with("V"),
      names_to = "Replicate",
      values_to = "abundance")
  peptideAbundances <- peptideAbudances |>
    tidyr::unite("sample", group, Replicate, remove =  FALSE)
  return(peptideAbundances)

}

