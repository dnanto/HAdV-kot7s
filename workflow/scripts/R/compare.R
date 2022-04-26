#!/usr/bin/env Rscript --vanilla

main <- function(snakemake) {
  library(tidyverse)
  library(BactDating)

  thresh.burn <- snakemake@params$burn
  thresh.ess <- snakemake@params$ess

  key <- c("likelihood", "mu", "sigma", "alpha", "prior")
  df.diagnostic <- (
    lapply(snakemake@input, function(path) {
      run <- qs::qread(path)
      key <- c("likelihood", "mu", "sigma", "alpha", "prior")
      rec <- with(run, record[max(1, round(nrow(record) * thresh.burn)):nrow(record), key])
      est <- (
        apply(rec[, key], 2, summary) %>%
          t() %>%
          apply(1, function(row) sprintf("%.2e [%.2e;%.2e]", row["Median"], row["1st Qu."], row["3rd Qu."])) %>%
          setNames(paste0("est.", key))
      )
      ess <- setNames(coda::effectiveSize(rec), paste0("ess.", key))
      c(
        path = path,
        name = basename(path),
        model = as.character(run$model),
        seed = run$seed,
        dic = run$dic,
        mu = run$mu,
        rootprob = run$rootprob,
        est,
        ess
      )
    }) %>%
      bind_rows() %>%
      mutate(
        pass = ess.likelihood >= thresh.ess & ess.mu >= thresh.ess & ess.alpha >= thresh.ess & ess.prior >= thresh.ess,
        pass = pass & ((model %in% c("poisson", "strictgamma")) | (ess.sigma >= thresh.ess))
      ) %>%
      arrange(desc(pass), desc(dic), across(starts_with("ess."), desc))
  )

  select(df.diagnostic, -path) %>% write_tsv(snakemake@output[["tsv"]])
  pull(df.diagnostic, path) %>%
    head(1) %>%
    stringr::str_replace(".qs$", ".json") %>%
    file.copy(snakemake@output[["obj"]]) ->
  devnull
  pull(df.diagnostic, path) %>%
    head(1) %>%
    stringr::str_replace(".qs$", ".tree") %>%
    file.copy(snakemake@output[["tree"]]) ->
  devnull
}

capture.output(
  capture.output(main(snakemake), file = snakemake@log[["out"]], type = "output"),
  file = snakemake@log[["msg"]], type = "message"
)
