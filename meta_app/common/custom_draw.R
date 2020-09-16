#edited version of draw_forest() from MetaVolcanoR

draw_forest2 <- function (remres, gene = "MMP9", genecol = "Symbol", foldchangecol = "Log2FC", 
          llcol = "CI.L", rlcol = "CI.R", jobname = "MetaVolcano", 
          outputfolder = ".", draw = "") 
{
  if (!draw %in% c("PDF", "HTML", "")) {
    stop("Oops! Seems like you did not provide a right 'draw' parameter. \n              Try 'PDF' or 'HTML'")
  }
  if (is(remres) != "MetaVolcano") {
    stop("Oops! Please, provide a MetaVolcano object as input")
  }
  rem <- sremres <- merge(remres@metaresult, remres@input, 
                          by = genecol) %>% dplyr::filter(!!rlang::sym(genecol) == 
                                                            gene)
  if (nrow(sremres) == 0) {
    stop(paste("Oops! Seems that", gene, "is not in the", 
               "provided REM result"))
  }
  stds <- unique(unlist(regmatches(colnames(sremres), regexec("_\\d+$", 
                                                              colnames(sremres)))))
  if (is.null(remres@inputnames)) {
    message("We recomend providing a character vector with the names\n\t        of the input studies")
    stds <- setNames(stds, paste("study_", seq_along(stds)))
  }
  else {
    stds <- setNames(stds, remres@inputnames)
  }
  edat <- Reduce(rbind, lapply(names(stds), function(sn) {
    std <- dplyr::select(sremres, dplyr::matches(paste0(genecol, 
                                                        "|", stds[sn], "$")))
    colnames(std) <- gsub("_\\d+$", "", colnames(std))
    std[["group"]] <- sn
    std
  }))
  if (!all(c(genecol, foldchangecol, llcol, rlcol) %in% colnames(edat))) {
    stop("Oops! Please, check the match among the provided parameters\n              and the colnames of the remres@metaresult and remres@input")
  }
  edat <- dplyr::select(edat, c(!!rlang::sym(genecol), !!rlang::sym(foldchangecol), 
                                !!rlang::sym(llcol), !!rlang::sym(rlcol), group))
  sdat <- data.frame(genecol = unique(edat[[genecol]]), foldchangecol = sremres[["randomSummary"]], 
                     llcol = sremres[["randomCi.lb"]], rlcol = sremres[["randomCi.ub"]], 
                     group = "FoldChange summary")
  colnames(sdat) <- c(genecol, foldchangecol, llcol, rlcol, 
                      "group")
  dat <- rbind(edat, sdat)
  dat[["class"]] <- ifelse(grepl("summary", dat[["group"]]), 
                           "FoldChange summary", "Study")
  sumfc <- dplyr::filter(dat, grepl("summary", class))[[foldchangecol]]
  maxfc <- max(dat[[rlcol]])
  minfc <- min(dat[[llcol]])
  if (sumfc > 0) {
    sumcol <- "#E41A1C"
    minlim <- -maxfc
    maxlim <- maxfc
  }
  else {
    sumcol <- "#377EB8"
    minlim <- minfc
    maxlim <- -minfc
  }
  gg <- ggplot(dat, aes(x = group, y = !!rlang::sym(foldchangecol), 
                        color = class)) + geom_point() + geom_errorbar(aes(ymin = !!rlang::sym(llcol), 
                                                                           ymax = !!rlang::sym(rlcol), width = 0.1, color = class)) + 
    scale_color_manual(values = c(sumcol, "#bdbdbd")) + 
    scale_x_discrete(limits = rev(dat[["group"]])) + theme_classic() + 
    ggtitle(unique(edat[[genecol]])) + geom_hline(yintercept = 0, 
                                                  linetype = "solid", size = 0.05, color = "#969696") + 
    geom_hline(yintercept = sumfc, linetype = "dashed", 
               size = 0.1, color = sumcol) + theme(legend.position = "none") + 
    scale_y_continuous(limits = c(minlim, maxlim)) + coord_flip()
  if (draw == "PDF") {
    pdf(paste0(normalizePath(outputfolder), "/Forestplot_", 
               unique(edat[[genecol]]), "_", jobname, ".pdf"), 
        width = 4, height = 5)
    plot(gg)
    dev.off()
  }
  else if (draw == "HTML") {
    htmlwidgets::saveWidget(as_widget(ggplotly(gg)), paste0(normalizePath(outputfolder), 
                                                            "/Forestplot_", unique(edat[[genecol]]), jobname, 
                                                            ".html"))
  }
  return(gg)
}
