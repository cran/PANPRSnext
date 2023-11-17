PlinkLD_transform <- function(plinkLD, all_keep_snps) {
    wkeep <- which(plinkLD$SNP_A %in% all_keep_snps & plinkLD$SNP_B %in% all_keep_snps)
    plinkLD <- plinkLD[wkeep, ]

    if (length(which(is.na(plinkLD))) > 0) {
        stop("")
    }

    ld_J <- plinkLD[, c("SNP_B", "SNP_A", "R")]
    names(ld_J) <- c("SNP_A", "SNP_B", "R")

    ld_J <- rbind(plinkLD[, c("SNP_A", "SNP_B", "R")], ld_J)

    return(ld_J)
}
