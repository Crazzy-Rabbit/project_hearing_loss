#------------- QC for GWAS summary ----//
sums_col_treate <- function(gwas, SNP, A1, A2, freq, P,
                            CHR = NULL, POS = NULL,
                            BETA = NULL, SE = NULL,
                            OR = NULL, Z = NULL, N=NULL) {
    gwas_col = c(
        CHR   =   CHR,
        POS   =   POS, 
        SNP   =   SNP,
        A1    =    A1, 
        A2    =    A2,
        freq  =  freq,
        b     =  BETA,
        se    =    SE,
        OR    =    OR,
        Z     =     Z,
        p     =     P,
        N     =     N)
    gwas_col = gwas_col[!sapply(gwas_col, is.null)]

    missing_required = setdiff(c("SNP", "A1", "A2", "freq", "p"), names(gwas_col))
    if (length(missing_required) > 0) {
        stop("Missing required columns: ", paste(missing_required, collapse = ", "))
    }
    has_b_se = all(c("b", "se") %in% names(gwas_col))
    has_or   = "OR" %in% names(gwas_col)
    has_z    = "Z" %in% names(gwas_col)
    if (!(has_b_se || has_or || has_z)) {
        stop("Must provide at least one of the following: (BETA + SE), OR, or Z.")
    }

    gwas = gwas[, unname(gwas_col), with=FALSE]
    setnames(gwas, old = unname(gwas_col), new = names(gwas_col))
    
    valid_p <- !is.na(gwas$p) & gwas$p >= 0 & gwas$p <= 1
    if (has_b_se){    
        gwas = gwas[valid_p & !is.na(b) & !is.na(se) & b != 0, ]
    } else if (has_or){
        gwas = gwas[valid_p & !is.na(OR), ]
    } else if (has_z){
        gwas = gwas[valid_p & !is.na(Z), ]
    }
    
    return(gwas)
}

library(data.table)
gwas = fread("Ever_smoke.fastGWA.gz")
clean = sums_col_treate(gwas, "variant_id", "effect_allele", "other_allele", "effect_allele_frequency", "p_value", BETA="beta", SE="standard_error" N="N")
clean = sums_col_treate(gwas, "SNP", "A1", "A2", "AF1", "P", BETA="BETA", SE="SE", N="N")
clean = clean[, .(SNP, A1, A2, freq, b, se, p, N)]
fwrite(clean, "../cleanGWAS/clean_ever_smoke.fastGWA.gz", sep="\t")