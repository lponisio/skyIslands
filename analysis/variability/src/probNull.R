

calcProbNull <- function(web, N){
    probNull <- function(M) {
        randiag <- function(vdims) {
            mats <- diag(vdims)
            ## what cell number are the interactions
            posdiag <- 1:vdims + ((1:vdims - 1) * vdims)
            desp <- matrix(rep(1:vdims, vdims), vdims, vdims,
                           byrow = TRUE) - (1:vdims)
            fdesp <- sample(1:vdims)
            move <- desp[cbind(1:vdims, fdesp)]
            moved <- posdiag + move
            mdesp <- matrix(0, vdims, vdims)
            mdesp[moved] <- 1
            return(mdesp)
        }
        M <- as.matrix(M)
        values <- M[which(M > 0)]
        lvalues <- length(values)
        if (identical(dim(M)[1], dim(M)[2])) {
            vdims <- dim(M)[1]
            vdimb <- dim(M)[1]
        }
        if (!identical(dim(M)[1], dim(M)[2])) {
            dims <- which(dim(M) == min(dim(M)))
            vdims <- dim(M)[dims]
            dimb <- which(dim(M) == max(dim(M)))
            vdimb <- dim(M)[dimb]
        }
        MR <- matrix(0, vdims, vdimb)
        lMR <- vdims * vdimb
        sample1 <- sample(vdimb, vdims)
        diag1 <- randiag(vdims)
        MR[, sample1] <- diag1
        sample2 <- (1:vdimb)[-sample1]
        pos <- sample(vdims, length(sample2), replace = TRUE)
        MR[cbind(pos, sample2)] <- 2
        MRoccupied <- which(MR > 0)
        vleft <- lvalues - vdimb
        if (vleft > 0)
            MRoccupied <- c(MRoccupied, sample((1:lMR)[-MRoccupied],
                                               vleft))
        MR[MRoccupied] <- sample(values)
        if (dim(MR)[1] != dim(M)[1])
            MR <- t(MR)
        return(MR)
    }
    replicate(N, probNull(web), simplify = FALSE)
}
