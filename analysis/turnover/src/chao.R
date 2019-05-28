

chao <- function (data, datatype = c("abundance"),
                  k = 10, conf = 0.95){
    if (is.matrix(data) == T || is.data.frame(data) == T) {
        if (datatype != "incidence_raw") {
            if (ncol(data) == 1) {
                data <- data[, 1]
            }
            else {
                data <- data[1, ]
            }
        }
        else {
            t <- ncol(data)
            dat <- rowSums(data)
            dat <- as.integer(dat)
            t_infreq <- sum(colSums(data[which(dat < k), ]) >=
                1)
            data <- dat
            data <- c(t_infreq, t, data)
        }
    }
    if (datatype == "abundance_freq_count") {
        data <- as.integer(data)
        length4b <- length(data)
        data <- rep(data[seq(1, length4b, 2)], data[seq(2, length4b,
            2)])
        names(data) <- paste("x", 1:length(data), sep = "")
        datatype <- "abundance"
    }
    if (datatype == "incidence_freq_count") {
        t <- as.integer(data[1])
        data <- data[-c(1)]
        data <- as.integer(data)
        lengthdat <- length(data)
        data <- rep(data[seq(1, lengthdat, 2)], data[seq(2, lengthdat,
            2)])
        data <- c(t, data)
        names(data) <- c("T", paste("y", 1:(length(data) - 1),
            sep = ""))
        datatype <- "incidence_freq"
    }
    method <- "all"
    if (k != round(k) || k < 0)
        stop("Error: The cutoff t to define less abundant species must be non-negative integer!")
    if (is.numeric(conf) == FALSE || conf > 1 || conf < 0)
        stop("Error: confidence level must be a numerical value between 0 and 1, e.g. 0.95")
    if (datatype == "abundance") {
        f <- function(i, data) {
            length(data[which(data == i)])
        }
        if (f(1, data) == sum(data)) {
            stop("Error: The information of data is not enough.")
        }
        z <- (list(Basic_data_information = basicAbun(data, k)[[1]],
            Rare_species_group = RareSpeciesGroup(data, k), Species_table = round(SpecAbunOut(data,
                method, k, conf), 3)))
    }
    else if (datatype == "incidence_raw") {
        dat <- data[-1]
        Q <- function(i, data) {
            length(data[which(data == i)])
        }
        if (Q(1, dat) == sum(dat)) {
            stop("Error: The information of data is not enough.")
        }
        z <- (list(Basic_data_information = basicInci(data[-1],
            k)[[1]], Infreq_species_group = InfreqSpeciesGroup(data[-1],
            k), Species_table = round(SpecInciOut_raw(data, method,
            k, conf), 3)))
    }
    else if (datatype == "incidence_freq") {
        dat <- data[-1]
        Q <- function(i, data) {
            length(data[which(data == i)])
        }
        if (Q(1, dat) == sum(dat)) {
            stop("Error: The information of data is not enough.")
        }
        z <- (list(Basic_data_information = basicInci(data, k)[[1]],
            Infreq_species_group = InfreqSpeciesGroup(data, k),
            Species_table = round(SpecInciOut(data, method, k,
                conf), 3)))
    }
    else {
        stop("Error: The data type is wrong.")
    }
    class(z) <- c("ChaoSpecies")
    z
}
