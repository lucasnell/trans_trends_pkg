
#
# Copyright (C)       2011-2019 Douglas Bates, Martin Maechler, Ben Bolker,
#                               Steven Walker, Rune Haubo Bojesen Christensen,
#                               Henrik Singmann, Bin Dai, Fabian Scheipl,
#                               Gabor Grothendieck, Peter Green, and John Fox
#
# The code below is part of the lme4 package (version 1.1-21).
# No significant changes have been made since its incorporation into this package on
# 4 Dec 2019.
#




# -------------------*
# lme4::nobars
# -------------------*


anyBars <- function (term) {
    any(c("|", "||") %in% all.names(term))
}


isBar <- function (term) {
    if (is.call(term)) {
        if ((term[[1]] == as.name("|")) || (term[[1]] == as.name("||"))) {
            return(TRUE)
        }
    }
    FALSE
}

isAnyArgBar <- function (term) {
    if ((term[[1]] != as.name("~")) && (term[[1]] != as.name("("))) {
        for (i in seq_along(term)) {
            if (isBar(term[[i]]))
                return(TRUE)
        }
    }
    FALSE
}


nobars_ <- function (term) {
    if (!anyBars(term))
        return(term)
    if (isBar(term))
        return(NULL)
    if (isAnyArgBar(term))
        return(NULL)
    if (length(term) == 2) {
        nb <- nobars_(term[[2]])
        if (is.null(nb))
            return(NULL)
        term[[2]] <- nb
        return(term)
    }
    nb2 <- nobars_(term[[2]])
    nb3 <- nobars_(term[[3]])
    if (is.null(nb2))
        return(nb3)
    if (is.null(nb3))
        return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
}


nobars <- function (term) {
    nb <- nobars_(term)
    if (is(term, "formula") && length(term) == 3 && is.symbol(nb)) {
        nb <- reformulate("1", response = deparse(nb))
    }
    if (is.null(nb)) {
        nb <- if (is(term, "formula")) {
            ~1
        } else 1
    }
    nb
}





# -------------------*
# lme4::findbars
# -------------------*



expandDoubleVerts <- function (term) {
    expandDoubleVert <- function(term) {
        frml <- formula(substitute(~x, list(x = term[[2]])))
        newtrms <- paste0("0+", attr(terms(frml), "term.labels"))
        if (attr(terms(frml), "intercept") != 0) newtrms <- c("1", newtrms)
        as.formula(paste("~(", paste(vapply(newtrms,
                                            function(trm) paste0(trm, "|",
                                                                 deparse(term[[3]])), ""),
                                     collapse = ")+("),
                         ")"))[[2]]
    }
    if (!is.name(term) && is.language(term)) {
        if (term[[1]] == as.name("(")) {
            term[[2]] <- expandDoubleVerts(term[[2]])
        }
        stopifnot(is.call(term))
        if (term[[1]] == as.name("||"))
            return(expandDoubleVert(term))
        term[[2]] <- expandDoubleVerts(term[[2]])
        if (length(term) != 2) {
            if (length(term) == 3)
                term[[3]] <- expandDoubleVerts(term[[3]])
        }
    }
    term
}


findbars <- function (term) {
    fb <- function(term) {
        if (is.name(term) || !is.language(term)) return(NULL)
        if (term[[1]] == as.name("(")) return(fb(term[[2]]))
        stopifnot(is.call(term))
        if (term[[1]] == as.name("|")) return(term)
        if (length(term) == 2) return(fb(term[[2]]))
        c(fb(term[[2]]), fb(term[[3]]))
    }
    expandSlash <- function(bb) {
        makeInteraction <- function(x) {
            if (length(x) < 2) return(x)
            trm1 <- makeInteraction(x[[1]])
            trm11 <- if (is.list(trm1))
                trm1[[1]]
            else trm1
            list(substitute(foo:bar, list(foo = x[[2]], bar = trm11)),
                 trm1)
        }
        slashTerms <- function(x) {
            if (!("/" %in% all.names(x))) return(x)
            if (x[[1]] != as.name("/")) stop("unparseable formula for grouping factor",
                                             call. = FALSE)
            list(slashTerms(x[[2]]), slashTerms(x[[3]]))
        }
        if (!is.list(bb)) {
            expandSlash(list(bb))
        } else unlist(lapply(bb, function(x) {
            if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]]))) {
                lapply(unlist(makeInteraction(trms)),
                       function(trm) substitute(foo | bar, list(foo = x[[2]], bar = trm)))
            } else x
        }))
    }
    modterm <- expandDoubleVerts(if (is(term, "formula")) term[[length(term)]] else term)
    expandSlash(fb(modterm))
}

