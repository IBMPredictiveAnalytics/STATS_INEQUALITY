#/***********************************************************************
# * Licensed Materials - Property of IBM 
# *
# * IBM SPSS Products: Statistics Common
# *
# * (C) Copyright IBM Corp. 2014
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp. 
# ************************************************************************/

# author__ = "SPSS, JKP"
# version__ = "1.0.1"

# History
# 25-apr-2014 Original Version


helptext="STATS INEQUALITY VARIABLE=variable name THEILP=variable name
MEASURE=GINI THEIL ATKINSON KOLM COEFVAR COEFVAR2 ENTROPY
ATKINSONPAR=value 
/OPTIONS GINICITYPE=PERC BCA GINICILEVEL=percentile
/OUTPUT LORENZPLOT=YES or NO.

Example:
STATS INEQUALITY VARIABLE=salary MEASURE=GINI ATKINSON
    /OPTIONS GINICILEVEL=95
    /OUTPUT LORENZPLOT=YES.

Compute various distribution inequality measures.
This command supports split files.

VARIABLE specifies the variable whose distribution is used.  Values
should be nonnegative for most measures.  Missing values are always
deleted.

MEASURE specifies one or more measures to compute.
GINI is the Gini index and is the default.
THEIL is the Theil U statistic.  It requires the THEILP keyword
ATKINSON is the Atkinson index.
KOLM is the Kolm measure.
COEFVAR is the coefficient of variation.
COEFVAR2 is the squared coefficient of variation.
ENTROPY is the entropy measure.

ATKINSONPAR is the parameter for the Atkinson index
and defaults to .5.  The value must be between 0 and 1.
Larger values increase the sensitivity of the index to changes
at the lower end of the distribution while, similarly, smaller
values make the index more sensitive to changes at the upper end.

GINICILEVEL specifies a confidence interval as a percent.  If
specified, the output includes lower and upper confidence intervals
for the Gini index.

GINICITYPE specifies how to compute the interval.  BCA specifies
bias-corrected bootstrap, and PERC specifies percentile boostrap.
BCA is the default.

LORENZPLOT specifies whether or not to plot the Lorenz curve for
the distribution.  The default is NO.

STATS INEQUALITY /HELP.  prints this information and does nothing else.
"

### MAIN ROUTINE ###
doineq = function(variable, theilp=NULL, measure=list("gini"), atkinsonpar=.5, ginicitype="bca",
    ginicilevel=NA, lorenzplot=FALSE) {
    # Calculate inequality measures and lorenz curve
    
    setuplocalization("STATS_INEQUALITY")
    
    # A warnings proc name is associated with the regular output
    # (and the same omsid), because warnings/errors may appear in
    # a separate procedure block following the regular output
    procname=gtxt("Inequality Measures")
    warningsprocname = gtxt("Inequality Measures: Warnings")
    omsid="STATSINEQUALITY"
    warns = Warn(procname=warningsprocname,omsid=omsid)

    # note Atkinson, Gini, and Lc functions exist in both modules
    tryCatch(library(DescTools), error=function(e){
        warns$warn(gtxtf("The R %s package is required but could not be loaded.", "DescTools"),dostop=TRUE)
        }
    )
    tryCatch(library(ineq), error=function(e){
        warns$warn(gtxtf("The R %s package is required but could not be loaded.", "ineq"),dostop=TRUE)
    }
    )
    if ("theilu" %in% measure && is.null(theilp)) {
        warns$warn(gtxt("The Theil measure requires a prediction variable"), dostop=TRUE)
    }
    ginicilevel = ginicilevel /100.  # rescaled for the R functions
    allargs = as.list(environment())

    # functions to use in a uniform manner via lapply using closures
    flist = list(
        gini = function(x) {
            DescTools::Gini(x, conf.level=ginicilevel,
                               type=ginicitype)
        },
        atkinson = function(x) {
            DescTools::Atkinson(x, parameter=atkinsonpar)
        },
        theilu = function(x) {
            DescTools::TheilU(x, dta[[2]], method=2)
        },
        # Kolm can take a parameter with any positive value
        # but not currently supported.
        kolm = function (x) {ineq::Kolm(x,  parameter=1)
        },
        coefvar = function(x) {ineq::var.coeff(x, square=FALSE)
        },
        coefvar2 = function(x) {ineq::var.coeff(x, square=TRUE)
        },
        entropy = function(x) {ineq::entropy(x)
        }
    )
    # reorder requested measures so output doesn't depend on input order
    # the list must match the MEASURE vallist list in Run
    knownmeasures = list("gini", "atkinson", "theilu", "kolm",
                         "coefvar", "coefvar2", "entropy")
    knownmeasures = Filter(function(m) m %in% measure, knownmeasures)
    measure = measure[match(measure, knownmeasures)]
    hdr = buildlabels(measure, ginicilevel, theilp, atkinsonpar)

    StartProcedure(procname, omsid)
    while (!spssdata.IsLastSplit()) {
        dta = spssdata.GetSplitDataFromSPSS(c(variable, theilp), missingValueToNA=TRUE)
        # with only 1 col, output needs to be coerced back to data frame :-()
        dta=data.frame(dta[complete.cases(dta),])
        res = computemeasures(measure, dta, warns)
        displayresults(allargs, res, hdr, warns)
        if (lorenzplot) {
            lcc = DescTools::Lc(dta[[1]], na.rm=TRUE, plot=FALSE)
            plot(lcc, main=gtxtf("Lorenz Curve of Variable %s", variable))
        }
    }
    spssdata.CloseDataConnection()
    warns$display(inproc=TRUE)
    ###spsspkg.EndProcedure()

}

buildlabels = function(meas, ginicilevel, theilp, atkinsonpar) {
    # return list of labels for pivot table
    ginicilevel = ginicilevel * 100  # express in user format
    labels = list(
        gini=if(!is.na(ginicilevel)) c(gtxt("Gini"), 
            gtxtf("Lower %s%% CI", ginicilevel), gtxtf("Upper %s%% CI", ginicilevel)) else
            gtxt("Gini"),
        atkinson = gtxtf("Atkinson(%s)", atkinsonpar),
        theilu = gtxtf("Theil U (%s)", theilp),
        kolm = gtxt("Kolm"),
        coefvar = gtxt("Coef. of Variation"),
        coefvar2 = gtxt("Squared Coef. of Variation"),
        entropy = gtxt("Entropy")
    )

    used=labels[meas]
    return(used)
}


computemeasures = function(measure, dta, warns) {
    # return a row of selected inequality measures

    # get wrapped function definitions from parent
    flist=get('flist', envir=parent.frame())
    # apply functions for the selected measures picking up
    # additional parameters from parent function definitions
    res = tryCatch({
       lapply(measure, function(f) flist[[f]](dta[[1]]))
    }, error=function(e) {warns$warn(e$message, dostop=TRUE)},
       warning=function(e) {warns$warn(e$message, dostop=FALSE)}
    )
    return(unlist(res))
}

displayresults = function(allargs, res, hdr, warns) {
    # Produce pivot tables and charts

    spsspivottable.Display(as.data.frame(rbind(res)), title=gtxtf("Inequality Measures for Variable %s",
        allargs[["variable"]]), 
       templateName="STATSINEQUALITYMEAS",
       hiderowdimlabel=TRUE,
       collabels= hdr,
       caption=gtxt("Results computed by R\nDescTools or ineq package"),
       isSplit=TRUE,
    )
}

Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = list2env(lcl) # makes this list into an environment

    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.

        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 

        if (is.null(msg) || dostop) {
            lcl$display(inproc)  # display messages and end procedure state
            if (dostop) {
                stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
            }
        }
    }
    
    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any

        if (lcl$msgnum == 0) {   # nothing to display
            if (inproc) {
                spsspkg.EndProcedure()
            }
        } else {
            if (!inproc) {
                procok =tryCatch({
                    StartProcedure(lcl$procname, lcl$omsid)
                    TRUE
                    },
                    error = function(e) {
                        FALSE
                    }
                )
            }
            if (procok) {  # build and display a Warnings table if we can
                table = spss.BasePivotTable("Warnings ","Warnings") # do not translate this
                rowdim = BasePivotTable.Append(table,Dimension.Place.row, 
                    gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)

                for (i in 1:lcl$msgnum) {
                    rowcategory = spss.CellText.String(as.character(i))
                    BasePivotTable.SetCategories(table,rowdim,rowcategory)
                    BasePivotTable.SetCellValue(table,rowcategory, 
                        spss.CellText.String(lcl$msglist[[i]]))
                }
                spsspkg.EndProcedure()   # implies display
            } else { # can't produce a table
                for (i in 1:lcl$msgnum) {
                    print(lcl$msglist[[i]])
                }
            }
        }
    }
    return(lcl)
}

# localization initialization
setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 
# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
        spsspkg.StartProcedure(procname, omsid)
    }
    else {
        spsspkg.StartProcedure(omsid)
    }
}

gtxt <- function(...) {
    return(gettext(...,domain="STATS_INEQUALITY"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_INEQUALITY"))
}


Run = function(args) {
    #Execute the STATS INEQUALITY command

    cmdname = args[[1]]
    args = args[[2]]
    oobj = spsspkg.Syntax(list(
        spsspkg.Template("VARIABLE", subc="",  ktype="existingvarlist", var="variable"),
        spsspkg.Template("THEILP", subc="", ktype="existingvarlist", var="theilp"),
        spsspkg.Template("MEASURE", subc="", ktype="str", var="measure", islist=TRUE,
            vallist=list("gini", "theilu", "atkinson", "kolm", 
            "coefvar", "coefvar2", "entropy")),
        spsspkg.Template("ATKINSONPAR", subc="", ktype="float", var="atkinsonpar",
            vallist=list(0,1)),
        
        spsspkg.Template("GINICITYPE", subc="OPTIONS", ktype="str", var="ginicitype",
            vallist=list("perc", "bca")),
        spsspkg.Template("GINICILEVEL", subc="OPTIONS", ktype="float", var="ginicilevel",
            vallist=list(50, 99.9999)),
        
        spsspkg.Template("LORENZPLOT", subc="OUTPUT", ktype="bool", var="lorenzplot")
    ))

    # A HELP subcommand overrides all else
    if ("HELP" %in% attr(args,"names")) {
        #writeLines(helptext)
        helper(cmdname)
    }
    else {
        res <- spsspkg.processcmd(oobj, args, "doineq")
    }
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}