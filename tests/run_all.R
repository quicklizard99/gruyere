#!/usr/bin/env Rscript
options(warn=2)
library(gruyere)

RunTests <- function(tests)
{
    # tests should be a vector of function names
    if(0==length(tests))
    {
        stop('No tests to run! Is the working directory gruyere/tests ?')
    }
    else
    {
        failed <- NULL

        for(test in tests)
        {
            cat(paste('Running [', test, ']\n', sep=''))
            res <- tryCatch(do.call(test, args=list()), error=function(e) e)
            if(!is.null(res))
            {
                # Strip whitespace from error message
                res <- gsub('\n', '', res)
                cat(paste('[', test, '] raised unexpected error [', res, ']\n', 
                    sep=''))
                traceback()
                failed <- c(failed, test)
            }
            else
            {
                
            }
        }

        cat(paste(length(tests), 'tests ran.\n'))
        cat(paste(length(tests) - length(failed), 'passed.\n'))
        if(!is.null(failed))
        {
            cat(paste(length(failed), 'failed.\n'))
            stop()
        }
    }
}

# Source all files in this dir except this one
files <- list.files(getwd(), pattern='*R$')
files <- setdiff(files, 'run_all.R')
junk <- sapply(file.path(getwd(), files), source)
RunTests(ls(pattern=glob2rx('Test*')))

