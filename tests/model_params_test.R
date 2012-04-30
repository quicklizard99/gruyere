TestModelParamsSpec <- function()
{
    # TODO
}

TestLoadModelParamsSpec <- function(dir)
{
    # TODO
}

TestIntermediateModelParams <- function()
{
    # We should get the same smallest primary producer, s, regardless of 
    # species order
    c1 <- Community(nodes=data.frame(node=paste('Species', 1:5), 
                                     category=c(rep('producer', 3), 
                                                rep('invertebrate', 2)),
                                     M=seq(10e-6, 10e-3, length.out=5), 
                                     N=10:6), 
                    properties=list(title='c1', M.units='kg', N.units='m^-2'))
    c2 <- Community(nodes=data.frame(node=paste('Species', 5:1), 
                                     category=c(rep('invertebrate', 2), 
                                                rep('producer', 3)),
                                     M=seq(10e-3, 10e-6, length.out=5), 
                                     N=6:10), 
                    properties=list(title='c2', M.units='kg', N.units='m^-2'))

    s1 <- IntermediateModelParams(c1, ModelParamsSpec())[['s']]
    s2 <- IntermediateModelParams(c2, ModelParamsSpec())[['s']]

    stopifnot(1==s1)
    stopifnot(5==s2)
    stopifnot(names(s1)==names(s2))
}

TestBuildModelParams <- function()
{
}

