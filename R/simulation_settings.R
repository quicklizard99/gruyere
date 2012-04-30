SaveSimulationSettings <- function(simulation.settings, dir)
{
    write.csv(do.call('cbind.data.frame', simulation.settings), 
              file.path(dir, 'simulation.settings.csv'),
              row.names=FALSE)
}

LoadSimulationSettings <- function(dir)
{
    sim.settings <- read.csv(file.path(dir, 'simulation.settings.csv'))
    stopifnot(1==nrow(sim.settings))
    sim.settings <- do.call(list, sim.settings[1,,drop=FALSE])
    return (sim.settings)
}

