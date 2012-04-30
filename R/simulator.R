RunSimulation <- function(initial.state, simulation, controller, 
                          visitors=NULL)
{
    # Runs 'simulation', starting at 'initial.state', in chunks until 
    # 'controller' indicates that the simulation should end.
    # 'visitors' are shown the state of the simulation at each chunk.
    # The return value is given by the controller.

    # simulation - an object that implements the S3 method RunChunk.
    # controller - an object that implements the S3 method ControlOnChunk.
    # visitors - a list of objects that implement the S3 methods VisitStart, 
    #            VisitChunk and VisitEnd.

    time <- 0
    current.state <- initial.state

    lapply(visitors, VisitStart, initial.state=initial.state)

    while(TRUE)
    {
        chunk <- RunChunk(simulation, time, current.state)

        # First row of this chunk is the last row of the previous chunk, 
        # or, if this is the first chunk, the same as initial.state.
        time <- tail(chunk, 1)[,1]
        current.state <- tail(chunk, 1)[,-1]

        lapply(visitors, VisitChunk, chunk=chunk)

        controller.result <- ControlOnChunk(controller, chunk)
        if(controller.result$terminate)
        {
            break
        }
    }

    lapply(visitors, VisitEnd, time=time, 
           final.state=controller.result$final.state)

    return (controller.result)
}

