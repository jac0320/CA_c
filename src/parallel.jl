if config.PARALLEL   # The whole script runs only when in parallel mode

	info("Reporting from the head node, I am head worker")

    if length(config.WORKERS) > 1
        if config.JOBPERWORKER > 1
            WORKERS = []
            for i = 1:config.JOBPERWORKER
                WORKERS = [WORKERS;config.WORKERS]
            end
        else
            WORKERS  = config.WORKER
        end
        addprocs(WORKERS)
    else
        NUMWORKERS = config.WORKERS[1]
		addprocs(NUMWORKERS)
    end

    # Spwan all workers to the current dir
    function sendto(p::Int; args...)
        for (nm, val) in args
            @spawnat(p, eval(Main, Expr(:(=), nm, val)))
        end
    end

    for (idx, pid) in enumerate(workers())
        sendto(pid, SourceDir = SourceDir)
    end

    @everywhere cd(SourceDir)

    # Report workers amount
    totalWorkers = length(workers())

	for worker in workers()
        remotecall(include, worker, "src/main.jl")
    end

	sleep(20) #Hard waiting for other nodes to load all packages
end
