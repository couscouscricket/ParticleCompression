mutable struct Box
    events::Vector{Event}
    t::Float64

    particles::Vector{Particle}
    collisionCount::Vector{Int}
    positionInGrid::Vector{Vector{Int}}

    grid::Matrix{Vector{Int}}
    nx::Int
    lx::Float64

    function Box(nParticles::Int, growthRate::Vector{Float64})
        # Random particles inside a unit box.
        particles::Vector{Particle} = Vector{Particle}(undef, nParticles)
        maxGrowth::Float64 = maximum(growthRate)
        for i = 1:nParticles
            # Random positiion
            rx::Float64 = rand()
            ry::Float64 = rand()

            # Random velocity of magnitude 1
            vx::Float64 = -1 + 2*rand()
            vy::Float64 = -1 + 2*rand()
            norm::Float64 = sqrt(vx^2 + vy^2)
            vx = vx/norm
            vy = vy/norm

            # Random growth rate
            D::Float64 = rand(growthRate)
            mass::Float64 = D/maxGrowth

            particles[i] = Particle(rx, ry, vx, vy, D, 0.0, mass)
        end

        width::Float64 = 1.0
        height::Float64 = 1.0

        # Create space partitioning data structure.
        numSpecies::Int = length(growthRate)
        nx::Int = floor(sqrt(nParticles/numSpecies))
        lx::Float64 = width / nx

        # Each bucket will have an integer position (x,y) mapped
        # to a linear array representing the grid.
        # There will be nx*nx cells in the grid.
        n::Int = nx*nx

        # Position particles on grid and init collision counter
        grid = [Int[] for row in 1:nx, col in 1:nx]
        tc = Vector{Float64}(undef, nParticles)
        collisionCount = Vector{Int}(undef, nParticles)
        positionInGrid = Vector{Vector{Int}}(undef, nParticles)
        for i in eachindex(particles)
            p = particles[i]
            col::Int = cld(p.rx, lx)
            row::Int = cld(p.ry, lx)
            tc[i] = 0.0
            collisionCount[i] = 0
            positionInGrid[i] = grid[row, col]
            push!(grid[row, col], i)
        end

        box = new(Event[], 0.0, particles, collisionCount, positionInGrid,
                  grid, nx, lx)

        for i in eachindex(box.particles)
            predict(i, box)
        end

        return box
    end
end

function predict(i::Int, box::Box)
    t = box.t
    # Current particle, its collision count
    # and its position in the grid.
    p::Particle = box.particles[i]
    pcc::Int = box.collisionCount[i]
    col::Int = cld(p.rx, box.lx)
    row::Int = cld(p.ry, box.lx)

    # Calculate the time it takes for the particle
    # to escape its bucket.
    eps = 1e-6
    xmin = (col - 1) * box.lx - eps
    xmax = col * box.lx + eps
    ymin = (row - 1) * box.lx - eps
    ymax = row * box.lx + eps
    escapeTime = t + time2escapeCell(p, xmin, xmax, ymin, ymax)

    # Queue as grid event.
    queue!(box.events, Event(escapeTime, i, i, pcc, pcc))

    # If particle is near wall then check for wall collisions.
    if row == 1 || row == box.nx
        tt = t + time2hitHorizontalWall(p)
        if tt < escapeTime
            queue!(box.events, Event(tt, 0, i, -1, pcc))
        end
    end
    if col == 1 || col == box.nx
        tt = t + time2hitVerticalWall(p)
        if tt < escapeTime
            queue!(box.events, Event(tt, i, 0, pcc, -1))
        end
    end

    # Check collisions with neighbouring particles.
    # and create every particle-particle event.
    for n in col-1:col+1, m in row-1:row+1
        if 1 <= m <= box.nx && 1 <= n <= box.nx
            bucket::Vector{Int} = box.grid[m, n]
            for idx in bucket
                q::Particle = box.particles[idx]
                if i == idx; continue; end
                tt = t + time2collide(p, q)
                if tt < escapeTime
                    qcc = box.collisionCount[idx]
                    queue!(box.events, Event(tt, i, idx, pcc, qcc))
                end
            end
        end
    end
end

function resolveCollisionEvent(event::Event, box::Box)
    i = event.i
    j = event.j

    dt = event.t - box.t
    for p in box.particles
        move!(p, dt)
    end
    box.t = event.t

    if i > 0 && j > 0
        p = box.particles[i]
        q = box.particles[j]
        collide!(p, q)
        box.collisionCount[i] += 1
        box.collisionCount[j] += 1
        predict(i, box)
        predict(j, box)
    elseif i > 0
        p = box.particles[i]
        collideVerticalWall!(p)
        box.collisionCount[i] += 1
        predict(i, box)
    elseif j > 0
        p = box.particles[j]
        collideHorizontalWall!(p)
        box.collisionCount[j] += 1
        predict(j, box)
    end
end

function resolveGridEvent(event::Event, box::Box)
    # Particle in question is either pi or pj since,
    # by construction, i = j on grid events.
    i = event.i
    p = box.particles[i]

    # First, locate particle's bucket.
    bucket::Vector{Int} = box.positionInGrid[i]

    # Remove particle's reference from current bucket.
    n::Int = findfirst(x -> x == i, bucket)
    deleteat!(bucket, n)

    # Advance particles in time.
    dt = event.t - box.t
    for q in box.particles
        move!(q, dt)
    end
    box.t = event.t

    # Update particle's position inside grid.
    col::Int = cld(p.rx, box.lx)
    row::Int = cld(p.ry, box.lx)
    push!(box.grid[row, col], i)
    box.positionInGrid[i] = box.grid[row, col]

    # Update event queue by predicting new
    # collisions on current particle.
    predict(i, box)
end


function update!(box::Box)
    # Discard invalid events and resolve grid events beforehand.
    # Grid events move the particles along the grid without collisions.
    local event
    while true
        event = dequeue!(box.events)
        if isValid(event, box.collisionCount)
            if event.i > 0 && event.i == event.j;
                resolveGridEvent(event, box)
            else
                break
            end
        end
    end

    # Current event is a particle-wall or particle-particle collision event.
    # Move the particles in time until their time is event's time and then
    # resolve the collision event.
    resolveCollisionEvent(event, box)

    #queue!(box.events, Event(t + 10, 0, 0 , 0, 0))
end
