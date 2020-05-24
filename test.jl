include("./particle.jl")
include("./event.jl")
include("./priorityQueue.jl")
include("./box.jl")

# Density inside a unit box.
function getPackingFraction(box::Box)
    areaParticles::Float64 = 0
    for p in box.particles
        areaParticles += π * p.radius^2
    end
    return areaParticles
end

function plotBox(box::Box, growthRate)
    phi::Float64 = getPackingFraction(box)
    gp = open(`gnuplot -p`, "w")
    println(gp, "
            set terminal qt font 'Helvetica,11';
            set colorsequence podo; unset colorbox;
            set style fill solid 0.5 border 1;
            set grid back lc 1 lw 2;
            set xtics $(box.lx);
            set ytics $(box.lx);
            set size ratio -1; set tics out nomirror;
            set xrange [0:1]; set yrange [0:1];
            ")
    println(gp, "set title 'φ = $(getPackingFraction(box)),   growthRate = $(growthRate)'")
    println(gp, "plot '-' w circles lc palette notitle;")
    for p in box.particles
        r = p.radius
        x = p.rx
        y = p.ry
        println(gp, x, " ", y, " ", r, " ", r)
    end
    println(gp, "e")
    close(gp)
end

function test()
    nParticles = 100
    growthRate = [8e-3, 6e-3]
    box = Box(nParticles, growthRate)

    @time for i in 1:100000
        update!(box)
    end

    plotBox(box, growthRate)
end

test()
