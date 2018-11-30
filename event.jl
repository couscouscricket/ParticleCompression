struct Event
    t::Float64
    i::Int
    j::Int
    iCount::Int
    jCount::Int
end

function isValid(e::Event, collisionCount::Vector{Int})
    if e.i > 0 && e.iCount != collisionCount[e.i]
        return false
    end
    if e.j > 0 && e.jCount != collisionCount[e.j]
        return false
    end
    return true
end

import Base.<
function <(a::Event, b::Event)
    return a.t < b.t
end
