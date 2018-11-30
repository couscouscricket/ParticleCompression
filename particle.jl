mutable struct Particle
    rx::Float64
    ry::Float64
    vx::Float64
    vy::Float64
    D::Float64 # growth rate D << 1
    radius::Float64
    mass::Float64
end

function move!(p::Particle, dt::Float64)
    p.rx += p.vx * dt
    p.ry += p.vy * dt
    p.radius += p.D * dt
end

function time2collide(p::Particle, q::Particle)
    # Particle never collides with itself
    if p === q; return Inf; end

    dx::Float64 = q.rx - p.rx
    dy::Float64 = q.ry - p.ry
    dvx::Float64 = q.vx - p.vx
    dvy::Float64 = q.vy - p.vy

    eps::Float64 = 1e-15

    D::Float64 = p.D + q.D
    s::Float64 = p.radius + q.radius
    dvdr::Float64 = dvx * dx + dvy * dy
    if (dvdr > -eps)
        return Inf
    end

    dvdv::Float64 = dvx^2 + dvy^2
    a::Float64 = dvdv - D^2
    if (-eps < a < eps)
        return Inf
    end

    drdr::Float64 = dx^2 + dy^2
    b::Float64 = dvdr - D * s
    d::Float64 = b^2 - a * (drdr - s^2)
    if (d < eps)
        return Inf
    end

    dt::Float64 = -(b + sqrt(d)) / a
    return dt < 0.0 ? 0.0 : dt
end

function time2hitVerticalWall(p::Particle)
    if p.vx > 0.0
        return (1.0 - p.radius - p.rx) / (p.vx + p.D)
    elseif p.vx < 0.0
        return (p.radius - p.rx) / (p.vx - p.D)
    else
        return Inf
    end
end
function time2hitHorizontalWall(p::Particle)
    if p.vy > 0.0
        return (1.0 - p.radius - p.ry) / (p.vy + p.D)
    elseif p.vy < 0.0
        return (p.radius - p.ry) / (p.vy - p.D)
    else
        return Inf
    end
end

function time2escapeCell(p::Particle, xmin::Float64, xmax::Float64,
        ymin::Float64, ymax::Float64)
    dtx = p.vx > 0.0 ? (xmax-p.rx)/p.vx : (xmin-p.rx)/p.vx
    dty = p.vy > 0.0 ? (ymax-p.ry)/p.vy : (ymin-p.ry)/p.vy
    return min(dtx,dty)
end

function collide!(p::Particle, q::Particle)
    dx = q.rx - p.rx
    dy = q.ry - p.ry

    dvx = q.vx - p.vx
    dvy = q.vy - p.vy
    drdr = dx*dx + dy*dy
    dvdr = dvx*dx + dvy*dy

    c = 2*dvdr/(drdr*(p.mass + q.mass))
    c1 = q.mass * c
    c2 = p.mass * c

    p.vx += c1 * dx
    p.vy += c1 * dy
    q.vx -= c2 * dx
    q.vy -= c2 * dy
end

collideVerticalWall!(p::Particle) = p.vx *= -1
collideHorizontalWall!(p::Particle) = p.vy *= -1
