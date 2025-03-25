function ProgressShower(T)
    p = Progress(100)
    stat(t, _) = update!(p, floor(Int, 100*t/T))
end

function TimeFilter(callbacks...; times)
    next::Int = 1
    function stats(t, s)
        while next <= lastindex(times) && times[next] <= t
            for cb in callbacks
                cb(times[next], s)
            end
            next += 1
        end
    end
end

function hexagon(x, y, r)
    polygon([[(x[i]+r*cos(π/3*(k+3/2)),
               y[i]+r*sin(π/3*(k+3/2))) for k in 0:6] 
                    for i in eachindex(x,y)])
end

function plot(M::Model, s::State)
    L = floor(Int, sqrt(length(M)))
    color = RGB.(s.nA ./ 30, s.nB ./ 30, 0)
    (x0,x1),(y0,y1) = extrema(M.posx), extrema(M.posy)
    compose(
        context(units=UnitBox(x0-2, y0-2, x1-x0+4, y1-y0+4)), 
        hexagon(M.posx, M.posy, 1),
        fill(color))
end

function Plotter(M::Model)
    stats(_, s) = display(plot(M, s))
end