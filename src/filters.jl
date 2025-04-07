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
    posx, posy = mm .* M.posx, mm .* M.posy
    (x0,x1),(y0,y1) = extrema(posx), extrema(posy)
    set_default_graphic_size(x1-x0+3mm, y1-y0+3mm)
    compose(context(), hexagon(posx, posy, 1mm), fill(RGB.(0.0, s.nB ./ 30, s.nA ./ 30))    )
end

function Plotter(M::Model)
    stats(_, s) = display(plot(M, s))
end