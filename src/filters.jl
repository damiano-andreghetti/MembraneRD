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

function plot(M::Model, s::State)
    L = floor(Int, sqrt(length(M)))
    color = RGB.(s.nA ./ 30, s.nB ./ 30, 0)
    x,y = M.posx./maximum(abs, M.posx), M.posy/maximum(abs, M.posy)
    compose(context(units=UnitBox(-1.2, -1.4, 2.4, 2.8)), ngon(x, y, fill(2/L/sqrt(3), length(M.posx)), fill(6, length(M.posx))), fill(color))
end

function Plotter(M::Model)
    stats(_, s) = display(plot(M, s))
end