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