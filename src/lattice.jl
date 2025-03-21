using LinearAlgebra, SparseArrays, Graphs, Compose, Colors

Open(n) = spdiagm(1=>trues(n-1))
Closed(n) = spdiagm(1=>trues(n-1), -n+1 => trues(1))
chain(n; boundary) = boundary(n) + boundary(n)'
rect(L1, L2; boundary=Closed) = kron(chain(L1; boundary),I(L2)) + kron(I(L1), chain(L2; boundary))
hexa(L1, L2; boundary=Closed) = rect(L1,L2; boundary) + kron(boundary(L1),boundary(L2)) + kron(boundary(L1)',boundary(L2)')
gen_square_lattice(L) = SimpleGraph(rect(L, L)), mod1.(1:L^2, L), fld1.(1:L^2, L)
function gen_hex_lattice(L)
    g = SimpleGraph(hexa(L, L))
    x, y = mod1.(1:L^2, L), fld1.(1:L^2, L)
    g, mod1.((y .- 1)/2 + x .- 1, L), (y .- 1) .* sqrt(3) ./ 2 
end


