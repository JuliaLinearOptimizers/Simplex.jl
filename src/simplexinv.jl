export simplexinv

function simplexinv(c, A, b, 𝔹=0, E = 0, ups = 0; max_iter = 20000)
  m, n = size(A)
  iter = 0
  λ = Array{Float64,1}(undef, m); d = Array{Float64,1}(undef, m)
  if 𝔹 == 0 # construct artificial problem
    artificial = true
    signb = sign.(b)
    A = [A sparse(Diagonal(signb))]
    𝔹 = collect(n+1:n+m) # indexes of basic variables
    ℕ = collect(1:n)
    co = copy(c)
    c = [zeros(n); ones(m)]
    E = zeros(A.m + 1, 3*A.m)
    for (i, j) in enumerate(b) # adjust for negative values in b
      if j < 0
        ups += 1
        E[i, ups] = -1; E[end, ups] = i
      end
    end
    xB = abs.(b) # solution in current basis
  else
    artificial = false
    if E == 0
      E = zeros(A.m + 1, 3*A.m)
      getPFI!(A, 𝔹, E); ups = A.m
    end
    ℕ = setdiff(1:n, 𝔹)
    xB = copy(b); solvePFI!(xB, E, ups)
  end
  getλ!(λ,c,𝔹)
  solvePFI!(λ, E, ups, true)
  q = getq(c, λ, A, ℕ)

  status = :Optimal

  while !(q == nothing || iter >= max_iter) # relative variable changes to direction
    iter += 1
    getAcol!(d,A,ℕ[q])
    solvePFI!(d, E, ups)

    xq = Inf
    for k in 1:m # find min xB/d such that d .> 0
      if d[k] >= eps(Float64)
        dfrac = xB[k]/d[k]
        if dfrac < xq
          xq = dfrac
          p = k
        end
      end
    end
    if xq == Inf
      status = :Unbounded; break
    end

    subdot!(xB, d, xq) # update solution
    xB[p] = xq
    𝔹[p], ℕ[q] = ℕ[q], 𝔹[p] # update indexes
    # update of inverse of B
    if ups < 3*A.m
      updatePFI!(E, d, p, ups); ups += 1
    else
      getPFI!(A, 𝔹, E); ups = A.m
      xB .= b; solvePFI!(xB, E, ups)
    end
    getλ!(λ, c, 𝔹)
    solvePFI!(λ, E, ups, true)
    q = getq(c, λ, A, ℕ)
  end

  if iter >= max_iter
    status = :UserLimit
  end

  x = zeros(n)
  if !artificial
    x[𝔹] = xB
    z = dot(c, x)
  else
    if dot(xB, c[𝔹]) > 1e-6
      status = (iter >= max_iter) ? :UserLimit : :Infeasible
      I = findall(𝔹 .<= n - m)
      x[𝔹[I]] = xB[I]
      z = dot(c[𝔹], x)
    elseif maximum(𝔹) > n # check for artificial variables in basis
      ℕ = setdiff(ℕ,n+1:n+m)
      Irows = collect(1:m)
      p, pind = findmax(𝔹)
      Ap = Array{Float64, 1}(undef, m)
      while p > n
        q = 1
        getAcol!(Ap,A,p,Irows)
        PivotAp = findfirst(Ap .!= 0) #findfirst(Ap)
        while q <= length(ℕ)
          getAcol!(d,A,ℕ[q],Irows)
          solvePFI!(d, E, ups)
          (abs(d[PivotAp]) >= eps(Float64)) ? break : q += 1
        end
        if q > length(ℕ)
          deleteat!(Irows, findfirst(A[Irows,p] .!= 0))
          deleteat!(𝔹, pind); deleteat!(xB, pind)
          deleteat!(Ap, pind); deleteat!(d, pind)
          E = zeros(length(Irows) + 1, 3*length(Irows))
          getPFI!(A[Irows, :], 𝔹, E); ups = length(Irows)
          xB .= b[Irows]; solvePFI!(xB, E, ups)
        else
          d = A[Irows, ℕ[q]]
          solvePFI!(d, E, ups)
          updatePFI!(E, d, p, ups); ups += 1
          if ups < 3*length(Irows)
            updatePFI!(E, d, p, ups); ups += 1
          else
            getPFI!(A[Irows, :], 𝔹, E); ups = length(Irows)
          end
          𝔹[p] = ℕ[q]
          deleteat!(ℕ, q)
        end
        p, pind = findmax(𝔹)
      end
      x, z, status = simplexinv(co, A[Irows,1:n], b[Irows], 𝔹, E, ups)
    else
      x, z, status = simplexinv(co, A[:,1:n], b, 𝔹, E, ups)
    end
  end
  return x, z, status
end
