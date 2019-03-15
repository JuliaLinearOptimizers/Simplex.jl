export simplexinv

function simplexinv(c, A, b, 𝔹=0, invB=0; max_iter = 20000)
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
    invB = sparse(Diagonal(signb))
    xB = abs.(b) # solution in current basis
  else
    artificial = false
    if (invB == 0) invB = sparse(inv(Matrix(A[:,𝔹]))) end
    ℕ = setdiff(1:n, 𝔹)
    xB = invB*b
  end
  getλ!(λ,c,𝔹)
  λ = invB'*λ # FIX
  q = getq(c, λ, A, ℕ)

  status = :Optimal

  while !(q == nothing || iter >= max_iter) # relative variable changes to directioner >= max_iter)
    iter += 1
    getAcol!(d,A,ℕ[q])
    d = invB * d # viable direction

    xq = Inf
    for k in 1:m # find min xB/d s.t. d .> 0
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

    subdot!(xB,d,xq) # update solution
    xB[p] = xq
    𝔹[p], ℕ[q] = ℕ[q], 𝔹[p] # update indexes
    #update of inverse of B
    E = one(zeros(m,m))
    dp = d[p]
    d[p] = -1
    E[:, p] = -d / dp
    invB = E*invB # STOP THIS
    getλ!(λ,c,𝔹)
    λ = invB'*λ # FIX
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
    if dot(xB, c[𝔹]) > eps(Float64)
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
          d = invB * d
          (abs(d[PivotAp]) >= eps(Float64)) ? break : q += 1
        end
        if q > length(ℕ)
          deleteat!(Irows, findfirst(A[Irows,p] .!= 0))
          deleteat!(𝔹, pind); deleteat!(xB, pind)
          deleteat!(Ap, pind); deleteat!(d, pind)
          invB = sparse(inv(Matrix(A[Irows,𝔹])))
        else
          𝔹[p] = ℕ[q]
          deleteat!(ℕ, q)
          d = invB * A[Irows,ℕ[q]]
          E = one(zeros(m,m))
          dp = d[p]
          d[p] = -1
          E[:, p] = -d / dp
          invB = E*invB # STOP THIS
        end
        p, pind = findmax(𝔹)
      end
      x, z, status = simplexinv(co, A[Irows,1:n], b[Irows], 𝔹, invB)
    else
      x, z, status = simplexinv(co, A[:,1:n], b, 𝔹, invB)
    end
  end
  return x, z, status
end
