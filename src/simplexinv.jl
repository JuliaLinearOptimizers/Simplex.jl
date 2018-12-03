export simplexinv

function simplexinv(c, A, b, IB=0, invB=0; max_iter = 4000)
  m, n = size(A)
  iter = 0
  if IB == 0 # construct artificial problem
    artificial = true
    signb = sign.(b)
    A = [A diagm(signb)]
    Im = speye(m)
    IB = collect(n+1:n+m) # indexes of basic variables
    IN = collect(1:n)
    co = collect(c)
    c = [zeros(n); ones(m)]
    invB = collect(diagm(signb))
    r = -A[:,IN]'*signb # artificial relative costs
    xB = collect(abs.(b)) # solution in current basis
  else
    artificial = false
    Im = speye(m)
    IN = setdiff(1:n, IB)
    r = c[IN] - A[:,IN]' * (invB' * c[IB])
    xB = invB*b
  end
  q = findfirst(r .< 0) # Bland's Rule

  status = :Optimal

  while !(q == 0 || iter >= max_iter)
    iter += 1
    @assert all(xB .>= 0)
    d = invB * A[:,IN[q]] # viable direction
    apfrac = xB ./ d # relative variable changes to direction
    indpos = find(d .> 0) # variables that decrease in d direction
    if length(indpos) == 0
      status = :Unbounded
      break
    end
    indxq = indmin(apfrac[indpos])
    xq = apfrac[indpos[indxq]]
    @assert xq >= 0
    @assert xq < Inf

    p = findfirst(apfrac, xq) # Bland's Rule
    xB -= xq * d; xB[p] = xq # update solution
    IB[p], IN[q] = IN[q], IB[p] # update indexes
    #update of inverse of B
    Q = collect(Im)
    dp = d[p]
    d[p] = -1
    Q[:, p] = -d / dp
    invB = Q*invB #no need to create Q
    r = c[IN] - A[:,IN]' * (invB' * c[IB])
    q = findfirst(r .< 0) # Bland's Rule
  end

  if iter >= max_iter
    status = :UserLimit
  end

  x = zeros(n)
  if !artificial
    x[IB] = xB
    z = dot(c, x)
  else
    if dot(xB, c[IB]) > 0
      status = :Infeasible
      I = find(IB .<= n - m)
      x[I] = xB[I]
      z = dot(co, x)
    elseif maximum(IB) > n # check for artificial variables in basis
      deleteat!(IN, find(IN .> n))
      Irows = collect(1:m)
      p = findfirst(IB .> n)
      while p != 0
        q = findfirst(invB[p,:]' * A[Irows,IN] .!= 0)
        if q == 0
          deleteat!(Irows, p)
          deleteat!(IB, p)
          invB = inv(A[Irows,IB])
          Im = speye(length(Irows))
        else
          IB[p] = IN[q]
          deleteat!(IN, q)
          d = invB * A[Irows,IN[q]]
          Q = collect(Im)
          dp = d[p]
          d[p] = -1
          Q[:, p] = -d / dp
          invB = Q*invB
        end
        p = findfirst(IB .> n)
      end
      @assert length(Irows) > 0
      x, z, status = simplexinv(co, A[Irows,1:n], b[Irows], IB, invB)
    else
      x, z, status = simplexinv(co, A[:,1:n], b, IB, invB)
    end
  end

  return x, z, status, IB
end
