export primalbasic!

"""
    primalbasic!(x, A, b, c)

Solves the linear program

     min  cᵀx
    s.to  Ax = b, x ≥ 0

Using the primal simplex method.
"""
function primalbasic!(
  x :: AbstractVector{T},
  A :: AbstractMatrix{T},
  b :: AbstractVector{T},
  c :: AbstractVector{T};
  IB = Int[],
  callback = (args...) -> nothing,
  max_iter = 100,
  max_time = 60.0
) where T

  ncon, nvar = size(A)

  t₀ = time()
  iter = 0

  if length(IB) == 0
    # @info("Phase I")
    xs = [x; zeros(T, ncon)]
    σ = [b[j] > 0 ? one(T) : -one(T) for j = 1:ncon]
    status, IB = primalbasic!(
      xs,
      [A  diagm(0 => σ)],
      b,
      [zeros(T, nvar); ones(T, ncon)],
      IB=collect(nvar .+ (1:ncon))
    )
    x .= xs[1:nvar]
    # @info("Phase II", IB, xs)
  else
    x[IB] .= A[:,IB] \ b
  end
  IN = setdiff(1:nvar, IB)

  Δt = time() - t₀
  solved = false
  tired = (iter > max_iter > 0 || Δt > max_time > 0)
  unbounded = false
  infeasible = any(IB .> nvar)

  while !(solved || tired || unbounded || infeasible)
    # Basic point
    B = A[:,IB]

    # Lagrange multiplier
    y = B' \ c[IB]

    # Reduced costs
    z = c[IN] - A[:,IN]' * y
    if all(z .≥ 0)
      solved = true
      continue
    end

    # Entering base
    Ientering = findall(z .< 0)
    q = rand(Ientering)
    d = B \ A[:,IN[q]]
    if all(d .≤ 0)
      unbounded = true
    end

    # Leaving
    xq⁺, p = findmin([d[j] > 0 ? x[IB[j]] / d[j] : Inf for j = 1:ncon])
    x[IB] .-= d * xq⁺
    x[IN[q]] = xq⁺

    @assert eltype(xq⁺) == T

    # Set update
    IB[p], IN[q] = IN[q], IB[p]

    iter += 1
    Δt = time() - t₀
    tired = (iter > max_iter > 0 || Δt > max_time > 0)
  end

  status = if solved
    :solved
  elseif tired
    if iter > max_iter > 0
      :max_iter
    elseif Δt > max_time > 0
      :max_time
    end
  elseif unbounded
    :unbounded
  elseif infeasible
    :infeasible
  end

  return status, IB
end