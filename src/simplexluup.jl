export simplexluup

function simplexluup(c::Vector{Float64}, A::SparseMatrixCSC{Float64,Int64}, b::Vector{Float64},
                     𝔹=0, L=0, U = 0,
                     prow=Vector{Int64}(A.m), Rs=Vector{Float64}(A.m),
                     xB::Vector{Float64}=Float64[]; max_iter::Int64 = 20000, maxups::Int64 = 10)

  ϵ = norm(c)*1e-13
  m, n = A.m, A.n # preparations
  iter = 0; ups = 0; maxed = false
  P, MP = Vector{Vector{Int64}}(maxups), Vector{SparseVector{Float64,Int64}}(maxups)
  lnz = Ref{Int64}(); unz = Ref{Int64}(); nz_diag = Ref{Int64}()
  n_row = Ref{Int64}(); n_col = Ref{Int64}()
  Lp = Vector{Int64}(m + 1); Up = Vector{Int64}(m + 1)
  pcol = Vector{Int64}(m); tempperm = Vector{Int64}(m)
  signb = sign.(b)
  for (i,j) in enumerate(signb)
    (j >= 0) ? signb[i] = 1 : signb[i] = -1
  end
  if (U == 0) U = spdiagm(signb) end

  w = Array{Float64,1}(m); d = Vector{Float64}(m)
  λ = Array{Float64,1}(m); Ucolp = Array{Float64,1}(m)

  if 𝔹 == 0 # construct artificial problem
    artificial = true
    Ao = A
    A = [A U]
    Utri = UpperTriangular(U)
    Ut = transpose(U); Uttri = LowerTriangular(Ut)
    𝔹 = collect(n+1:n+m); ℕ = collect(1:n) # artificial indexes
    ca = [zeros(n); ones(m)];
    λ .= signb
    xB = abs.(b) # solution in current basis
  else
    artificial = false
    ℕ = setdiff(1:n, 𝔹)
    getλ!(λ,c,𝔹)
    if L == 0
      F = lufact(A[:,𝔹]) # (Rs.*A)[prow,pcol] * x[pcol] = b[prow]
      xB = F\b
      L, U, prow, pcol, Rs = F[:(:)]
      savepermute!(tempperm, pcol, 𝔹)
      permute!!(xB, pcol) # pcol is lost
    end
    Utri = UpperTriangular(U); Ltri = LowerTriangular(L)
    Ut = transpose(U); Lt = transpose(L);
    Uttri = LowerTriangular(Ut); Lttri = UpperTriangular(Lt)
    Lt = transpose(L); Lttri = UpperTriangular(Lt)
    A_ldiv_B!(Uttri,λ); A_ldiv_B!(Lttri,λ)
    savepermute!(tempperm, prow, λ, true)
    λ .= λ.*Rs
  end
  artificial? q = getq(ca,λ,A,ℕ,ϵ) : q = getq(c,λ,A,ℕ,ϵ)
  status = :Optimal

  # simplex search
  while !(q == 0 || iter > max_iter)

      iter += 1
      getAcol!(w,A,ℕ[q])
      if L != 0
        w .= w.*Rs
        savepermute!(tempperm, prow, w)
        A_ldiv_B!(Ltri,w)
      end
      for k in 1:ups
        savepermute!(tempperm, P[k], w)
        w[end] -= dot(MP[k], w)
      end
      d .= w
      A_ldiv_B!(Utri,d)
      xq = Inf
      for k in 1:m # find min xB/d s.t. d .> 0
        if d[k] >= ϵ
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
    𝔹[p], ℕ[q] = ℕ[q], 𝔹[p] # column change: update indexes
    sortN!(ℕ,q)

    if ups >= maxups # reset LU
      maxed = true
      F = lufact(A[:,𝔹])
      ccall(("umfpack_dl_get_lunz",:libumfpack), Int64,(Ptr{Int64},Ptr{Int64},
            Ptr{Int64},Ptr{Int64},Ptr{Int64},Ptr{Void}),
            lnz,unz,n_row,n_col,nz_diag,F.numeric)
      Lj = Vector{Int64}(lnz[]); Lx = Vector{Float64}(lnz[])
      Ui = Vector{Int64}(unz[]); Ux = Vector{Float64}(unz[])
      ccall(("umfpack_dl_get_numeric",:libumfpack),Int64, (Ptr{Int64},
             Ptr{Int64},Ptr{Float64},Ptr{Int64},Ptr{Int64},Ptr{Float64},
             Ptr{Int64},Ptr{Int64},Ptr{Void},Ref{Int64},Ptr{Float64},
             Ptr{Void}),Lp,Lj,Lx,Up,Ui,Ux,prow,pcol,C_NULL,0, Rs, F.numeric)
      if L == 0
        Lt = SparseMatrixCSC(m, m, increment!(Lp), increment!(Lj), Lx)
        L = transpose(Lt)
        Ltri = LowerTriangular(L); Lttri = UpperTriangular(Lt)
      else
        copy!(Lt, SparseMatrixCSC(m, m, increment!(Lp), increment!(Lj), Lx))
        copy!(L,Lt); transpose!(L,Lt)
      end
      copy!(U, SparseMatrixCSC(m, m, increment!(Up), increment!(Ui), Ux))
      increment!(prow); increment!(pcol)
      savepermute!(tempperm, pcol, 𝔹)
      permute!!(xB, pcol) # pcol is lost
      ups = 0
    else # update LU
      insertAcol!(U,w,p)
      if findlast(w) > p
        ups += 1
        if maxed
          P[ups] .= 1:m; reverse!(reverse!(P[ups],p,m),p,m-1)
        else
          P[ups] = 1:m; reverse!(reverse!(P[ups],p,m),p,m-1)
        end
        maxed? MP[ups] .= spzeros(m) : MP[ups] = spzeros(m)
        if length(Ut.rowval) < nnz(U)
          copy!(Ut, U)
        end
        halfperm!(Ut, U, P[ups])
        getAcol!(Ucolp,Ut,p)
        for i in p:m-1
          if Ucolp[i] != 0.
            (MP[ups])[i] = Ucolp[i]/Ut[i,i+1]
            for j in nzrange(Ut,i+1)
              Ucolp[Ut.rowval[j]] -= (MP[ups])[i]*Ut.nzval[j]
            end
          end
        end
        Ut[end,p] = Ucolp[end]
        halfperm!(U, Ut, P[ups])
        savepermute!(tempperm, P[ups], 𝔹)
        savepermute!(tempperm, P[ups], xB)
      end
    end

    # check optimality and choose variable to leave basis if necessary
    artificial? getλ!(λ,ca,𝔹) : getλ!(λ,c,𝔹)
    Uend = length(U.rowval)
    if length(Ut.rowval) < nnz(U)
      copy!(Ut,U)
    end
    transpose!(Ut,U)
    A_ldiv_B!(Uttri,λ)
    for j in 1:ups
      subdot!(λ,MP[ups-j+1],λ[end])
      savepermute!(tempperm, P[ups-j+1], λ, true)
    end
    if L != 0
      A_ldiv_B!(Lttri,λ)
      savepermute!(tempperm, prow, λ, true)
      λ .= λ.*Rs
    end
    artificial? q = getq(ca,λ,A,ℕ,ϵ) : q = getq(c,λ,A,ℕ,ϵ)
  end

  if iter >= max_iter
    status = :UserLimit
  end

  x = zeros(n) # finalization
  if !artificial
    x[𝔹] = xB
    z = dot(c, x)
  else
    Irows = collect(1:m)
    ℕ = setdiff(ℕ,n+1:n+m)
    artificial? getλ!(λ,ca,𝔹) : getλ!(λ,c,𝔹)
    if dot(xB, λ)/norm(xB) > 1e-6
      status = (iter >= max_iter) ? :UserLimit : :Infeasible
      I = find(𝔹 .<= n - m)
      x[𝔹[I]] = xB[I]
      z = dot(c, x)
    elseif maximum(𝔹) > n # check for artificial variables in basis
      # remove artificial variables from basis
      p, pind = findmax(𝔹)
      Ap = Array{Float64, 1}(m)
      while p > n
        q = 1
        getAcol!(Ap,A,p,Irows)
        if L != 0 savepermute!(tempperm, prow, Ap) end
        for j in 1:ups
          savepermute!(tempperm, P[j], Ap)
        end
        PivotAp = findfirst(Ap)
        while q <= length(ℕ) # searching for columns to substitute artificials ℕ basis
          getAcol!(d,A,ℕ[q],Irows)
          if L != 0
            d .= d.*Rs
            savepermute!(tempperm, prow, d)
            A_ldiv_B!(Ltri,d)
          end
          for j in 1:ups
            savepermute!(tempperm, P[j], d)
            d[end] -= dot(MP[j], d)
          end
          A_ldiv_B!(Utri,d)
          (abs(d[PivotAp]) >= ϵ) ? break : q += 1
        end
        if q > length(ℕ)
          deleteat!(Irows, findfirst(A[Irows,p]))
          deleteat!(𝔹, pind); deleteat!(xB, pind)
          deleteat!(Ap, pind); deleteat!(d, pind)
        else
          𝔹[pind] = ℕ[q]
          deleteat!(ℕ, q)
        end
        F = lufact(A[Irows,𝔹])
        L, U, prow, pcol, Rs = F[:(:)]
        Utri = UpperTriangular(U); Ltri = LowerTriangular(L)
        ups = 0
        savepermute!(tempperm, pcol, 𝔹)
        permute!!(xB, pcol) # pcol is lost
        p, pind = findmax(𝔹)
      end
      return simplexluup(c, A[Irows,1:n], b[Irows], 𝔹, L, U, prow, Rs, xB)
    else
      return simplexluup(c, Ao, b, 𝔹, L, U, prow, Rs, xB)
    end
  end
  return x, z, status
end
