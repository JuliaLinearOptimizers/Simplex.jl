# getλ! is equivalent to l .= c[X]
# avoids creating a vector when calling c[X]
function getλ!(λ::Vector{Float64}, c::Vector{Float64}, X::Vector{Int64})
  λ .= 0.
  for (i,j) in enumerate(X)
    @inbounds λ[i] = c[j]
  end
end

# getAcol! is equivalent to l .= A[:,col]
# avoids creating a vector when calling A[:,col]
function getAcol!(l::Vector{Float64},
                  A::SparseMatrixCSC{Float64,Int64}, col::Int64)
  l .= 0.
  for i in nzrange(A,col)
    l[A.rowval[i]] = A.nzval[i]
  end
end
function getAcol!(l::Vector{Float64}, A::SparseMatrixCSC{Float64,Int64},
                  col::Int64, Irows::Vector{Int64})
  l .= 0.; j = 0
  for i in nzrange(A,col)
    j = findfirst(Irows,A.rowval[i])
    if j != 0
      l[j] = A.nzval[i]
    end
  end
end

# same as l -= z*a where z is a number
function subdot!(l::Vector{Float64},a::SparseVector{Float64,Int64},z::Float64)
  for (i, j) in enumerate(a.nzind)
    l[j] -= z*a.nzval[i]
  end
end
function subdot!(l::Vector{Float64},a::Vector{Float64},z::Float64)
  for (i,j) in enumerate(a)
    if j != 0
      l[i] -= z*j
    end
  end
end

# returns l[perm] without allocating unnecessary memory
# returns l[invperm(perm)] if inv = true
function savepermute!(tempperm, perm, l, inv::Bool = false)
  copy!(tempperm, perm)
  !inv ? permute!!(l, tempperm) : ipermute!!(l, tempperm)
end

function shiftN!(N, q, backwards = false)
  endN = length(N)
  if backwards
    while q != 1 && N[q-1] > N[q]
      N[q], N[q-1] = N[q-1], N[q]
      q -= 1
    end
  else
    while q != endN && N[q+1] < N[q]
      N[q], N[q+1] = N[q+1], N[q]
      q += 1
    end
  end
end

# buble sort to keep N ordered
function sortN!(N::Vector{Int64}, q::Int64)
  if q == 1
    shiftN!(N,q)
  elseif q == length(N)
    shiftN!(N,q,true)
elseif N[q+1] < N[q]
    shiftN!(N,q)
  else
    shiftN!(N,q,true)
  end
end

# Bland's Rule: determine elements of r = c[N] - (l'*A[:,N])' one by one,
# q will be the first s.t. r[q] < -ϵ
function getq(c::Vector{Float64}, λ::Vector{Float64},
              A::SparseMatrixCSC{Float64,Int64}, N::Vector{Int64}, ϵ::Float64)
  for (colN,colA) in enumerate(N)
    rq = c[colA]
    for i in nzrange(A,colA)
      if λ[A.rowval[i]] != 0
        rq -= λ[A.rowval[i]]*A.nzval[i]
      end
    end
    if rq < -ϵ
      return colN
    end
  end
  return 0
end

# does A[:,p] .= a without allocating much memory
function insertAcol!(A::SparseMatrixCSC{Float64,Int64},a::Vector{Float64},p::Int64)
  for i in nzrange(A,p)
    A.nzval[i] = 0.
  end
  for i in find(a)
    A[i,p] = a[i]
  end
end
