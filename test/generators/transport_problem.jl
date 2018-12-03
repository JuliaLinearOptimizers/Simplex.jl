# Simple transport problem:
# Given m points of departure with m[i] availability
# and n points of arrival with n[i] demand, define
# the routes that minize costs of transport.

function transport_instance(m, n)
  c = rand(100.:1000, m * n)

  availArr, demandArr = rand(100.:1000, m), rand(100.:1000, n)
  diff = sum(availArr) - sum(demandArr)
  if diff > 0
    demandArr[rand(1:size(demandArr)[1])] += diff
  else
    availArr[rand(1:size(availArr)[1])] -= diff
  end

  b = [availArr ; demandArr]

  A = spzeros(m + n, m * n)
  for i = 1:m # Sum of resources leaving from each starting point
    A[i, (i - 1) * n + 1 : i * n] = 1
  end

  for i = 1:n # Sum of resources arriving to each end point
    for j = 1:m
      A[m+i, (j - 1) * n + i] = 1
    end
  end

  return A, b, c
end
