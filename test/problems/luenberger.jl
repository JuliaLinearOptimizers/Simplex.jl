function luenberger1()
  A = [2. 1 2;
       3  3 1]
  b = [4.; 3]
  c = [4.; 1; 1]
  sol = [0; 2//5; 9//5]
  return A, b, c, sol
end

function luenberger2()
  c = [4.; 1; 1]
  A = [-2. -1 -2; -3 -3 -1]
  b = [-4.; -3]
  sol = [0; 2//5; 9//5]
  return A, b, c, sol
end

function luenberger3()
  c = [-1.; 1; 0; 0; 0]
  A = [-2. 1 1 0 0; 1 -3 0 1 0; 1 1 0 0 1]
  b = [2.; 2; 4]
  sol = [7//2; 1//2; 17//2; 0; 0]
  return A, b, c, sol
end