function henon(θ,x,y)
  return(θ[1] + θ[2] * x + θ[3] * y + θ[4] * x * y + θ[5] * x ^ 2 + θ[6] * y ^ 2)
end


function noiseMix3(p1::Float64, lam1::Float64, lam2::Float64)
  if rand() < p1
    z = rand(Normal(0, sqrt(1 / lam1)))
  else
    z = rand(Normal(0, sqrt(1 / lam2)))
  end
  return z
end

function genData(npar, n, Θ, x0, x00, p1, lam1, lam2, seed)
  srand(seed)
  x = zeros(n)

  if npar .== 0
    x[1] = henon(Θ,x0,x00) + rand(Normal(0, sqrt(1 / lam1)))
    x[2] = henon(Θ,x[1],x0) + rand(Normal(0, sqrt(1 / lam1)))
    for i in 3:1:n
      x[i] = henon(Θ,x[i-1],x[i-2]) + rand(Normal(0, sqrt(1 / lam1)))
    end
  elseif npar .== 1
    x[1] = henon(Θ,x0,x00) + noiseMix1()
    x[2] = henon(Θ,x[1],x0) + noiseMix1()
    for i in 3:1:n
      x[i] = henon(Θ,x[i-1],x[i-2]) + noiseMix1()
    end
  elseif npar .== 2
    x[1] = henon(Θ,x0,x00) + noiseMix2(4, 1e-03)
    x[2] = henon(Θ,x[1],x0) + noiseMix2(4, 1e-03)
    for i in 3:1:n
      x[i] = henon(Θ,x[i-1],x[i-2]) + noiseMix2(4, 1e-03)
    end
  else
    x[1] = henon(Θ,x0,x00) + noiseMix3(p1, lam1, lam2)
    x[2] = henon(Θ,x[1],x0) + noiseMix3(p1, lam1, lam2)
    for i in 3:1:n
      x[i] = henon(Θ,x[i-1],x[i-2]) + noiseMix3(p1, lam1, lam2)
    end
  end

  return x
end

function unifmixrnd(x)
  local U::Float64
  n = length(x)
  nc = 0
  for i in 1:2:n
    nc += x[i + 1] - x[i]
  end

  U = 0
  u = rand()
  prob = 0
  for i in 1:2:n
    prob += (x[i + 1] - x[i]) / nc
    if u .< prob
      U = x[i] + (x[i + 1] - x[i]) * rand()
      break
    end
  end
  return(U)
end

function rangeIntersection(A::Array{Float64}, B::Array{Float64})
#=
Purpose: Range/interval intersection

 A and B two ranges of closed intervals written
 as vectors [lowerbound1 upperbound1 lowerbound2 upperbound2]
 or as matrix [lowerbound1, lowerbound2, lowerboundn;
               upperbound1, upperbound2, upperboundn]
 A and B have to be sorted in ascending order

 out is the mathematical intersection A n B

 EXAMPLE USAGE:
   >> out=rangeIntersection([1 3 5 9],[2 9])
   	out =  [2 3 5 9]
   >> out=rangeIntersection([40 44 55 58], [42 49 50 52])
   	out =  [42 44]
=#

# Allocate, as we don't know yet the size, we assume the largest case
  out1 = zeros(length(B)+(length(A)-2))
  k = 1

  while isempty(A) .== 0 && isempty(B) .== 0
  # make sure that first is ahead second
    if A[1] .> B[1]
      temp = copy(B)
      B = copy(A)
      A = copy(temp)
    end

    if A[2] .< B[1]
      A = copy(A[3:end])
      continue
    elseif A[2] .== B[1]
      out1[k] = B[1]
      out1[k + 1] = B[1]
      k = k + 2

      A = A[3:end]
      continue
    else
      if A[2] .== B[2]
        out1[k] = B[1]
        out1[k+1] = B[2]
        k = k + 2

        A = A[3:end]
        B = B[3:end]

      elseif A[2] .< B[2]
        out1[k] = B[1]
        out1[k+1] = A[2]
        k = k + 2

        A = A[3:end]
      else
        out1[k] = B[1]
        out1[k+1] = B[2]
        k = k + 2

        B = B[3:end]
      end
     end
    end

  # Remove the tails
  out = copy(out1[1:k-1])
  return(out)
end

function etgeornd(p::Float64, k::Vector{Int64})
local z::Vector{Int64}
  z = floor.(log.(rand(length(k))) ./ log.(1 - p)) + k
  return(z)
end


function cfy(yi, ypp, yp, yn, ynn, theta)
  return ((yi - henon(theta, yp, ypp)) ^ 2 + (yn - henon(theta, yi, yp)) ^ 2 + (ynn - henon(theta, yn, yi)) ^ 2)
end

