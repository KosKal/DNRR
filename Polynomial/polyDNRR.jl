function polyDNRR(x::Array{Float64}, degree::Int64, k::Int64, gibbsIter::Int64, burnIn::Int64, dpr1::Float64, dpr2::Float64, wdist::Float64, gEps1::Float64, gEps2::Float64, thLow::Float64,
  thUp::Float64, zLow::Float64, zUp::Float64, papr::Float64, pbpr::Float64, thin::Int64, samplerSeed::Int64)


  srand(samplerSeed)

  filename = "/Results"
  if k .> 0
    savelocation = string(pwd(), filename, "/seed$samplerSeed")
  else 
    savelocation = string(pwd(), filename, "/seed$samplerSeed/Reco")
  end
  mkpath(savelocation)
  # initialization
  ss = length(x)

  degree += 1 # {n-deg => (n+1)-parameters}

  nn = Int((gibbsIter - burnIn) / thin)

  theta = zeros(degree) # polynomial coeeficients vector
  sampledTheta = zeros(nn, degree) # matrix to store sample thetas
  thetaprev = copy(theta)

  N = collect(1:1:ss)
  d = collect(1:1:ss)

  noise = zeros(nn)

  clusters = zeros(nn)

  ps = zeros(nn)
  p = 1.

  x0s = zeros(nn)
  sx0 = 0.5


  if k .> 0 #introduce k-length denoised variables

    deltas = zeros(nn)
    delta = 1e-04

    ys = zeros(nn, k)
    sacost = zeros(nn,2)

    dens = zeros(k)
    for j = 1:k
        dens[j] = x[j]#rand(Uniform(zLow, zUp))
    end
   aprobs = zeros(nn)

  end

  toler = 1e-06


  display("Starting MCMC ...")

  ##########################################################################
  # MCMC
  ##########################################################################

  for iter in 1:gibbsIter

      # Compute weights up to ss
      ##########################################################################

      M = maximum(N)
      w = zeros(M)

      for i in 1:M
        w[i] = p * (1 - p)^(i - 1)
      end


      ## Sample precisions
      ##########################################################################

      Lambdas = zeros(M)

      for j in 1:M
        count = 0.0
        temp = 0.0
        if d[1] .== j
          count = count + 1.0
          temp = (x[1] - polyMap(theta, sx0)) ^ 2
        end
        for i in 2:ss
          if d[i] .== j
            count = count + 1.0
            temp = temp + (x[i] - polyMap(theta, x[i-1])) ^ 2
          end
        end
        shapel = gEps1 + 0.5 * count
        ratel = gEps2 + 0.5 * temp
        Lambdas[j] = rand(Gamma(shapel, 1. / ratel))
      end


      ## Sample the indicator variables d(i) for i=1,...,ss
      ##########################################################################

      nc1 = 0.0
      prob1 = 0.0
      for j in 1:1:N[1]
        nc1 = nc1 + Lambdas[j]^0.5 * exp(-0.5 * Lambdas[j] * (x[1] - polyMap(theta, sx0))^2)
      end
      rd = rand()
      for j in 1:1:N[1]
        prob1 = prob1 + Lambdas[j]^0.5 * exp(-0.5 * Lambdas[j] * (x[1] - polyMap(theta, sx0))^2) / nc1
        if rd .< prob1
          d[1] = j
          break
        end
      end

      for i in 2:1:ss
        nc = 0.0
        prob = 0.0
        for j in 1:1:N[i]
          nc = nc + Lambdas[j]^0.5 * exp(-0.5 * Lambdas[j] * (x[i] - polyMap(theta, x[i-1]))^2)
        end
        rd = rand()
        for j in 1:1:N[i]
          prob = prob + Lambdas[j]^0.5 * exp(-0.5 * Lambdas[j] * (x[i] -polyMap(theta, x[i-1]))^2) / nc
          if rd .< prob
            d[i] = j
            break
          end
        end
      end

      # Sample N_i auxilliary variables
      ##########################################################################

      N = etgeornd(p, d)

      ## Number of clusters
      ##########################################################################

      ucl = length(unique(d))
      # println("clusters: $ucl")
      # @printf("uCl: %s \n", ucl)

      ## Sample geometricx probability
      ##########################################################################

      p = rand(Beta(papr + 2 * ss, pbpr + sum(N) - ss))


        ## Sample the vector with the coefficients
        ##########################################################################

        for j in 1:1:length(theta)
          thetaj = copy(theta)
          thetaj[j] = 0.

          tauj = Lambdas[d[1]] * sx0 ^ (2(j-1))
          meanj = Lambdas[d[1]] * sx0 ^ (j-1) * (x[1] - polyMap(thetaj, sx0))
          for i in 2:1:ss
            tauj += Lambdas[d[i]] * x[i-1] ^ (2(j-1))
            meanj += Lambdas[d[i]] * x[i-1] ^ (j-1) * (x[i] - polyMap(thetaj, x[i-1]))
          end
          meanj /= tauj

          vj = -(2. / tauj) * log.(rand()) + (thetaprev[j] - meanj) ^ 2
          theta[j] = rand(Uniform(max(-sqrt(vj)+meanj, thLow), min(sqrt(vj)+meanj, thUp)))

        end

        thetaprev = copy(theta)


        ## Sample initial condition
        ##########################################################################

        # aux = - 2 / Lambdas[d[1]] * log.(rand()) + (x[1] - polyMap(theta, sx0)) ^ 2
        # polyRight = copy(theta)
        # polyRight[1] = polyRight[1] - (x[1] - sqrt(aux))
        # polyRight = polyRight ./ polyRight[end]
        # rootsRight = Polynomials.roots(Poly(polyRight))

        # polyLeft = copy(theta)
        # polyLeft[1] = polyLeft[1] - (x[1] + sqrt(aux))
        # polyLeft = polyLeft ./ polyLeft[end]
        # rootsLeft = Polynomials.roots(Poly(polyLeft))
        # # println("allRoots = $rootsLeft,$rootsRight")
        # allRoots = [rootsLeft rootsRight]
        # #println("allRoots = $allRoots")
        # sel = abs.(imag(allRoots)) .< toler # treat as real the roots with imagine part < ϵ
        # realRoots = sort(real(allRoots[sel]))
        # intervals = rangeIntersection(realRoots, [zLow zUp])
        # sx0 = unifmixrnd(intervals)

        # sample x0
        aux = -(2.0/Lambdas[d[1]])*log(rand()) + (x[1]-polyMap(theta,sx0))^2
        poly_right = copy(theta)
        poly_right[1] = poly_right[1] - (x[1] - sqrt(aux))
        poly_right = poly_right./poly_right[end]
        roots_right = roots(poly_right)
        poly_left = copy(theta)
        poly_left[1] = poly_left[1] - (x[1] + sqrt(aux))
        poly_left = poly_left./poly_left[end]
        roots_left = roots(poly_left)
        allroots = [roots_left roots_right]
        idx = abs.(imag(allroots)).<1e-08
        realroots = sort(real(allroots[idx]))
        intervals = rangeIntersection(realroots,[-20. 20.])
        sx0 = unifmixrnd(intervals)

        if k.> 0

          pvar = 1e-06

          mAccProbs = zeros(k) # vector of acceptance probabilities

          ## Sample xden[1]
          ##########################################################################

          yprop = rand(Normal(dens[1], sqrt(pvar)))

          accProb = min(1, exp(- 0.5 * delta * ( cfy(yprop, sx0, dens[2], theta) - cfy(dens[1], sx0, dens[2], theta) ) - 0.5 * wdist *
                                              ( (yprop - x[1])^2 - (dens[1] - x[1])^2 ) ) )

          if rand() .< accProb
            dens[1] = yprop
          end

          mAccProbs[1] = accProb

          ## Sample xden[2] - xden[ss-1]
          ##########################################################################

          for j in 2:(k - 1)

            yprop = rand(Normal(dens[j], sqrt(pvar)))

            accProb = min(1, exp(- 0.5 * delta * ( cfy(yprop, dens[j-1], dens[j+1], theta) - cfy(dens[j], dens[j-1], dens[j+1], theta) ) - 0.5 * wdist *
                                              ( (yprop - x[j])^2 - (dens[j] - x[j])^2 ) ) )

            if rand() .< accProb
              dens[j] = yprop
            end
          
            mAccProbs[j] = accProb
          
          end

          ## Sample xden[k]
          ##########################################################################

          yprop = rand(Normal(dens[k], sqrt(pvar)))

          accProb = min(1, exp(- 0.5 * delta * ( (yprop - polyMap(theta, dens[k-1]))^2 - (dens[k] - polyMap(theta, dens[k-1]))^2 ) - 0.5 * wdist *
                                              ( (yprop - x[k])^2 - (dens[k] - x[k])^2 ) ) )

          if rand() .< accProb
            dens[k] = yprop
          end

          mAccProbs[k] = accProb
          
          ## Sample denoising parameter delta
          ##########################################################################

          dtemp = (dens[1] - polyMap(theta, sx0)) ^ 2
          for i in 2:ss
            dtemp += (dens[i] - polyMap(theta, dens[i-1])) ^ 2
          end

          drate = dpr2 + 0.5 * dtemp
          delta = rand(Gamma(dpr1 + 0.5 * ss, 1./drate))


        end


        ## After Burn-In period
        ##########################################################################

      if (iter .> burnIn) & ((iter-burnIn) % thin .== 0)

        ii = Int((iter - burnIn)/thin)

        if k .> 0

          ## Calculate cost at each iteration (E_dyn)
          ##########################################################################
          
          tt = (dens[1] - polyMap(theta, sx0)) ^ 2 

          for i in 2:length(dens)
            tt +=  (dens[i] - polyMap(theta, dens[i-1])) ^ 2 
          end

          sacost[ii,1] = sqrt(1/length(dens) *tt)

          ## Calculate cost at each iteration (E0)
          ##########################################################################
          
          tt = 0.

          for i in 1:length(dens)
            tt +=  (dens[i] - x[i]) ^ 2 
          end

         sacost[ii,2] = sqrt(1/length(dens) *tt)
       end        
          ## Sample noise predictive
          ##########################################################################

          cW = cumsum(w)
          flag = rand()
          if cW[end] .< flag
            pred = rand(Normal(0, sqrt(1 / rand(Gamma(gEps1, 1 / gEps2))))) # draw from the prior
          else
            for j in 1:1:length(cW)
              if flag .< cW[j]
                 pred = rand(Normal(0, sqrt(1 / Lambdas[j])))
                 break
              end
            end
          end

          # Store values
          sampledTheta[ii, :] = thetaprev
          x0s[ii, :] = sx0
          noise[ii] = pred
          clusters[ii] = ucl
          ps[ii] = p
          if k .> 0
            aprobs[ii] = mean(mAccProbs)
            deltas[ii] = delta
            # distws[ii] = wdist
            ys[ii, :] = dens
          end

        end

        if iter % 10000 .== 0
           println("MCMC Iterations: $iter")
        end


  end

  display("MCMC finished !")

  ## Write values in .txt files - specific path
  ##########################################################################

  writedlm(string(savelocation,"/data.txt"), x, '\n')
  writedlm(string(savelocation,"/thetas.txt"), sampledTheta)
  writedlm(string(savelocation,"/x0s.txt"), x0s, '\n')
  writedlm(string(savelocation,"/noise.txt"), noise, '\n')
  writedlm(string(savelocation,"/ucl.txt"), clusters, '\n')
  writedlm(string(savelocation,"/ps.txt"), ps, '\n')

  if k.> 0
    writedlm(string(savelocation,"/probs.txt"), aprobs, '\n')
    writedlm(string(savelocation,"/ys.txt"), ys)
    writedlm(string(savelocation,"/deltas.txt"), deltas, '\n')
    writedlm(string(savelocation,"/cost.txt"), sacost)
  end

  return 0

  end
