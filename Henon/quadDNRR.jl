function quadDNRR(x::Array{Float64}, k::Int64, gibbsIter::Int64, burnIn::Int64, dela::Float64, delb::Float64, wdist::Float64, gEps1::Float64, gEps2::Float64, thLow::Float64,
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

  nn = Int((gibbsIter - burnIn) / thin)

  theta = zeros(6) # polynomial coefficients vector
  sampledTheta = zeros(nn, 6) # matrix to store sample thetas

  x0s = zeros(nn, 2)
  sx01, sx02 = 0.5, 0.8

  N = collect(1:1:ss)
  d = collect(1:1:ss)
  noise = zeros(nn)
  clusters = zeros(nn)
  ps = zeros(nn)
  p = 0.5

  # Lambdas = 10000*ones(ss)

  deltas = zeros(nn)
  delta = 1e07

  dista, distb = 0., 0.

  if k .> 0 #introduce k-length denoised variables

    ys = zeros(nn, k)
    dens = zeros(k)
    for j = 1:k
        dens[j] = x[j]
    end

    aprobs = zeros(nn)
    sacost = zeros(nn,2)

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
          temp = (x[1] - henon(theta, sx01, sx02)) ^ 2
        end
        if d[2] .== j
          count = count + 1.0
          temp += (x[2] - henon(theta, x[1], sx01)) ^ 2
        end
        for i in 3:ss
          if d[i] .== j
            count = count + 1.0
            temp = temp + (x[i] - henon(theta, x[i-1], x[i-2])) ^ 2
          end
        end
        shapel = gEps1 + 0.5 * count
        ratel = gEps2 + 0.5 * temp
        # println("shape: $shapel")
        # println("rate: $ratel")
        Lambdas[j] = rand(Gamma(shapel, 1. / ratel))
      end


      ## Sample the indicator variables d(i) for i=1,...,ss
      ##########################################################################

      nc1 = 0.0
      prob1 = 0.0
      for j in 1:1:N[1]
        nc1 = nc1 + Lambdas[j]^0.5 * exp(-0.5 * Lambdas[j] * (x[1] - henon(theta, sx01, sx02))^2)
      end
      rd = rand()
      for j in 1:1:N[1]
        prob1 = prob1 + Lambdas[j]^0.5 * exp(-0.5 * Lambdas[j] * (x[1] - henon(theta, sx01, sx02))^2) / nc1
        if rd .< prob1
          d[1] = j
          break
        end
      end

      nc2 = 0.0
      prob2 = 0.0
      for j in 1:1:N[2]
        nc2 = nc2 + Lambdas[j]^0.5 * exp(-0.5 * Lambdas[j] * (x[2] - henon(theta, x[1], sx01))^2)
      end
      rd = rand()
      for j in 1:1:N[2]
        prob2 = prob2 + Lambdas[j]^0.5 * exp(-0.5 * Lambdas[j] * (x[2] - henon(theta, x[1], sx01))^2) / nc2
        if rd .< prob2
          d[2] = j
          break
        end
      end

      for i in 3:1:ss
        nc = 0.0
        prob = 0.0
        for j in 1:1:N[i]
          nc = nc + Lambdas[j]^0.5 * exp(-0.5 * Lambdas[j] * (x[i] - henon(theta, x[i-1], x[i-2]))^2)
        end
        rd = rand()
        for j in 1:1:N[i]
          prob = prob + Lambdas[j]^0.5 * exp(-0.5 * Lambdas[j] * (x[i] - henon(theta, x[i-1], x[i-2]))^2) / nc
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


      ## Sample θ₀ = theta[1]
      ##########################################################################

      muth = Lambdas[d[1]] * (x[1] - (theta[5] * sx01^2 + theta[4] * sx01 * sx02 + theta[6] * sx02^2 + theta[2] * sx01 + theta[3] * sx02)) +
             Lambdas[d[2]] * (x[2] - (theta[5] * x[1]^2 + theta[4] * x[1] * sx01 + theta[6] * sx01^2 + theta[2] * x[1] + theta[3] * sx01))

      tauth = Lambdas[d[1]] + Lambdas[d[2]]

      for j = 3:ss
        muth += Lambdas[d[j]] * (x[j] - (theta[5] * x[j-1]^2 + theta[4] * x[j-1] * x[j-2] + theta[6] * x[j-2]^2 + theta[2] * x[j-1] + theta[3] * x[j-2]))
        tauth += Lambdas[d[j]]
      end
      muth = muth / tauth

      temp = -2. / tauth * log(rand()) + (theta[1] - muth) ^ 2

      theta[1] = rand(Uniform(max(thLow, muth - temp ^ 0.5), min(thUp, muth + temp ^ 0.5)))

      #
      # ## Sample θ₁ = theta[2]
      # ##########################################################################
      #
      muth =  Lambdas[d[1]] * (x[1] * sx01 - (theta[5] * sx01^3 + theta[4] * sx01^2 * sx02 + theta[6] * sx01 * sx02^2 +
                              theta[3] * sx01 * sx02 + theta[1] * sx01)) +
              Lambdas[d[2]] * (x[2] * x[1] - (theta[5] * x[1]^3 + theta[4] * x[1]^2 * sx01 + theta[6] * x[1] * sx01^2 +
                              theta[3] * x[1] * sx01 + theta[1] * x[1]))

      tauth = Lambdas[d[1]] * sx01 ^ 2 + Lambdas[d[2]] * x[1] ^ 2

      for j = 3:ss
        muth += Lambdas[d[j]] * (x[j] * x[j-1] - (theta[5] * x[j-1]^3 + theta[4] * x[j-1]^2 * x[j-2] + theta[6] * x[j-1] * x[j-2]^2 +
                                theta[3] * x[j-1] * x[j-2] + theta[1] * x[j-1]))
        tauth += Lambdas[d[j]] * x[j - 1] ^ 2
      end
      muth = muth / tauth

      temp = -2. / tauth * log(rand()) + (theta[2] - muth) ^ 2

      theta[2] = rand(Uniform(max(thLow, muth - temp ^ 0.5), min(thUp, muth + temp ^ 0.5)))

      #
      # ## Sample θ₂ = theta[3]
      # ##########################################################################
      #
      muth =  Lambdas[d[1]] * (x[1] * sx02 - (theta[5] * sx01^2 * sx02 + theta[4] * sx01 * sx02^2 + theta[6] * sx02^3 +
                              theta[2] * sx01 * sx02 + theta[1] * sx02)) +
              Lambdas[d[2]] * (x[2] * sx01 - (theta[5] * x[1]^2 * sx01 + theta[4] * x[1] * sx01^2 + theta[6] * sx01^3 +
                              theta[2] * x[1] * sx01 + theta[1] * sx01))

      tauth = Lambdas[d[1]] * sx02 ^ 2 + Lambdas[d[2]] * sx01 ^ 2

      for j = 3:ss
        muth += Lambdas[d[j]] * (x[j] * x[j-2] - (theta[5] * x[j-1]^2 * x[j-2] + theta[4] * x[j-1] * x[j-2]^2 + theta[6] * x[j-2]^3 +
                                theta[2] * x[j-1] * x[j-2] + theta[1] * x[j-2]))
        tauth += Lambdas[d[j]] * x[j - 2] ^ 2
      end

      muth = muth / tauth

      temp = -2. / tauth * log(rand()) + (theta[3] - muth) ^ 2

      theta[3] = rand(Uniform(max(thLow, muth - temp ^ 0.5), min(thUp, muth + temp ^ 0.5)))

      #
      # ## Sample θ₃ = theta[4]
      # ##########################################################################
      #
      muth = Lambdas[d[1]] * (x[1] * sx01 * sx02 - (theta[5] * sx01^3 * sx02 + theta[6] * sx01 * sx02^3 + theta[2] * sx01^2 * sx02 +
                              theta[3] * sx01 * sx02^2 + theta[1] * sx01 * sx02)) +
             Lambdas[d[2]] * (x[2] * x[1] * sx01 - (theta[5] * x[1]^3 * sx01 + theta[6] * x[1] * sx01^3 + theta[2] * x[1]^2 * sx01 +
                              theta[3] * x[1] * sx01^2 + theta[1] * x[1] * sx01))

      tauth = Lambdas[d[1]] * sx01 ^ 2 * sx02 ^ 2 + Lambdas[d[2]] * x[1] ^ 2 *sx01 ^ 2

      for j = 3:ss
        muth += Lambdas[d[j]] * (x[j] * x[j-1] * x[j-2] - (theta[5] * x[j-1]^3 * x[j-2] + theta[6] * x[j-1] * x[j-2]^3 + theta[2] * x[j-1]^2 * x[j-2] +
                                theta[3] * x[j-1] * x[j-2]^2 + theta[1] * x[j-1] * x[j-2]))
        tauth += Lambdas[d[j]] * x[j - 1] ^ 2 * x[j - 2] ^ 2
      end

      muth = muth / tauth

      temp = -2. / tauth * log(rand()) + (theta[4] - muth) ^ 2

      theta[4] = rand(Uniform(max(thLow, muth - temp ^ 0.5), min(thUp, muth + temp ^ 0.5)))

      #
      # ## Sample θ₄ = theta[5]
      # ##########################################################################
      #
      muth = Lambdas[d[1]] * (x[1] * sx01^2 - (theta[4] * sx01^3 * sx02 + theta[6] * sx01^2 * sx02^2 + theta[2] * sx01^3 +
                              theta[3] * sx01^2 * sx02 + theta[1] * sx01^2)) +
             Lambdas[d[2]] * (x[2] * x[1]^2 - (theta[4] * x[1]^3 * sx01 + theta[6] * x[1]^2 * sx01^2 + theta[2] * x[1]^3 +
                              theta[3] * x[1]^2 * sx01 + theta[1] * x[1]^2))

      tauth = Lambdas[d[1]] * sx01 ^ 4 + Lambdas[d[2]] * x[1] ^ 4

      for j = 3:ss
        muth += Lambdas[d[j]] * (x[j] * x[j-1]^2 - (theta[4] * x[j-1]^3 * x[j-2] + theta[6] * x[j-1]^2 * x[j-2]^2 + theta[2] * x[j-1]^3 +
                                theta[3] * x[j-1]^2 * x[j-2] + theta[1] * x[j-1]^2))
        tauth += Lambdas[d[j]] * x[j - 1] ^ 4
      end

      muth = muth / tauth

      temp = -2. / tauth * log(rand()) + (theta[5] - muth) ^ 2

      theta[5] = rand(Uniform(max(thLow, muth - temp ^ 0.5), min(thUp, muth + temp ^ 0.5)))

      #
      # ## Sample θ₅ = theta[6]
      # ##########################################################################
      #
      muth = Lambdas[d[1]] * (x[1] * sx02^2 - (theta[5] * sx01^2 * sx02^2 + theta[4] * sx01 * sx02^3 + theta[2] * sx01 * sx02^2 +
                              theta[3] * sx02^3 + theta[1] * sx02^2)) +
             Lambdas[d[2]] * (x[2] * sx01^2 - (theta[5] * x[1]^2 * sx01^2 + theta[4] * x[1] * sx01^3 + theta[2] * x[1] * sx01^2 +
                              theta[3] * sx01^3 + theta[1] * sx01^2))

      tauth = Lambdas[d[1]] * sx02 ^ 4 + Lambdas[d[2]] * sx01 ^ 4

      for j = 3:ss
        muth += Lambdas[d[j]] * (x[j] * x[j-2]^2 - (theta[5] * x[j-1]^2 * x[j-2]^2 + theta[4] * x[j-1] * x[j-2]^3 + theta[2] * x[j-1] * x[j-2]^2 +
                                theta[3] * x[j-2]^3 + theta[1] * x[j-2]^2))
        tauth += Lambdas[d[j]] * x[j - 2] ^ 4
      end

      muth = muth / tauth

      temp = -2. / tauth * log(rand()) + (theta[6] - muth) ^ 2

      theta[6] = rand(Uniform(max(thLow, muth - temp ^ 0.5), min(thUp, muth + temp ^ 0.5)))

      ## Sample initial conditions x₀,y₀ = x0,x00
      ##########################################################################

      # x₀

      a4 = Lambdas[d[1]] * theta[5] ^ 2 + Lambdas[d[2]] * theta[6] ^ 2

      a3 = 2. * Lambdas[d[1]] * (theta[4] * theta[5] * sx02 + theta[2] * theta[5]) + 2. * Lambdas[d[2]] * (theta[4] * theta[6] * x[1] + theta[3] * theta[6])

      a2 = Lambdas[d[1]] * (theta[4]^2 * sx02^2 + 2. * theta[5] * theta[6] * sx02^2 + 2. * theta[2] * theta[4] * sx02 + 2. * theta[3] * theta[5] * sx02 -
                            2. * theta[5] * x[1] + 2. * theta[1] * theta[5] + theta[2]^2) +
           Lambdas[d[2]] * (theta[4]^2 * x[1]^2 + 2. * theta[5] * theta[6] * x[1]^2 + 2 * theta[2] * theta[6] * x[1] + 2 * theta[3] * theta[4] * x[1] -
                            2. * theta[6] * x[2] + 2. * theta[1] * theta[6] + theta[3]^2)

      a1 = 2. * (Lambdas[d[1]] * (theta[4] * theta[6] * sx02^3 + (theta[2] * theta[6] + theta[3] * theta[4]) * sx02^2 + (-theta[4] * x[1] + theta[1] * theta[4] + theta[2] * theta[3]) * sx02 -
                                  theta[2] * x[1] + theta[1] * theta[2]) +
                Lambdas[d[2]] * (theta[4] * theta[5] * x[1]^3 + (theta[2] * theta[4] + theta[3] * theta[5]) * x[1]^2 -
                                theta[4] * x[1] * x[2] + (theta[1] * theta[4] + theta[2] * theta[3]) * x[1] - theta[3] * x[2] + theta[1] * theta[3]))

      aux = -2. * log.(rand()) + a1 * sx01 + a2 * sx01 ^ 2 + a3 * sx01 ^ 3 + a4 * sx01 ^ 4

      poly = [-aux / a4; a1 / a4; a2 / a4; a3 / a4; 1.0]
      #allRoots = roots(poly, polish=true, epsilon=1e-20)
      allRoots = roots(poly)

      #println("allRoots = $allRoots")
      # sel = imag(allRoots) .== 0.0
      # realRoots = sort(real(allRoots[sel]))

      sel = abs.(imag(allRoots)) .< toler # treat as real the roots with imagine part < ϵ
      realRoots = sort(real(allRoots[sel]))
      intervals1 = rangeIntersection(realRoots, [zLow; zUp])
      sx01 = unifmixrnd(intervals1)
      #
      # x00 = y₀

      a4  = Lambdas[d[1]] * theta[6] ^ 2

      a3 = 2. * theta[6] * Lambdas[d[1]] * (theta[4] * sx01 + theta[3])

      a2 = Lambdas[d[1]] * ( -2. * theta[6] * (x[1] - theta[1] - theta[2] * sx01 - theta[5] * sx01^2) + (theta[4] * sx01 + theta[3])^2 )

      a1 = -2. * Lambdas[d[1]] * (x[1] - theta[1] - theta[5] * sx01^2 - theta[2] * sx01) * (theta[4] * sx01 + theta[3])

      aux = -2. * log.(rand()) + a1 * sx02 + a2 * sx02 ^ 2 + a3 * sx02 ^ 3 + a4 * sx02 ^ 4

      poly = [-aux / a4; a1 / a4; a2 / a4; a3 / a4; 1.0]
      allRoots = roots(poly)
      sel = abs.(imag(allRoots)) .< toler # treat as real the roots with imagine part < ϵ
      # sel = imag(allRoots) .== 0.0 # treat as real the roots with imagine part < ϵ
      realRoots = sort(real(allRoots[sel]))
      intervals2 = rangeIntersection(realRoots, [zLow; zUp])
      sx02 = unifmixrnd(intervals2)


      if k.> 0

        pvar = 1e-06#1e-06

        ## Sample xden[1]
        ##########################################################################

        mAccProbs = zeros(k)

        yprop = rand(Normal(dens[1], sqrt(pvar)))

        accProb = min(1, exp(- 0.5 * delta * ( cfy(yprop, sx02, sx01, dens[2], dens[3], theta) - cfy(dens[1], sx02, sx01, dens[2], dens[3], theta) ) - 0.5 * wdist *
                                              ( (yprop - x[1])^2 - (dens[1] - x[1])^2 ) ) )

        if rand() .< accProb
          dens[1] = yprop
        end

        mAccProbs[1] = accProb

        ## Sample xden[2]
        ##########################################################################

        yprop = rand(Normal(dens[2], sqrt(pvar)))

        accProb = min(1, exp(- 0.5 * delta * ( cfy(yprop, sx01, dens[1], dens[3], dens[4], theta) - cfy(dens[2], sx01, dens[1], dens[3], dens[4], theta) ) - 0.5 * wdist *
                                              ( (yprop - x[2])^2 - (dens[2] - x[2])^2 ) ) )

        if rand() .< accProb
          dens[2] = yprop
        end

        mAccProbs[2] = accProb

        ## Sample xden[3] - xden[ss-2]
        ##########################################################################

        for j in 3:(k - 2)

          yprop = rand(Normal(dens[j], sqrt(pvar)))

          accProb = min(1, exp(- 0.5 * delta * ( cfy(yprop, dens[j-2], dens[j-1], dens[j+1], dens[j+2], theta) - cfy(dens[j], dens[j-2], dens[j-1], dens[j+1], dens[j+2], theta) ) - 0.5 * wdist *
                                                ( (yprop - x[j])^2 - (dens[j] - x[j])^2 ) ) )

          if rand() .< accProb
            dens[j] = yprop
          end

          mAccProbs[j] = accProb

        end

        
        ## Sample xden[k-1]
        ##########################################################################

        yprop = rand(Normal(dens[k-1], sqrt(pvar)))

        accProb = min(1, exp(- 0.5 * delta * ( (yprop - henon(theta, dens[k-2], dens[k-3])) ^ 2 + (dens[k] - henon(theta, yprop, dens[k-2])) ^ 2
                              - (dens[k-1] - henon(theta, dens[k-2], dens[k-3])) ^ 2 - (dens[k] - henon(theta, dens[k-1], dens[k-2])) ^ 2  ) - 0.5 * wdist *
                               ( (yprop - x[k-1])^2 - (dens[k-1] - x[k-1])^2 ) ) )
       
        if rand() .< accProb
          dens[k-1] = yprop
        end

        mAccProbs[k-1] = accProb

        ## Sample xden[k]
        ##########################################################################

        yprop = rand(Normal(dens[k], sqrt(pvar)))

        accProb = min(1, exp(- 0.5 * delta * ( (yprop - henon(theta, dens[k-1], dens[k-2])) ^ 2 - (dens[k] - henon(theta, dens[k-1], dens[k-2])) ^ 2 ) - 0.5 * wdist *
                                                ( (yprop - x[k])^2 - (dens[k] - x[k])^2 ) ) )

        if rand() .< accProb
          dens[k] = yprop
        end

        mAccProbs[k] = accProb

        ## Sample denoising parameter delta
        ##########################################################################

        dtemp = (dens[1] - henon(theta, sx01, sx02)) ^ 2 + (dens[2] - henon(theta, dens[1], sx01)) ^ 2
        for i in 3:ss
          dtemp += (dens[i] - henon(theta, dens[i-1], dens[i-2])) ^ 2
        end

        drate = delb + 0.5 * dtemp
        delta = rand(Gamma(dela + 0.5 * ss, 1./drate))

      end


      ## After Burn-In period
      ###############################wwwwwwwww###########################################

      if (iter .> burnIn) & ((iter-burnIn) % thin .== 0)

        ii = Int((iter - burnIn)/thin)

        if k .> 0

          ## Calculate cost at each iteration (E_dyn)
          ##########################################################################
          
          tt = (dens[1] - henon(theta, sx01, sx02)) ^ 2 + (dens[2] - henon(theta, dens[1], sx01)) ^ 2

          for i in 3:length(dens)
            tt +=  (dens[i] - henon(theta, dens[i-1], dens[i-2])) ^ 2 
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
        sampledTheta[ii, :] = theta
        x0s[ii, :] = [sx01, sx02]
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

	display("... MCMC finished !")

	## Write values in .txt files - specific path
	##########################################################################

	writedlm(string(savelocation,"/data.txt"), x, '\n')
	writedlm(string(savelocation, "/thetas.txt"), sampledTheta)
	writedlm(string(savelocation,"/x0s.txt"), x0s)
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
