# TASK 1, Basil

function solve_basil(
    N::Int,       
    X::Float64,   
    c::Float64,  
    q::Float64,   
    pmin::Float64,
    pmax::Float64,
    step::Float64 = 0.1  
)

    pvals = collect(pmin:step:pmax)  
    M = length(pvals)                
    prob_price = 1.0 / M            

    v   = zeros(Float64, N+1)             
    vT  = zeros(Float64, N+1)             
    vA  = zeros(Float64, N+1)            
    sigma_approach = zeros(Int, N+1)      
    sigma_buy      = zeros(Int, N+1, M)   
    
    vT[N+1] = -c * N   
    
    vA[N+1] = -c + q * prob_price * sum( 
                  max( X - p - c*N, -c*N ) for p in pvals 
               ) + (1-q)*(-c*N)
    v[N+1] = max(vT[N+1], vA[N+1])
    sigma_approach[N+1] = vA[N+1] > vT[N+1] ? 1 : 0

    
    if sigma_approach[N+1] == 1
        for (pidx, p) in enumerate(pvals)
            sigma_buy[N+1,pidx] = (X - p - c*N >= -c*N) ? 1 : 0
        end
    end

    # Backward induction for n = N-1, N-2, ..., 0
    for n in (N-1):-1:0
        vT[n+1] = -c * n

       
        expected_if_orchid = prob_price * sum(
            max( X - p - c*(n+1), v[n+2] )  for p in pvals
        )
        vA[n+1] = -c + q * expected_if_orchid + (1-q)*v[n+2]

        v[n+1] = max(vT[n+1], vA[n+1])

        sigma_approach[n+1] = (vA[n+1] > vT[n+1]) ? 1 : 0

        if sigma_approach[n+1] == 1
            for (pidx, p) in enumerate(pvals)
                sigma_buy[n+1,pidx] = ( X - p - c*(n+1) >= v[n+2] ) ? 1 : 0
            end
        end
    end

    return v, sigma_approach, sigma_buy, pvals
end



function main()
    N    = 50      
    X    = 50.0    
    c    = 0.5     
    q    = 0.15   
    pmin = 10.0
    pmax = 100.0

    v, sigma_approach, sigma_buy, pvals = solve_basil(N, X, c, q, pmin, pmax)

    

    
    prob_state = zeros(Float64, N+1)
    prob_state[1] = 1.0  

    buy_prob = 0.0
    price_times_prob = 0.0  

    for n in 0:N
        idx = n+1
        if prob_state[idx] == 0.0
            continue
        end

        if sigma_approach[idx] == 1
            

            frac_see_orchid = prob_state[idx] * q
            frac_no_orchid  = prob_state[idx] * (1.0 - q)

            for (pidx,p) in enumerate(pvals)
                if sigma_buy[idx,pidx] == 1
                    buy_here = frac_see_orchid * (1.0/length(pvals))
                    buy_prob += buy_here
                    price_times_prob += p * buy_here
                else
                    
                    frac_reject = frac_see_orchid * (1.0/length(pvals))
                    if n < N
                        prob_state[idx+1] += frac_reject
                    end
                end
            end

            if n < N
                prob_state[idx+1] += frac_no_orchid
            end

        else
            
        end
    end

    println("Probability Basil buys the orchid = ", buy_prob)

    if buy_prob > 0
        println("Expected price paid (given a purchase) = ", price_times_prob / buy_prob)
    else
        println("No purchase occurs under this policy (buy_prob=0).")
    end

    

    prob_state .= 0.0
    prob_state[1] = 1.0
    expected_n_approached = 0.0

    for n in 0:N
        idx = n+1
        if prob_state[idx] == 0.0
            continue
        end
        if sigma_approach[idx] == 1
            expected_n_approached += prob_state[idx]  
            frac_see_orchid = prob_state[idx] * q
            frac_no_orchid  = prob_state[idx] * (1.0 - q)
            
            buy_fraction = frac_see_orchid * sum(sigma_buy[idx,:]) / length(pvals)
            not_buy_fraction = frac_see_orchid - buy_fraction
            if n < N
                prob_state[idx+1] += frac_no_orchid + not_buy_fraction
            end
        end
    end

    println("Expected number of vendors approached = ", expected_n_approached)

    
    return
end

main()


#Task 2 :(
    using Plots

    beta = 0.95
    c = 2.0
    wgrid = 5.0:1.0:20.0
    π = fill(1.0/length(wgrid), length(wgrid))
    pvals = 0:0.01:1.0
    wstar_vals = similar(pvals)
    q_vals = similar(pvals)
    dur_vals = similar(pvals)
    
    function solve_model(p, wgrid, π, c, β)
        barVU = c/(1-β)
        diff = 1.0
        tol = 1e-8
        while diff>tol
            VE = (wgrid .+ β*p*barVU)./(1.0 - β*(1.0-p))
            VU = max.(VE, barVU)
            new_barVU = c + β*dot(VU, π)
            diff = abs(new_barVU - barVU)
            barVU = new_barVU
        end
        VE = (wgrid .+ β*p*barVU)./(1.0 - β*(1.0-p))
        wstar_index = findfirst(i->VE[i]>=barVU, eachindex(wgrid))
        if wstar_index == nothing
            wstar = maximum(wgrid)
        else
            wstar = wgrid[wstar_index]
        end
        q = sum(π[i] for i in eachindex(wgrid) if wgrid[i]>=wstar)
        return wstar, q, 1/q
    end
    
    for (i,pp) in enumerate(pvals)
        wstar_vals[i], q_vals[i], dur_vals[i] = solve_model(pp, wgrid, π, c, beta)
    end
    
    plot(pvals, wstar_vals, xlabel="p", ylabel="w*", legend=false)
    plot(pvals, q_vals, xlabel="p", ylabel="q", legend=false)
    plot(pvals, dur_vals, xlabel="p", ylabel="Expected Unemployment Duration", legend=false)


# TASK 3, NGM

using Plots

function f1_kstar(β,α,δ) #function f had to be renamed f1 for all 4 tasks to run back to back
    return ((α*β)/(1-β*(1-δ)))^(1/(1-α))
end

function f1_f(k,α)
    return k^α
end

function f1_simulate(β,α,δ,γ,k₀,T)
    K=zeros(T+1)
    C=zeros(T+1)
    I=zeros(T+1)
    Y=zeros(T+1)
    K[1]=k₀
    for t in 1:T
        Y[t]=f1_f(K[t],α)
        if γ==1
            C[t]=Y[t]+(1-δ)*K[t] - (β*α*Y[t])
        else
            C[t]=Y[t]+(1-δ)*K[t] - (β*α*Y[t])^(1/γ)
        end
        I[t]=Y[t]+(1-δ)*K[t]-C[t]
        K[t+1]=Y[t]+(1-δ)*K[t]-C[t]
    end
    Y[T+1]=f1_f(K[T+1],α)
    if γ==1
        C[T+1]=Y[T+1]+(1-δ)*K[T+1] - (β*α*Y[T+1])
    else
        C[T+1]=Y[T+1]+(1-δ)*K[T+1] - (β*α*Y[T+1])^(1/γ)
    end
    I[T+1]=Y[T+1]+(1-δ)*K[T+1]-C[T+1]
    return K,C,I,Y
end

function f1_convergence_table(β,α,δ,γvals)
    kst=f1_kstar(β,α,δ)
    tbl=[]
    for γ in γvals
        K,C,I,Y=f1_simulate(β,α,δ,γ,0.5*kst,200)
        m=0
        for t in 1:length(K)
            if kst-K[t]<0.5*(kst-0.5*kst)
                m=t-1
                break
            end
        end
        push!(tbl,(γ,m))
    end
    return tbl
end

function f1_plot_panels(β,α,δ,γvals)
    p=plot(layout=(2,2))
    for γ in γvals
        K,C,I,Y=f1_simulate(β,α,δ,γ,0.5*f1_kstar(β,α,δ),50)
        plot!(p[1],K,label="γ=$γ")
        plot!(p[2],Y,label="γ=$γ")
        plot!(p[3],I./Y,label="γ=$γ")
        plot!(p[4],C./Y,label="γ=$γ")
    end
    display(p)
end

β=0.95
α=0.3
δ=0.05
γvals=[0.5,1.0,2.0]
display(f1_convergence_table(β,α,δ,γvals))
f1_plot_panels(β,α,δ,γvals)




#Task 4, Markov 

using LinearAlgebra

function σ(x,z)
    z==1 ? 0 : z==2 ? x : x<5 ? x+1 : 3
end

function build_transition_matrix(P)
    s=[]
    for x in 0:5, z in 1:3
        push!(s,(x,z))
    end
    T=zeros(18,18)
    for i in 1:18
        (x,z)=s[i]
        x′=σ(x,z)
        for z′ in 1:3
            j=findall(u->u==(x′,z′),s)[1]
            T[i,j]=P[z,z′]
        end
    end
    s,T
end

function stationary_distribution(T)
    λ,V=eigen(transpose(T))
    i=argmin(abs.(λ.-1))
    p=real(V[:,i])
    p./sum(p)
end

function marginalX(p,s)
    m=zeros(6)
    for i in 1:length(p)
        (x,z)=s[i]
        m[x+1]+=p[i]
    end
    m
end

function expectedX(m)
    sum(x*m[x+1] for x in 0:5)
end

function markov_analysis(P)
    s,T=build_transition_matrix(P)
    p=stationary_distribution(T)
    m=marginalX(p,s)
    e=expectedX(m)
    T,p,m,e
end

P=[0.5 0.3 0.2;0.2 0.7 0.1;0.3 0.3 0.4]
T,p,m,e=markov_analysis(P)





    
