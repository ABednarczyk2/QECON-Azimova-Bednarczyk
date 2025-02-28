
# hello world 

using Plots, Distributions, Statistics

# 1. Model Parameters & Calibration

β    = 0.96             # discount factor
γ    = 2.0              # CRRA coefficient
ρ    = 0.9              # persistence of productivity
σ    = 0.4              # std dev of productivity shock
ϕ    = 0.0              # borrowing limit
labor_share = 0.6291631 # from Penn World Table
α    = 1.0 - labor_share  # capital share in production
println("alpha = $α")

# We want the baseline with r=0.04, w=1, G/Y=0.2, I/Y=0.2
r_target = 0.04
w_target = 1.0

# Derive K, δ, A from conditions in the project:
# 1) I/Y = δK / Y = 0.2  => δK = 0.2Y
# 2) r = α A K^(α-1)*L^(1-α) - δ  (with L=1)
# 3) w = (1-α)*A*K^α
# 4) w=1 => A*K^α = 1/(1-α)
# Let Y = A*K^α * L^(1-α) = A*K^α (since L=1).
# Then G/Y=0.2 => T = 0.2Y (in baseline, λ=0 => T = τ*y => τ=0.2)
Y_target = 1.0/(1-α)               # from w=1 => A*K^α=1/(1-α)
K_guess  = ((α - 0.2)*Y_target)/r_target
δ_guess  = 0.2*Y_target / K_guess
A_guess  = Y_target/(K_guess^α)

println("Initial guesses:")
println("   K = $(round(K_guess,3)), δ = $(round(δ_guess,3)), A = $(round(A_guess,3))")

τ_baseline = 0.2
λ_baseline = 0.0
λ_reform   = 0.15   


# 2. Grids & Tauchen Discretization

a_min   = 0.0
a_max   = 10.0
num_a   = 100
a_grid  = range(a_min, a_max, length=num_a)

num_z   = 5

function tauchen(n, ρ, σ, m=3)
    s = σ/sqrt(1-ρ^2)
    z_max = m*s; z_min = -z_max
    z_points = range(z_min, z_max, length=n)
    step = z_points[2] - z_points[1]
    P = zeros(n,n)
    for i in 1:n
        for j in 1:n
            midpoint_low  = (z_points[j] - ρ*z_points[i] - step/2)/σ
            midpoint_high = (z_points[j] - ρ*z_points[i] + step/2)/σ
            if j == 1
                P[i,j] = cdf(Normal(0,σ), midpoint_high)
            elseif j == n
                P[i,j] = 1.0 - cdf(Normal(0,σ), midpoint_low)
            else
                P[i,j] = cdf(Normal(0,σ), midpoint_high) - cdf(Normal(0,σ), midpoint_low)
            end
        end
    end
    z_shifted = z_points .- mean(z_points)
    z_levels = exp.(z_shifted)
    # we scale so that average is ~1
    mean_z = sum(z_levels)/n
    z_levels ./= mean_z
    return z_levels, P
end

z_grid, P_z = tauchen(num_z, ρ, σ)


# 3. Utility & Tax Functions
function utility(c)
    return c>0 ? (c^(1-γ) - 1)/(1-γ) : -1e10
end

# T(y) = y - (1-τ)*(y/ȳ)^(1-λ)* ȳ, but we’ll treat ȳ=1 for simplicity
function tax(y, τ, λ)
    return y - (1-τ)* (y^(1-λ))
end

function after_tax(y, τ, λ)
    return (1-τ)* (y^(1-λ))
end


# 4. Household VFI

function solve_household(λ, τ; max_iter=500, tol=1e-6)
    # Value function array
    V = zeros(num_a, num_z)
    policy_idx = fill(1, (num_a, num_z))
    
    for it in 1:max_iter
        V_new = similar(V)
        max_diff = 0.0
        
        for iz in 1:num_z
            z = z_grid[iz]
            # Pre-tax labor income with w_target:
            y = z * w_target
            # After-tax labor income piece => y - tax(y)
            # But we'll compute inside the loop
            for ia in 1:num_a
                a_now = a_grid[ia]
                
                best_val  = -1e10
                best_aidx = 1
                
                for ja in 1:num_a
                    a_next = a_grid[ja]
                    c = (y - tax(y, τ, λ)) + (1 + r_target)*a_now - a_next
                    u = utility(c)
                    
                    # expected future
                    EV = 0.0
                    for izp in 1:num_z
                        EV += P_z[iz, izp]*V[ja, izp]
                    end
                    val = u + β*EV
                    if val>best_val
                        best_val = val
                        best_aidx = ja
                    end
                end
                
                V_new[ia, iz] = best_val
                diff_local = abs(V_new[ia, iz] - V[ia, iz])
                if diff_local>max_diff
                    max_diff = diff_local
                end
                policy_idx[ia, iz] = best_aidx
            end
        end
        
        V .= V_new
        if it % 10 == 0
            println("VFI iteration $it, error = $max_diff")
        end
        if max_diff<tol
            println("VFI converged at iteration $it, error = $max_diff")
            break
        end
    end
    
    policy = [a_grid[policy_idx[ia, iz]] for ia in 1:num_a, iz in 1:num_z]
    return V, policy, policy_idx
end


# 5. Stationary Distribution

function stationary_dist(policy_idx; tol=1e-8, max_iter=500)
    dist = fill(1.0/(num_a*num_z), (num_a, num_z))
    
    for it in 1:max_iter
        dist_new = zeros(num_a, num_z)
        
        for iz in 1:num_z
            for ia in 1:num_a
                a_next_idx = policy_idx[ia, iz]
                mass = dist[ia, iz]
                for izp in 1:num_z
                    dist_new[a_next_idx, izp] += mass*P_z[iz, izp]
                end
            end
        end
        
        err = maximum(abs.(dist_new .- dist))
        if it % 10 == 0
            println("Dist iteration $it, error = $err")
        end
        dist .= dist_new
        if err<tol
            println("Distribution converged at iteration $it, error=$err")
            return dist
        end
    end
    
    println("WARNING: Dist did not converge fully after $max_iter iterations.")
    return dist
end


# 6. Outer Equilibrium Solve

# We hold r_target, w_target fixed in the VFI. Then we try to find τ so that
# the ratio of government revenue to output is 0.2. We'll do a small dampening
# to avoid overshoot.

function solve_equilibrium(λ; max_iter=50, tol=1e-4)
    # Start from baseline guess
    τ = τ_baseline
    θ = 0.2   # dampening factor for updating τ
    
    for it in 1:max_iter
        #  household problem at current τ
        println("\n--- Outer iteration $it, current τ=$τ ---")
        V, policy, p_idx = solve_household(λ, τ)
        dist = stationary_dist(p_idx)
        
        # Aggregate assets:
        aggA = sum(dist .* repeat(collect(a_grid),1,num_z))
        # Aggregate labor:
        L    = sum([sum(dist[:,iz])*z_grid[iz] for iz in 1:num_z])
        # Production:
        Y    = A_guess*aggA^α * L^(1-α)
        
        # Government revenue:
        revenue = 0.0
        for iz in 1:num_z
            mass_z = sum(dist[:,iz])
            y_iz   = z_grid[iz]*w_target
            revenue += tax(y_iz, τ, λ)*mass_z
        end
        rev_ratio = revenue/Y
        err = abs(rev_ratio - 0.2)
        println("   Government rev_ratio = $rev_ratio, error=$err")
        
        τ_new = τ + θ*(0.2 - rev_ratio)*τ  # scale the tax by (0.2 - rev_ratio)
        
        if err<tol
            println("Equilibrium in tax rate found at iteration $it, τ=$τ_new, error=$err")
            return (τ=τ_new, V=V, policy=policy, policy_idx=p_idx, dist=dist,
                    K=aggA, L=L, Y=Y, rev_ratio=rev_ratio)
        end
        
        τ = 0.5*τ + 0.5*τ_new  
    end
    
    println("WARNING: Equilibrium loop ended without hitting tolerance.")
    return (τ=τ, V=zeros(num_a,num_z), policy=zeros(num_a,num_z),
            policy_idx=zeros(Int,num_a,num_z), dist=zeros(num_a,num_z),
            K=0.0, L=0.0, Y=0.0, rev_ratio=999.0)
end


# 7. Gini & Lorenz

function gini(x, w)
    order = sortperm(x)
    x_s   = x[order]
    w_s   = w[order]
    cumw  = cumsum(w_s)
    cumxw = cumsum(x_s .* w_s)
    B     = sum(w_s .* cumxw)/(sum(w_s)*sum(x_s.*w_s))
    return 1.0 - 2.0*B
end

function lorenz_curve(x, w)
    order = sortperm(x)
    x_s   = x[order]
    w_s   = w[order]
    cum_pop  = cumsum(w_s)./sum(w_s)
    cum_xval = cumsum(x_s.*w_s)./sum(x_s.*w_s)
    return cum_pop, cum_xval
end

# 8. Baseline (λ=0) and Reform (λ=0.15)

println("========== Solving Baseline (lambda=0) ==========")
eq_base = solve_equilibrium(λ_baseline)

println("\n========== Solving Reform (lambda=0.15) ==========")
eq_reform = solve_equilibrium(λ_reform)

println("\n--- Baseline Results ---")
println("   τ = ", eq_base.τ)
println("   K = ", eq_base.K, ", Y = ", eq_base.Y, ", K/Y=", eq_base.K/eq_base.Y)
println("   rev_ratio = ", eq_base.rev_ratio)

println("\n--- Reform Results ---")
println("   τ = ", eq_reform.τ)
println("   K = ", eq_reform.K, ", Y = ", eq_reform.Y, ", K/Y=", eq_reform.K/eq_reform.Y)
println("   rev_ratio = ", eq_reform.rev_ratio)

# 9. Gini Coefficients

function compute_ginis(eq, λ)
    dist = eq.dist
    # Flatten
    nstates = num_a*num_z
    inc_vec = similar(dist)
    ast_vec = similar(dist)
    idx = 1
    for iz in 1:num_z
        for ia in 1:num_a
            inc_vec[ia, iz] = after_tax(z_grid[iz]*w_target, eq.τ, λ)
            ast_vec[ia, iz] = a_grid[ia]
        end
    end
    wghts = dist[:]
    incv  = inc_vec[:]
    astv  = ast_vec[:]
    gini_inc = gini(incv, wghts)
    gini_ast = gini(astv, wghts)
    return (gini_income=gini_inc, gini_assets=gini_ast)
end

gini_base   = compute_ginis(eq_base, λ_baseline)
gini_reform = compute_ginis(eq_reform, λ_reform)

println("\n--- Baseline Gini ---")
println("   After-tax labor income: ", gini_base.gini_income)
println("   Assets: ", gini_base.gini_assets)

println("\n--- Reform Gini ---")
println("   After-tax labor income: ", gini_reform.gini_income)
println("   Assets: ", gini_reform.gini_assets)

# I lost the initisl work 














