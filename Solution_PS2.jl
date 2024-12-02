# Problem 1

function iterative_solver(f, x0; α=0.0, ε=1e-6, maxiter=1000)
    # Initialize variables
    x_n = x0
    iterations = [x0]
    residuals = []
    
    # Iterate
    for i in 1:maxiter
        g_x = f(x_n) + x_n
        x_next = (1 - α) * g_x + α * x_n
        diff = abs(x_next - x_n)
        
        # Store residual and update current value
        push!(residuals, diff)
        push!(iterations, x_next)
        
        # Check for convergence
        if diff / (1 + abs(x_next)) < ε
            return 0, x_next, f(x_next), diff, iterations, residuals
        end
        
        # Update x_n for the next iteration
        x_n = x_next
    end
    
    # Return NaN if no solution found
    return 1, NaN, NaN, NaN, iterations, residuals
end

# Define the function f(x)
f(x) = x^3 - x - 1

# Test the function
flag, solution, f_value, diff, iterations, residuals = iterative_solver(f, 1.0; α=0.5)

# Print the results
println("Flag (0 if solution found, 1 otherwise): $flag")
println("Solution: $solution")
println("f(solution): $f_value")
println("Difference between iterations: $diff")
println("All iterations: $iterations")
println("All residuals: $residuals")


# Problem 2 

using LinearAlgebra

# Function to compute the exact solution
function exact_solution(α, β)
    x5 = 1  # From the last row of the matrix
    x4 = x5 + 0
    x3 = x4 + 0
    x2 = x3 + 0
    x1 = α - β * x5 + β
    return [x1, x2, x3, x4, x5]
end

# Function to compute all required results
function solve_system(α, β)
    # Define A and b
    A = [1 -1  0  α-β   β; 
         0  1 -1     0   0; 
         0  0  1    -1   0; 
         0  0  0     1  -1; 
         0  0  0     0   1]
    b = [α, 0, 0, 0, 1]
    
    # Compute exact solution
    x_exact = exact_solution(α, β)
    
    # Compute numerical solution using \
    x_numerical = A \ b
    
    # Compute relative residual
    residual = norm(A * x_numerical - b) / norm(b)
    
    # Compute condition number
    cond_number = cond(A)
    
    return x_exact[1], x_numerical[1], residual, cond_number
end

# Create the table for α = 0.1 and varying β
function create_table()
    α = 0.1
    β_values = [10.0^i for i in 0:12]
    results = []
    
    for β in β_values
        x1_exact, x1_num, residual, cond_num = solve_system(α, β)
        push!(results, (β, x1_exact, x1_num, residual, cond_num))
    end
    
    return results
end

# Generate and display the table
results = create_table()

println("Table of Results:")
println("β\tExact x1\tNumerical x1\tResidual\tCondition Number")
for row in results
    println(row)
end

# Problem 3 

using Roots 

# Function to calculate Net Present Value (NPV)
function NPV(r, C)
    T = length(C) - 1
    npv = sum(C[t + 1] / (1 + r)^t for t in 0:T)
    return npv
end

# Function to calculate the internal rate of return (IRR)
function internal_rate(C)
    # Check for valid input
    if all(c >= 0 for c in C) || all(c <= 0 for c in C)
        return "Warning: No IRR exists as all cash flows have the same sign."
    end
    
    # Define the function whose root we need to find
    f(r) = NPV(r, C)
    
    # Use a root-finding solver to find the IRR
    try
        root = find_zero(f, (0.0, 1.0))  # Search for IRR in the range (0, 1)
        return root
    catch e
        return "Warning: Unable to find a root for the IRR."
    end
end

# Example test
C = [-5, 0, 0, 2.5, 5]
irr = internal_rate(C)
println("The IRR is: ", irr)

# Problem 4 

using JuMP
using Ipopt

function cost_minimization(α, σ, w1, w2, y)
    # Create a model with the Ipopt optimizer
    model = Model(Ipopt.Optimizer)
    
    # Define decision variables x1 and x2
    @variable(model, x1 >= 0)
    @variable(model, x2 >= 0)
    
    # Define the objective function (minimizing cost)
    @objective(model, Min, w1 * x1 + w2 * x2)
    
    # Add the production constraint
    if σ == 1
        # Cobb-Douglas special case
        @constraint(model, x1^α * x2^(1 - α) == y)
    else
        # General CES production function
        @constraint(model, (α * x1^((σ - 1) / σ) + (1 - α) * x2^((σ - 1) / σ))^(σ / (σ - 1)) == y)
    end
    
    # Solve the optimization problem
    optimize!(model)
    
    # Return the optimal values of x1, x2, and the minimum cost
    return (value(x1), value(x2), objective_value(model))
end


using Plots

function plot_cost_and_inputs(α, σ_vals, w1_vals, w2, y)
    for σ in σ_vals
        costs = []
        x1_vals = []
        x2_vals = []
        
        for w1 in w1_vals
            x1, x2, cost = cost_minimization(α, σ, w1, w2, y)
            push!(costs, cost)
            push!(x1_vals, x1)
            push!(x2_vals, x2)
        end
        
        # Plot results for the current σ
        plot(w1_vals, costs, label="Cost (σ=$σ)", xlabel="w1", ylabel="Value", legend=:top)
        plot!(w1_vals, x1_vals, label="x1 (σ=$σ)")
        plot!(w1_vals, x2_vals, label="x2 (σ=$σ)")
    end
end


α = 0.5
σ_vals = [0.25, 1, 4]
w1_vals = 0.1:0.1:10
w2 = 1
y = 1

# Plot cost and input demand functions
plot_cost_and_inputs(α, σ_vals, w1_vals, w2, y)
