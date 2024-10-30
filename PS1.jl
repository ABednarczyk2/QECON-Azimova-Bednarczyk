# Problem 1: Odd_or_Even 
#= Using Type Annotation 
(I believe it to be more accurate cause it will identify only integers): =#

function odd_or_even(n :: Int)
    if n % 2 == 0 
        println("Even")
    else 
        println("Odd")
    end 
end 

odd_or_even(3)
odd_or_even(6)
odd_or_even(1)
odd_or_even(0)
odd_or_even(78)

# or without Type Annotation
function odd_or_even(n)
    if n % 2 == 0 
        println("Even")
    else 
        println("Odd")
    end 
end 

odd_or_even(3)
odd_or_even(67)

# Problem 2: Boolean operators 

function compare_three(a, b, c)
    if a == 0 && b == 0 && c == 0
        println("All numbers are zero")
    elseif a > 0 && b > 0 && c > 0 
        println("All numbers are positive")
    else 
        println("At least one number is not positive")
    end 
end 
compare_three(0, 0, 0)
compare_three(5, 7, -2)
compare_three(1000, -1000, 2)
compare_three(2, 3, 5)

# Problem 3: Factorial Calculation Using a Loop

function my_factorial(n :: Int)
    if n < 0 
        error("defined only for non-negative integers")
    end 

    result = 1 
    for i in 1:n 
        result *= i #Multiply result by i in each iteration
    end 
    return result 

end 
my_factorial(5)
my_factorial(-2)

# Problem 4: Count Positive Numbers Using a Loop

function count_positives(arr)
    count = 0 
    for num in arr
        if num > 0 
            count += 1 
        end 
    end 
    println("Output:  $count")
end

count_positives([1, 3, -9, 0, 0, -6])

# Problem 5: Plotting Powers of x Using a Loop

using Plots
function plot_powers(n)
    power_plots = plot(
        xlabel="x", ylabel="y", title="Powers of x",
        linewidth=3, line=:dash
    )

    x_values = -10:0.2:10
    for i in 1:n
        y_values = x_values .^ i
        plot!(power_plots, x_values, y_values, label="x^$i")
    end
    return power_plots
end

my_plot = plot_powers(5)

using Plots
function plot_powers(n)
    x_values = -10:0.2:10
    plot(
        xlim=(-10, 10), xlabel="x", ylabel="y", title="Powers of x",
        linewidth=3, line=:dash, legend=:topright
    )
    for i in 1:n
         y_values = x_values .^ i  
        plot!(x_values, y_values, label="x^$i")
    end
end
plot_powers(5)

# Problem 5: Count Positive Numbers Using Broadcasting

function count_positives_broadcasting(arr)
    positives = arr .> 0 
    count = sum(positives)
    println("Output: $count")
end 
count_positives_broadcasting([1, 3, -9 -100, 50, -2000])

# Problem 6: Standard Deviation Using Broadcasting and Vectorized Operations

    function standard_deviation(x)
        mean_x = sum(x) / length(x)
        d = x .- mean_x
        squared_d = d .^ 2
        variance = sum(squared_d) / (length(x) - 1)
        sd = sqrt(variance)
        return sd
    end
    
    
    println(standard_deviation([1, 2, 3, 4, 5]))  
    println(standard_deviation([5, 10, 15]))      
    println(standard_deviation(2:7))         
    
    # Problem 7: Import the Data, Plot It, and Calculate Correlations
    using Statistics, Plots, DelimitedFiles
    data_1 = readdlm("/Users/diyora/Desktop/Class_materials/Problem_Set_1/dataset.csv", ',',Float64)

plot_1 = scatter(data_1[:,2], data_1[:,1]; legend=false, color=:green, markersize = 5, opacity=0.7)
xaxis!(plot_1, "Education")
yaxis!(plot_1, "Earnings")
title!(plot_1, "Scatter plot of Earnings & Education")
display(plot_1)

plot_2 = scatter(data_1[:,2], data_1[:,3]; legend=false, color=:red, markersize = 5, opacity=0.7)
xaxis!(plot_2, "Hours worked")
yaxis!(plot_2, "Earnings")
title!(plot_2, "Scatter plot of Earnings & Hours Worked")
display(plot_2)

cor(data_1[:,2],data_1[:,1])
cor(data_1[:,2],data_1[:,3])
