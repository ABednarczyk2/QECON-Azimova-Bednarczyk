using Plots



function odd_or_even(n)
    return n % 2==0 ? "even" : "odd"
end
odd_or_even(137336146183)



compare_three(1, 2, 3)
function compare_three(a, b, c)
        if a > 0 && b > 0 && c > 0
            println("All numbers are positive")
        elseif a == 0 && b == 0 && c == 0
            println("All numbers are zero")
        else
            println("At least one number is not positive")
        end
    end
    
 compare_three(0, 0, -1)  

    variable_result=1
    
    """
    my_factorial(n::Int)

TBW
"""
function my_factorial(n::Int)
    result=1
    for i in 1:n
        result*=i
end
return result
end
my_factorial(2)
my_factorial(16)
my_factorial(2)
my_factorial(20)

function count_positives(arr::AbstractVector{<:Number})
return count(x->x>0,arr)
end

count_positives([1,2,3])
count_positives([1, -3, 4, 7, -2, 0])

rand(3,5)


using DelimitedFiles
using Plots
using Statistics

using Plots

function plot_powers(n)
    x = -10:0.2:10  
    plot(title="Powers of x", xlabel="x", ylabel="y")

    for i in 1:n
        y = x .^ i   
        plot!(x, y, label="x^$i")  
    end

    display(plot())  
end

plot_powers(50)  

function plot_powers(n)
    x = -10:0.2:10  
    plot2 = plot(title="Powers of x", xlabel="x", ylabel="y") 

    for i in 1:n
    
        power_function = x -> x^i

        
        y = power_function.(x)

        
        plot!(x, y, label="x^$i", linewidth=3, linestyle=:dash)
    end

    return plot2  
end

my_plot = plot_powers(3) 

function count_positives_broadcasting(arr)
   
    positives = arr .>0
    
   
    count = sum(positives)
    
    return count
end

arr = [1, 92, -1,-1,-0.1,-3, 4, 7, -2, 0]
count_positives2 = count_positives_broadcasting(arr)

function standard_deviation(x)
    n = length(x)                      
    mean_x = sum(x) / n                
    squared_diffs = (x .- mean_x) .^ 2 
    variance = sum(squared_diffs) / n  
    return sqrt(variance)              
end

x = [1, 2, 3, 4, 5,6,7,8,9,10]
y=[1,2,3,4,5]
std_dev = standard_deviation(y)
std_dev = standard_deviation(x)
println(std_dev)  


function standard_deviation(x)
    mean_x = sum(x) / length(x)  
    d = x .- mean_x              
    squared_diffs = d .^ 2       
    variance = sum(squared_diffs) / length(x)  
    return sqrt(variance)                      
end
standard_deviation([1, 2, 3, 4, 5])

standard_deviation([5, 10, 15])

Cmatrix_data = readdlm("C:\\Users\\2022\\Desktop\\Stuff\\dataset.csv", ',')

typeof(Cmatrix_data)
size(Cmatrix_data) 
earnings = Cmatrix_data[:, 1]      
education = Cmatrix_data[:, 2]     
hours_worked = Cmatrix_data[:, 3]  








scatter(education, earnings, xlabel="Education", ylabel="Earnings",
        title="Earnings vs Education", markercolor=:green, legend=false)

scatter(hours_worked, earnings, xlabel="Hours Worked", ylabel="Earnings",
        title="Earnings vs Hours Worked",markercolor=:red, legend=false)


        #Diyora's file 

        #= Problem 1: Odd_or_Even 

Write a function odd_or_even(n) that takes an integer n and prints: • "Odd" if the number is odd.
"Even" if the number is even.1

=# 

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

#= function compare_three(a, b, c)
    if a > 0 && b > 0 && c > 0
        println("All numbers are positive")
    elseif 
        println("At least one number is not positive")
    elseif 

    I can't use that ==> ORDER MATTERS
    
        =# 

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
#oh