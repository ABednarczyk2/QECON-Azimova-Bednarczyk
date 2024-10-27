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

Cmatrix_data = readdlm("C:\\Users\\2022\\Desktop\\Stuff\\dataset.csv", ',')


