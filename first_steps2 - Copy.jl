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



matrix_data = readdlm("C:\\Users\\2022\\Desktop\\Stuff\\data.csv", ',')


