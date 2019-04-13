#############################################################################
# Date created: 2018/10/06
# Description: Lessons learned the hard way.
#############################################################################
# Do not redefine functions
function foo1()
    x = 1
    if x != 1
        error("Wrong")
    end
    x = 2
end
foo1()

function foo()
    function p(t) return t + 1 end
    if p(1) != 2
        error("Wrong")
    end
    function p(t) return 1 end
end
foo() # throws the error

# Can use anonymous function instead
# https://stackoverflow.com/questions/52680527/redefining-functions-in-julia-function-gives-strange-behaviour/52680696#52680696
function foo2()
    p = t -> t + 1
    if p(1) != 2
        error("Wrong")
    end
    p = t -> 1
end
foo2()