M = [1 2; 3 5]
print(rank(M))

N = [2 3; 4 6]
print(rank(N))

M_N = hcat(M,N)
print(rank(M_N))

A = [2 3 1 9; 3 4 5 6; 4 5 6 7; 3 8 7 5]
print(rank(A))

combined = vcat(M_N, A)
print(rank(combined))

print(rank(M_N))
print(rank(vcat(M_N, A[1:1,:])))
print(rank(vcat(M_N, A[1:2,:])))
print(rank(vcat(M_N, A[1:3,:])))
print(rank(vcat(M_N, A[1:2,:], A[4:4,:])))

using RowEchelon
rref(A)
rref(M_N)
rref(combined')