function declu(A, ϵ = 1e-12)
    n = size(A, 1)
    for j = 1:n
        if abs(A[j,j]) < ϵ
            error("Matriz muito próxima de ser singular")
        end
        I = j+1:n
        @views A[I,j] ./= A[j,j]
        @views A[I,I] .-= A[I,j] * A[j,I]'
    end
    return A
end

function reslu(A, b)
    n = length(b)

    for i = 1:n
        for j = 1:i-1
            b[i] += - A[i,j] * b[j]
        end
    end

    for i = n:-1:1
        for j = i+1:n
            b[i] += - A[i,j] * b[j]
        end
        b[i] *= 1 / A[i,i]
    end

    return b
end

function reslu2(A, b)
    n = length(b)
    for i = 1:n
        for j = 1:i-1
            b[i] += - A[i,j] * b[j]
        end
    end

    for j = n:-1:1
        b[j] *= 1 / A[j,j]
        for i = j-1:-1:1
            b[i] += - A[i,j] * b[j]
        end
    end
    return b
end

function reslu3(A, b)
    n = length(b)
    for j = 1:n
        for i = j+1:n
            b[i] += - A[i,j] * b[j]
        end
    end

    for i = n:-1:1
        for j = i+1:n
            b[i] += - A[i,j] * b[j]
        end
        b[i] *= 1 / A[i,i]
    end
    return b
end

function reslu4(A, b)
    n = length(b)
    for j = 1:n
        for i = j+1:n
            b[i] *= - A[i,j] * b[j]
        end
    end

    for j = n:-1:1
        b[j] *= 1 / A[j,j]
        for i = j-1:-1:1
            b[i] += - A[i,j] * b[j]
        end
    end
    return b
end
