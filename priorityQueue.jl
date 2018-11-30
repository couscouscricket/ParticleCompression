function queue!(heap::Vector{T}, key::T) where {T}
    push!(heap,key)
    swim!(heap)
end

function dequeue!(heap::Vector{T}) where {T}
    if length(heap) > 1
        last = length(heap)
        heap[1], heap[last] = heap[last], heap[1]
    end
    removedKey = pop!(heap)
    sink!(heap)
    return removedKey
end

up(i::Int) = i >>> 1
left(i::Int) = i << 1
right(i::Int) = (i << 1) + 1

function swim!(heap::Vector{T}) where {T}
    node = length(heap)
    upper = up(node)
    while node > 1 && heap[node] < heap[upper]
        heap[node], heap[upper] = heap[upper], heap[node]
        node = upper
        upper = up(node)
    end
end

function sink!(heap::Vector{T}) where {T}
    node = 1
    lft = left(node)
    rght = right(node)
    last = length(heap)
    while lft <= last
        minChild = lft
        if lft < last && heap[rght] < heap[lft]; minChild = rght end
        if heap[node] < heap[minChild]; break; end
        heap[node], heap[minChild] = heap[minChild], heap[node]
        node = minChild
        lft = left(node)
        rght = right(node)
    end
end
