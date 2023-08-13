function _build_weightmx(
    ℒ, data::AbstractVector{D}, adjl::Vector{Vector{T}}, basis::B
) where {D<:AbstractArray,T<:Int,B<:AbstractRadialBasis}
    TD = eltype(first(data))
    dim = length(first(data)) # dimension of data
    nmon = binomial(dim + basis.poly_deg, basis.poly_deg)
    k = length(first(adjl))  # number of data in influence/support domain
    sizes = (k, nmon)

    # build monomial basis and operator
    mon = MonomialBasis(dim, basis.poly_deg)
    ℒmon = ℒ(mon)
    ℒrbf = ℒ(basis)

    # allocate arrays to build sparse matrix
    I = zeros(Int, k * length(adjl))
    J = reduce(vcat, adjl)
    V = zeros(TD, k * length(adjl))

    n_threads = Threads.nthreads()

    # create work arrays
    n = sum(sizes)
    A = Symmetric[Symmetric(zeros(TD, n, n), :U) for _ in 1:n_threads]
    b = Vector[zeros(TD, n) for _ in 1:n_threads]
    range = Vector{UnitRange{<:Int}}(undef, n_threads)
    d = Vector{Vector{eltype(data)}}(undef, n_threads)

    # build stencil for each data point and store in global weight matrix
    Threads.@threads for i in eachindex(adjl)
        range[Threads.threadid()] = ((i - 1) * k + 1):(i * k)
        @turbo I[range[Threads.threadid()]] .= i
        d[Threads.threadid()] = data[adjl[i]]
        V[range[Threads.threadid()]] = @views _build_stencil!(
            A, b, Threads.threadid(), ℒrbf, ℒmon, d, basis, mon, k
        )
    end

    return sparse(I, J, V, length(adjl), length(data))
end

function _build_weight_vec(
    ℒ, data::AbstractVector{D}, adjl::Vector{Vector{T}}, basis::B
) where {D<:AbstractArray,T<:Int,B<:AbstractRadialBasis}
    TD = eltype(first(data))
    dim = length(first(data)) # dimension of data
    nmon = binomial(dim + basis.poly_deg, basis.poly_deg)
    k = length(first(adjl))  # number of data in influence/support domain
    sizes = (k, nmon)

    # build monomial basis and operator
    mon = MonomialBasis(dim, basis.poly_deg)
    ℒmon = ℒ(mon)
    ℒrbf = ℒ(basis)

    # allocate arrays to build sparse matrix
    V = [zeros(TD, k) for _ in 1:length(adjl)]

    n_threads = Threads.nthreads()

    # create work arrays
    n = sum(sizes)
    A = Symmetric[Symmetric(zeros(TD, n, n), :U) for _ in 1:n_threads]
    b = Vector[zeros(TD, n) for _ in 1:n_threads]
    d = Vector{Vector{eltype(data)}}(undef, n_threads)

    # build stencil for each data point and store in global weight matrix
    Threads.@threads for i in eachindex(adjl)
        d[Threads.threadid()] = data[adjl[i]]
        V[i] = @views _build_stencil!(
            A, b, Threads.threadid(), ℒrbf, ℒmon, d, basis, mon, k
        )
    end

    return V
end

function _build_stencil!(
    A::Vector{<:Symmetric},
    b::Vector{<:Vector},
    id::Int,
    ℒrbf,
    ℒmon,
    data::AbstractVector{D},
    basis::B,
    mon::MonomialBasis,
    k::Int,
) where {D<:AbstractArray,B<:AbstractRadialBasis}
    _build_collocation_matrix!(A[id], data[id], basis, mon, k)
    _build_rhs!(b[id], ℒrbf, ℒmon, data[id], basis, k)
    return (A[id] \ b[id])[1:k]
end

function _build_collocation_matrix!(
    A::Symmetric, data::AbstractVector{D}, basis::B, mon::MonomialBasis, k::K
) where {D<:AbstractArray,B<:AbstractRadialBasis,K<:Int}
    # radial basis section
    @inbounds for j in 1:k, i in 1:j
        parent(A)[i, j] = basis(data[i], data[j])
    end

    # monomial augmentation
    if basis.poly_deg > -1
        @inbounds for i in 1:k
            parent(A)[i, (k + 1):end] .= mon(data[i])
        end
    end

    return nothing
end

function _build_rhs!(
    b::AbstractVector, ℒrbf, ℒmon, data::AbstractVector{D}, basis::B, k::K
) where {D<:AbstractArray,B<:AbstractRadialBasis,K<:Int}
    # radial basis section
    @inbounds for i in eachindex(data)
        b[i] = ℒrbf(first(data), data[i])
    end

    # monomial augmentation
    if basis.poly_deg > -1
        N = length(b)
        bmono = view(b, (k + 1):N)
        ℒmon(bmono, first(data))
    end

    return nothing
end

#=

function _build_weightmx(
    ℒ, data::AbstractVector{D}, adjl::Vector{Vector{T}}, basis::B
) where {D<:AbstractArray,T<:Int,B<:AbstractRadialBasis}
    TD = eltype(first(data))
    dim = length(first(data)) # dimension of data
    nmon = binomial(dim + basis.poly_deg, basis.poly_deg)
    k = length(first(adjl))  # number of data in influence/support domain
    sizes = (k, nmon)

    # build monomial basis and operator
    mon = MonomialBasis(dim, basis.poly_deg)
    ℒmon = ℒ(mon)
    #mon = nothing
    #ℒmon = nothing
    ℒrbf = ℒ(basis)

    # allocate arrays to build sparse matrix
    I = zeros(Int, k * length(data))
    J = reduce(vcat, adjl)
    V = zeros(TD, k * length(data))

    # create work arrays
    n = sum(sizes)
    A = Symmetric(zeros(TD, n, n), :U)
    b = zeros(TD, n)

    # build stencil for each data point and store in global weight matrix
    for i in eachindex(adjl)
        range = ((i - 1) * k + 1):(i * k)
        @turbo I[range] .= i
        #scale = sum([euclidean(d, data[adjl[i][1]]) for d in data[adjl[i]]]) / k
        #d = data[adjl[i]] / scale
        d = data[adjl[i]]
        V[range] = @views _build_stencil!(A, b, ℒrbf, ℒmon, d, basis, mon, k)
    end
    return sparse(I, J, V)
end

function _build_stencil!(
    A::Symmetric, b::Vector, ℒrbf, ℒmon, data::AbstractVector{D}, basis::B, mon, k::Int
) where {D<:AbstractArray,B<:AbstractRadialBasis}
    _build_collocation_matrix!(A, data, basis, mon, k)
    _build_rhs!(b, ℒrbf, ℒmon, data, basis, k)
    display(A)
    display(cond(A))
    return (A \ b)[1:k]
end

function _build_collocation_matrix!(
    A::Symmetric, data::AbstractVector{D}, basis::B, mon, k::K
) where {D<:AbstractArray,B<:AbstractRadialBasis,K<:Int}
    # radial basis section
    @inbounds for j in 1:k, i in 1:j
        parent(A)[i, j] = basis(data[i], data[j])
    end

    # monomial augmentation
    if basis.poly_deg > -1
        @inbounds for i in 1:k
            parent(A)[i, (k + 1):end] .= mon(data[i])
        end
    end

    return nothing
end

function _build_rhs!(
    b::AbstractVector, ℒrbf, ℒmon, data::AbstractVector{D}, basis::B, k::K
) where {D<:AbstractArray,B<:AbstractRadialBasis,K<:Int}
    # radial basis section
    @inbounds for i in eachindex(data)
        b[i] = ℒrbf(first(data), data[i])
    end

    # monomial augmentation
    if basis.poly_deg > -1
        N = length(b)
        bmono = view(b, (k + 1):N)
        ℒmon(bmono, first(data))
    end

    return nothing
end
=#
