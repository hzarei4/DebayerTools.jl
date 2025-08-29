export detect_bit_depth, connected_components_pixels, mean_except

"""
    detect_bit_depth(img::AbstractArray)

Detect the bit depth of the image based on its maximum pixel value.
Returns the bit depth as an integer.
"""
function detect_bit_depth(img::AbstractArray)
    maxval = maximum(img)
    return ceil(Int, log2(maxval + 1))  # e.g. 12 bits if max ~ 4095
end



"""
    cfa_color(y, x; pattern=:RGGB) -> Symbol

Return which CFA color (:R, :G, or :B) is at pixel (y,x)
for the given Bayer pattern.
"""
function cfa_color(y::Int, x::Int; pattern::Symbol=:RGGB)
    yo, xo = isodd(y), isodd(x)
    if pattern == :RGGB
        return (!yo && !xo) ? :R : (yo && xo) ? :B : :G
    elseif pattern == :BGGR
        return (!yo && !xo) ? :B : (yo && xo) ? :R : :G
    elseif pattern == :GRBG
        return (!yo && !xo) ? :G : (!yo && xo) ? :R : (yo && !xo) ? :B : :G
    elseif pattern == :GBRG
        return (!yo && !xo) ? :G : (!yo && xo) ? :B : (yo && !xo) ? :R : :G
    else
        error("Unsupported Bayer pattern: $pattern")
    end
end


"""
    _rb_parity(pattern) -> ((r_rowOdd, r_colOdd), (b_rowOdd, b_colOdd))

Parity (odd=true, even=false) of R and B anchor sites for each Bayer pattern.
"""
function _rb_parity(pattern::Symbol)
    if pattern === :RGGB
        # R at (even,even), B at (odd,odd)
        return ((false,false), (true,true))
    elseif pattern === :BGGR
        # B at (even,even), R at (odd,odd)
        return ((true,true), (false,false))
    elseif pattern === :GRBG
        # R at (even,odd),  B at (odd,even)
        return ((false,true), (true,false))
    elseif pattern === :GBRG
        # R at (odd,even),   B at (even,odd)
        return ((true,false), (false,true))
    else
        error("Unsupported Bayer pattern: $pattern")
    end
end



"""
    connected_components_pixels(img::AbstractMatrix{Bool}; connectivity::Int=8, return_linear::Bool=false)

Find connected components in a binary image and return, for each component,
the full list of pixel indices.

# Arguments
- `img`            : 2D Bool array (true = foreground)
- `connectivity`   : 4 or 8
- `return_linear`  : if true, return linear indices (Vector{Int}) per component
                     instead of CartesianIndex{2}.

# Returns
- `labels::Matrix{Int}` : labeled image (0 for background, 1..N for components)
- `pixels::Vector{Vector{CartesianIndex{2}}}`  if `return_linear == false`
  OR
  `pixels::Vector{Vector{Int}}`                if `return_linear == true`
"""
function connected_components_pixels(img::AbstractMatrix{Bool}; connectivity::Int=8, return_linear::Bool=false)
    @assert connectivity == 4 || connectivity == 8 "connectivity must be 4 or 8"
    nrows, ncols = size(img)
    labels = zeros(Int, nrows, ncols)

    # neighborhood offsets
    neigh = connectivity == 4 ?
        (CartesianIndex(-1,0), CartesianIndex(1,0), CartesianIndex(0,-1), CartesianIndex(0,1)) :
        (CartesianIndex(-1,0), CartesianIndex(1,0), CartesianIndex(0,-1), CartesianIndex(0,1),
         CartesianIndex(-1,-1), CartesianIndex(-1,1), CartesianIndex(1,-1), CartesianIndex(1,1))

    # worst-case stack
    stack = Vector{CartesianIndex{2}}(undef, nrows*ncols)
    top = 0

    comps = return_linear ? Vector{Vector{Int}}() : Vector{Vector{CartesianIndex{2}}}()
    current_label = 0

    @inbounds for r in 1:nrows, c in 1:ncols
        if img[r,c] && labels[r,c] == 0
            current_label += 1
            # start a new component bucket
            push!(comps, return_linear ? Int[] : CartesianIndex{2}[])

            # seed
            top += 1; stack[top] = CartesianIndex(r,c)
            labels[r,c] = current_label

            # DFS flood-fill
            while top > 0
                p = stack[top]; top -= 1
                if return_linear
                    push!(comps[end], LinearIndices(img)[p])
                else
                    push!(comps[end], p)
                end

                for off in neigh
                    q = p + off
                    qr, qc = q.I
                    if 1 <= qr <= nrows && 1 <= qc <= ncols &&
                       img[qr,qc] && labels[qr,qc] == 0
                        labels[qr,qc] = current_label
                        top += 1; stack[top] = q
                    end
                end
            end
        end
    end

    return labels, comps
end


"""
    mean_except(img::AbstractMatrix{<:Number}, skip::AbstractVector{CartesianIndex})

Compute the mean of all values in `img`, ignoring the pixels at indices in `skip`.
"""
function mean_except(img::AbstractMatrix{<:Number}, skip::AbstractVector{CartesianIndex{2}})
    nrows, ncols = size(img)
    total = zero(eltype(img))
    count = 0

    # Build a Set for fast membership checking
    skipset = Set(skip)

    @inbounds for r in 1:nrows, c in 1:ncols
        idx = CartesianIndex(r,c)
        if !(idx in skipset)
            total += img[r, c]
            count += 1
        end
    end
    return count == 0 ? NaN : total / count
end