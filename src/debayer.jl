export debayer_bilinear_cpu


"""
    debayer_bilinear_cpu(img; pattern=:RGGB) -> Matrix{RGB{N0f16}}

Bilinear demosaicing for Bayer CFAs. Supported patterns: :RGGB, :BGGR, :GRBG, :GBRG.
"""
function debayer_bilinear_cpu(img::AbstractMatrix{<:Unsigned}; pattern::Symbol = :RGGB)
    h, w = size(img)
    R = zeros(UInt16, h, w)
    G = zeros(UInt16, h, w)
    B = zeros(UInt16, h, w)

    # parity of (row_is_odd, col_is_odd) for R and B anchor sites
    rpar, bpar = _rb_parity(pattern)

    @inline function is_R_site(y::Int, x::Int)
        (isodd(y) == rpar[1]) & (isodd(x) == rpar[2])
    end
    @inline function is_B_site(y::Int, x::Int)
        (isodd(y) == bpar[1]) & (isodd(x) == bpar[2])
    end
    @inline function is_G_site(y::Int, x::Int)
        !is_R_site(y,x) & !is_B_site(y,x)
    end
    # Distinguish the two green positions: G on "R row" vs "B row"
    @inline function is_G_on_Rrow(y::Int, x::Int)
        is_G_site(y,x) & (isodd(y) == rpar[1])
    end
    @inline function is_G_on_Brow(y::Int, x::Int)
        is_G_site(y,x) & (isodd(y) == bpar[1])
    end

    # Skip 1-pixel border (uses 4-neighborhood/diagonals)
    for y in 2:h-1, x in 2:w-1
        if is_R_site(y,x)
            R[y,x] = img[y,x]
            G[y,x] = sum(UInt32(img[y+i, x+j]) for (i,j) in ((-1,0),(1,0),(0,-1),(0,1))) ÷ 4
            B[y,x] = sum(UInt32(img[y+i, x+j]) for (i,j) in ((-1,-1),(-1,1),(1,-1),(1,1))) ÷ 4

        elseif is_B_site(y,x)
            B[y,x] = img[y,x]
            G[y,x] = sum(UInt32(img[y+i, x+j]) for (i,j) in ((-1,0),(1,0),(0,-1),(0,1))) ÷ 4
            R[y,x] = sum(UInt32(img[y+i, x+j]) for (i,j) in ((-1,-1),(-1,1),(1,-1),(1,1))) ÷ 4

        elseif is_G_on_Rrow(y,x)
            # Green located on the same row-parity as R anchors:
            # R from left/right (horizontal), B from up/down (vertical)
            G[y,x] = img[y,x]
            R[y,x] = (UInt32(img[y, x-1]) + UInt32(img[y, x+1])) ÷ 2
            B[y,x] = (UInt32(img[y-1, x]) + UInt32(img[y+1, x])) ÷ 2

        elseif is_G_on_Brow(y,x)
            # Green located on the same row-parity as B anchors:
            # R from up/down (vertical), B from left/right (horizontal)
            G[y,x] = img[y,x]
            R[y,x] = (UInt32(img[y-1, x]) + UInt32(img[y+1, x])) ÷ 2
            B[y,x] = (UInt32(img[y, x-1]) + UInt32(img[y, x+1])) ÷ 2
        end
    end

    norm = typemax(UInt16)
    RGB_image = [RGB{N0f16}(R[y,x]/norm, G[y,x]/norm, B[y,x]/norm) for y in 1:h, x in 1:w]
    return colorview(RGB{N0f16}, RGB_image)
end
