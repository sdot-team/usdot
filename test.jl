using Images
using RadonKA

function sino( A, L )
    img = zeros( Float64, L, L )
    img[ 40:90, 50:80 ] .= 1.0
    return radon( img, A )
end

function points( L, N )
    return rand( Float64, N, 2 ) .* L
end

function unidimensionnal_shift( points, density )
    ndi = size( points )[ 1 ] / sum( density )
    nps = sortperm( points )
    num_point = 1
    acd = 0.0
    res = zeros( Float64, size( points ) )
    for ( i, d ) in enumerate( density )
        loc = d * ndi
        while acd + loc >= num_point
            res[ nps[ num_point ] ] = i + ( num_point - acd ) / loc
            num_point += 1
        end

        acd += loc
    end
    if num_point <= size( points )[ 1 ]
        res[ nps[ num_point ] ] = size( density )[ 1 ]
    end
    return res
end


L = 100
A = range( 0, stop = 2*pi, length = L )
s = sino( A, L )

p = points( L, 10 )

us = unidimensionnal_shift( [ 0, 1, 2 ], [ 0.1, 1.0, 0.1 ] )
println( us )

# Display the sinogram
# using ImageView
# imshow(sino(); canvassize=(800, 800))
