using Images, Colors, ImageIO, ImageMagick

imgSize = 400
r = imgSize/2
runs = 25

function xtorgb(x)

    #Three color values are set to 0 by default
    red = 0
    blue = 0
    green = 0

    #Determines the green color
    if x <= -1
        green = 255
    elseif x > -1 && x < 0
        green = -x*255
    elseif x >= 0 && x <= 2
        green = 0
    elseif x > 2 && x <= 3
        green = (x-2)*255
    elseif x > 3
        green = 255
    end

    #Determines the red color
    if x <= -2
        red = 0
    elseif x > -2 && x < -1
        red = (2+x)*255
    elseif x >= -1 && x <= 1
        red = 255
    elseif x > 1 && x < 2
        red = (2-x)*255
    elseif x >= 2
        red = 0
    end

    #Determines the blue color
    if x <= -2
        blue = (-x-2)*255
    elseif x >= -2 && x <= 0
        blue = 0
    elseif x > 0 && x < 1
        blue = x*255
    elseif x >= 1
        blue = 255
    end

    return [red; green; blue]
end

function picture(plane)
    resolutioni = size(plane)[1]
    resolutionj = size(plane)[2]
    white = [255.0; 255.0; 255.0]
    p = zeros(3, resolutionj, resolutioni)

    for j1 = 1:resolutionj
        line = zeros(3, resolutioni)

        for i1 = 1:resolutioni
            i = real.(plane[i1,j1])
            j = imag.(plane[i1,j1])
            a = xtorgb(3*atan(i,j)/pi)
            d = sqrt(i^2 + j^2)

            if d < 1
                c = white*(1-d)+d*a
            else
                c = a/d
            end

            for l = 1:3
                if c[l] < 0.0
                    c[l] = 0.0
                elseif c[l] > 255.0
                    c[l] = 255.0
                end
            end

            line[:,i1] = c
        end

        p[:,:,j1] = line
    end

    return p
end

function f(x)
    return x^4 + 2x - 1
end

function df(x)
    return 4x^3 + 2
end

function ddf(x)
    return 12x^2
end

function cubicmethod(x)
    return x - (f(x)/df(x))*(1 + (f(x)*ddf(x))/(2*df(x)^2))
end

function createFractal(cx, cy, fw, fn, numruns)
    tmpa = transpose(fw*collect(-r:r-1)/r);
    tmpb = ones(1,imgSize);
    plane = ((transpose(tmpa)*tmpb).+cy)*1im;
    plane = plane + ((transpose(tmpb)*tmpa).+cx);
    #pictureArray = zeros(imgSize, imgSize, numruns)*0im;
    #newPlane = zeros(imgSize, imgSize)*0im;

    colorwheel = copy(plane)

    #save("Cubic Fractals\\$(fn)cubic0.png", map(clamp01nan,colorview(RGB, picture(plane)/255)))

    #global plane, newPlane, pictureArray
    resolutioni = size(plane)[1]
    resolutionj = size(plane)[2]

    for i1 = 1:resolutioni
        for j1 = 1:resolutionj
            for i = 1:numruns
                plane[i1,j1] = cubicmethod(plane[i1,j1])
            end
        end
    end

    #pictureArray[:,:,i] = plane

    save("Cubic Fractals\\cubic$(fn).png", map(clamp01nan,colorview(RGB, picture(plane)/255)))
end

#save("Cubic Fractals\\cubicfractal.png", colorview(RGB, picture(pictureArray[:,:,runs]) / 255))
#{{x -> -1.39534}, {x -> 0.460355 - 1.13932 I}, {x -> 0.460355 + 1.13932 I}, {x -> 0.474627}}

#createFractal(-1.739674185463659, 1.3742230576441103, 0.020025062656641612, "z")
#createFractal(0.0, 0.0, 2.0, "a")

wa = [0.0 0.0 2.0]
wz = [-1.739674185463659 1.3742230576441103 0.020025062656641612]

nmax = 100
s = 30
oldsize = 1
alpha = (1/nmax)*log(wz[3]/wa[3])

for n = 0:nmax
    if n < 51
        numruns = 10
    else
        numruns = 20
    end

    theta = (wa[3]/(wz[3]-wa[3]))*(exp(alpha*n)-1)
    wn = (1-theta)*wa+theta*wz

    #println(oldsize/wn[3], " ", temptest, " ", wn)
    #global oldsize = wn[3]
    createFractal(wn[1], wn[2], wn[3], n, numruns)
end

#m(62) = 1.3742230576441103
#l(26) = -1.739674185463659
#width = 0.020025062656641612
