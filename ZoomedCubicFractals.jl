using Images, Colors, ImageIO, ImageMagick

imgSize = 400
r = imgSize/2
runs = 100

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

function l(x)
    return -2 + x*((1/200)+(2/399))
end

function m(x)
    return -(l(x)+(1/200))
end

tmpa1 = zeros(1,imgSize);
tmpb1 = zeros(1,imgSize);
for i=1:imgSize
    xval = l(25)+(i-1)*((l(27)-l(25))/imgSize);
    yval = m(61)+(i-1)*((m(63)-m(61))/imgSize);
    tmpa1[i] = xval;
    tmpb1[i] = yval;
end

tmpa2 = zeros(imgSize,imgSize);
tmpb2 = zeros(imgSize,imgSize);
for i=1:imgSize
    tmpa2[i,:] = tmpa1;
    tmpb2[:,i] = transpose(tmpb1);
end
tmpb2 = tmpb2*1im

plane = tmpa2 + tmpb2;
pictureArray = zeros(imgSize, imgSize, runs)*0im;
newPlane = zeros(imgSize, imgSize)*0im;

colorwheel = copy(plane)

save("cubicfractal0.png", map(clamp01nan,colorview(RGB, picture(plane)/255)))

for i = 1:runs
    global plane, newPlane, pictureArray
    resolutioni = size(plane)[1]
    resolutionj = size(plane)[2]

    for i1 = 1:resolutioni
        for j1 = 1:resolutionj
            newPlane[i1,j1] = cubicmethod(plane[i1,j1])
        end
    end

    plane = copy(newPlane)

    pictureArray[:,:,i] = newPlane
    save("Zoomed Cubic Fractals\\cubicfractal$i.png", map(clamp01nan,colorview(RGB, picture(pictureArray[:,:,i])/255)))
end

#save("Cubic Fractals\\cubicfractal.png", colorview(RGB, picture(pictureArray[:,:,runs]) / 255))
#{{x -> -1.39534}, {x -> 0.460355 - 1.13932 I}, {x -> 0.460355 + 1.13932 I}, {x -> 0.474627}}
