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
        red = 255
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
    if x <= -3
        blue = 255
    elseif x >=-3 && x <= -2
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
    p = zeros(3, resolutioni, resolutionj)

    for i1 = 1:resolutioni
        line = zeros(3, resolutioni)

        for j1 = 1:resolutionj
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

            line[:,j1] = c
        end

        p[:,:,i1] = line
    end

    return p
end

function f(x)
    return x^5 - x + 1
end

function df(x)
    return 5x^4 - 1
end

function newton(x)
    return x - f(x)/df(x)
end

tmpa = transpose(2*collect(-r:r-1)/r);
tmpb = ones(1,imgSize);
plane = (transpose(tmpa)*tmpb)*1im;
plane = plane + (transpose(tmpb)*tmpa);
pictureArray = zeros(imgSize, imgSize, runs)*0im;
newPlane = zeros(imgSize, imgSize)*0im;

for i = 1:runs
    global plane, newPlane, pictureArray
    resolutioni = size(plane)[1]
    resolutionj = size(plane)[2]

    for i1 = 1:resolutioni
        for j1 = 1:resolutionj
            newPlane[i1,j1] = newton(plane[i1,j1])
        end
    end

    plane = newPlane

    pictureArray[:,:,i] = newPlane
    save("Newton Fractals\\newtonfractal$i.png", colorview(RGB, picture(pictureArray[:,:,i]) / 255))
end

#save("newtonfractal.png", colorview(RGB, picture(pictureArray[:,:,runs]) / 255))
