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

function picture(plane,stripes)
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
            d = sqrt(i^2 + j^2)+(stripes[i1,j1]%5)/10

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
    tmpa = transpose(fw*collect(-xhalf:xhalf-1)/xhalf);
    tmpb = ones(1,xsize);
    plane = ((transpose(tmpa)*tmpb).+cy)*1im;
    plane = plane + ((transpose(tmpb)*tmpa).+cx);
    plane = plane[Int64(xhalf-yhalf):Int64(xhalf+yhalf-1),:]
    #tmpa = transpose(fw*collect(-r:r-1)/r);
    #tmpb = ones(1,imgSize);
    #plane = ((transpose(tmpa)*tmpb).+cy)*1im;
    #plane = plane + ((transpose(tmpb)*tmpa).+cx);
    #pictureArray = zeros(imgSize, imgSize, numruns)*0im;
    #newPlane = zeros(imgSize, imgSize)*0im;

    colorwheel = copy(plane)

    #save("Cubic Fractals\\$(fn)cubic0.png", map(clamp01nan,colorview(RGB, picture(plane)/255)))

    #global plane, newPlane, pictureArray
    resolutioni = size(plane)[1]
    resolutionj = size(plane)[2]
    global stripes = zeros(Int64,resolutioni, resolutionj)

    for i1 = 1:resolutioni
        for j1 = 1:resolutionj
            for i = 1:numruns
                t = cubicmethod(plane[i1,j1])
                if abs(t-plane[i1,j1]) < 1e-8*abs(t)
                    stripes[i1,j1] = i
                    break
                end
                plane[i1,j1] = t
            end
        end
    end

    #pictureArray[:,:,i] = plane

    save("Cubic Fractals\\cubic$(fn).png", map(clamp01nan,colorview(RGB, picture(plane,stripes)/255)))
end

#save("Cubic Fractals\\cubicfractal.png", colorview(RGB, picture(pictureArray[:,:,runs]) / 255))
#{{x -> -1.39534}, {x -> 0.460355 - 1.13932 I}, {x -> 0.460355 + 1.13932 I}, {x -> 0.474627}}

#createFractal(-1.739674185463659, 1.3742230576441103, 0.020025062656641612, "z")
#createFractal(0.0, 0.0, 2.0, "a")
function getnewrun()
    m1 = 0
    m2 = 0
    for i = 1:size(stripes)[1]
        for j = 1:size(stripes)[2]
            t = stripes[i,j]
            if m1 < t
                m2 = m1
                m1 = t
            elseif m2 < t
                m2 = t
            end
        end
    end
    return m2
end
#Center (x,y) Distance
#The y val gets flipped so we need to do 400-y
#84 274(126) to 86 276(124) so Center = 85 275(125)
#findx(x) = -2 + 4*(x/400)
#findy(y) = 2 - 4*(y/400)
#wa = [0.0 0.0 2.0]
#wz = [-1.15 0.75 0.01]
#201 67 to 203 69 so Center = 202 68
#findx(x) = -1.16 + 0.02*(x/400)
#findy(y) = 0.74 + 0.02*(y/400)
#wa2 = wz
#wz2 = [-1.1499 0.7434 0.00005]
#272 340 to 274 342 so Center = 273 341
#findx(x) = -1.14995 + 0.0001*(x/400)
#findy(y) = 0.74335 + 0.0001*(y/400)
#wa3 = wz2
#wz3 = [-1.14988175 0.74343525 0.00000025]
initialval = [0.0 0.0 2.0]
vals = zeros(10, 3)

open("keyframes.txt") do f
  line = 0

  while ! eof(f)
     s = readline(f)
     sprime = split(s, " ")
     vals[line+1,1] = parse(Float64,sprime[1])
     vals[line+1,2] = parse(Float64,sprime[2])
     vals[line+1,3] = parse(Float64,sprime[3])
     line += 1
  end
end

nmax = 100
alphas = zeros(10)
alphas[1] = (1/nmax)*log(vals[1,3]/initialval[3])
for i = 2:10
    alphas[i] = (1/nmax)*log(vals[i,3]/vals[i-1,3])
end
#alpha1 = (1/nmax)*log(wz[3]/wa[3])
#Threads.@threads before the for loop, but on the same line
global numruns = 60

for j = 0:nmax
    global numruns
    theta = (initialval[3]/(vals[1,3]-initialval[3]))*(exp(alphas[1]*j)-1)
    wn = (1-theta)*initialval[:]+theta*vals[1,:]
    createFractal(wn[1], wn[2], wn[3], j, numruns)
    numruns = 3*getnewrun()÷2 + 1

    println(numruns, "     ", j, "     ", wn)
end

n = nmax+1
counter = 2
while n%100 != 0
    theta = (vals[counter-1,3]/(vals[counter,3]-vals[counter-1,3]))*(exp(alphas[counter]*(n-(counter-1)*nmax))-1)
    wn = (1-theta)*vals[counter-1,:]+theta*vals[counter,:]
    createFractal(wn[1], wn[2], wn[3], n, numruns)
    if n%100 == 99
        global n += 1
        theta = (vals[counter-1,3]/(vals[counter,3]-vals[counter-1,3]))*(exp(alphas[counter]*(n-(counter-1)*nmax))-1)
        wn = (1-theta)*vals[counter-1,:]+theta*vals[counter,:]
        createFractal(wn[1], wn[2], wn[3], n, numruns)
        println(numruns, "     ", n, "     ", wn)
        global counter += 1
    end

    global numruns = 3*getnewrun()÷2 + 1

    println(numruns, "     ", n, "     ", wn)
    if n == 1000
        break
    else
        global n += 1
    end
end

    #println(oldsize/wn[3], " ", temptest, " ", wn)
    #global oldsize = wn[3]

#m(62) = 1.3742230576441103
#l(26) = -1.739674185463659
#width = 0.020025062656641612
