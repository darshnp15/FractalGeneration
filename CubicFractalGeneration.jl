using Images, Colors, ImageIO, ImageMagick, Printf

#imgSize = 400
#imgSize = 400
#r = imgSize/2
xsize = 1920
xhalf = xsize÷2
ysize = 1080
yhalf = ysize÷2

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
    p = zeros(3, resolutioni, resolutionj)

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

    colorwheel = copy(plane)
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

    save("Cubic Fractals/cubic"*@sprintf("%06d",fn)*".png",
		map(clamp01nan,colorview(RGB, picture(plane,stripes)/255)))
end

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

nmax = zeros(Int,10)
alphas = zeros(10)
#nmax[1] = 500
nmax[1] = 50
alphas[1] = (1/nmax[1])*log(vals[1,3]/initialval[3])
for i = 2:10
	nmax[i] = Int(trunc(log(vals[i,3]/vals[i-1,3])/alphas[i-1]+0.5))
    alphas[i] = (1/nmax[i])*log(vals[i,3]/vals[i-1,3])
end
global numruns = 60

# nstart -- start frame
# nend -- end frame
# a  -- alpha coefficient in exp(alpha*t)
# v  -- 2x3 vector consisiting of
#             initial_x   initial_y   initial_size
#               final_x     final_y     final_size
function dowork(nstart,nend,a,v)
	for j = nstart:nend
	    global numruns
		theta = v[1,3]/(v[2,3]-v[1,3])*(exp(a*(j-nstart))-1)
		wn = (1-theta)*v[1,:]+theta*v[2,:]
	    createFractal(wn[1],wn[2],wn[3],j,numruns)
		numruns=3*getnewrun()÷2+1
	    println(numruns, "     ", j, "     ", wn)
	end
end

println("Computing segment 0 with alpha=",alphas[1],
		" from 0 to ",nmax[1],"...")
dowork(0,nmax[1],alphas[1],[initialval; vals[1,:]'])
for i = 2:10
	println("Computing segment $i with alpha=",alphas[i],
		" from ",sum(nmax[1:i-1])-1," to ",sum(nmax[1:i]),"...")
	dowork(sum(nmax[1:i-1])-1,sum(nmax[1:i]),alphas[i],[vals[i-1,:]'; vals[i,:]'])
end
