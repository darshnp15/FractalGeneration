using Images, Colors, ImageIO, ImageMagick, ImageView

imgSize = 400

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

            line[:,resolutioni+1-i1] = c
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
	r = imgSize/2
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

    save("data/cubic$(fn).png", map(clamp01nan,colorview(RGB, picture(plane,stripes)/255)))
end

#save("Cubic Fractals\\cubicfractal.png", colorview(RGB, picture(pictureArray[:,:,runs]) / 255))
#{{x -> -1.39534}, {x -> 0.460355 - 1.13932 I}, {x -> 0.460355 + 1.13932 I}, {x -> 0.474627}}

#createFractal(-1.739674185463659, 1.3742230576441103, 0.020025062656641612, "z")

function main()
	ImageView.closeall()
    r = imgSize/2
	cx=0.0; cy=0.0; fw=2.0
	for i=1:10
		createFractal(cx,cy,fw,"tmp",200)
		img=load("data/cubictmp.png")
		guidict=imshow(img)
		s=readline()
		if s=="quit"
			break
		end
		println(guidict["roi"]["zoomregion"].value.currentview)
		xrange=guidict["roi"]["zoomregion"].value.currentview.x
		yrange=guidict["roi"]["zoomregion"].value.currentview.y
		ci=(xrange.left+xrange.right)/2
		cj=(yrange.left+yrange.right)/2
	    cw=(xrange.right-xrange.left)/2
	    delta=(yrange.right-yrange.left)/2
		if cw<delta
			cw=delta
		end
		println("ci=$(ci) cj=$(cj) cw=$(cw)")
		cx=cx+fw*(ci-r-1)/r
		cy=cy+fw*(400-cj-r)/r
#		cy=cy+fw*(cj-r-1)/r
		fw=fw*(cw+1)/r
		println("cx=$(cx) cy=$(cy) fw=$(fw)")
        open("keyframes.txt","a") do io
            println(io,cx," ",cy," ", fw)
        end
	end
end
