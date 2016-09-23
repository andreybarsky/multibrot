from PIL import Image
import sys, colorsys, cmath, numpy, images2gif

def make_mbrot(d, scale=200, itermax=30):
    filename = "mbrot\\mbrot%s.png" % str(d)
    escape = 2**2

    scale = scale
    itermax = itermax
    d = float(d)

    imgx = scale*4   # width of image, in pixels
    imgy = scale*3   # height of image, in pixels

    xaxis = imgy/2   # position of horizontal axis along vertical axis, in pixels
    yaxis = imgx/1.4 # position of vertical axis along horizontal axis, in pixels

    mbrot = Image.new('RGB',(int(imgx),int(imgy)))     # draws blank image

    def scaledx(x):     # converts pixel x-coordinate to complex plane x-coordinate
        newx = float(x) - float(yaxis)
        newx = newx / float(scale)
        return newx
        
    def scaledy(y):     # converts pixel y-coordinate to complex plane y-coordinate
        newy = float(y) - float(xaxis)
        newy = newy / float(scale)
        return 0.0 - newy
        
    # def sqabsol(c): # squared absolute value of complex number
        # return abs(c**2)

    def checkset(x,y,d):  # check if pixel is part of mandelbrot set    
        cx = scaledx(x) # scaled x coordinate on complex plane
        cy = scaledy(y) # scaled y coordinate on complex plane
        c = complex(cx, cy)
        iter = 0
        # zx = 0.0 # real part
        # zy = 0.0 # imaginary part
        z = complex(0.0, 0.0)
            
        while (absc(z**2) < escape) and (iter < itermax):   # while z has not yet escaped...
            z = z**d + c
            iter += 1
        
        return iter    # return number of iterations it took to escape, up to max
        
    def rgb(h,s,v):    # defines rgb color values (0-255) from hsv color values (0-255)
        starthue = -0.2 # arbitrary
        h,s,v = (float(h)/255)+starthue,float(s)/255,float(v)/255
        r,g,b = colorsys.hsv_to_rgb(h,s,v)
        r = int(r*255); g = int(g*255); b = int(b*255)
        return (r,g,b)
            
    pic = mbrot.load()
    colourscale = 255/itermax  # scales colour hue to max number of iterations

    print "Computing..."  
    onepercent = float(imgy) / 100 # progress bar stuff
    newcent = onepercent
    centval = 0

    for y in range(0,int(imgy)):          # draw the picture...
        for x in range(0,int(imgx)):
            check = checkset(x,y,d)  
            if check == itermax: # if pixel has never escaped;
                pic[x,y] = (255,255,255) # colour in white
            else:
                pixval = 255 - (check*colourscale)  # else, colour according to the hue gradient set, based on the number of iterations to escape
                #pixval = rgb(pixval, 255, 200) # but we're just leaving it in black and white for now
                pic[x,y] = (255-pixval,255-pixval,255-pixval)

        if y > newcent:         # progress bar
            if centval % 5 == 0:
                print str(centval) + "%",
            centval += 1
            newcent += onepercent        

    mbrot.save(filename)
    print "\nSaved as %s" % filename

drange = numpy.arange(0.0,8.1, 0.1)
# main loop:
for x in drange:
    make_mbrot(float(x))

gifrange = numpy.arange(1.0, 6.1, 0.1)
gifrange2 = numpy.arange(5.9, 1.0, -0.1)
files = [("mbrot\\mbrot%s.png" % x) for x in numpy.concatenate([[1.0, 1.0], gifrange, [6.0, 6.0], gifrange2])]
images = [Image.open(fn) for fn in files]

images2gif.writeGif("mbrot\\mbrot.gif", images, duration=0.05)
