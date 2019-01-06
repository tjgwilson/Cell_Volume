from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import math
from math import *

image_width = 250
image_height = 250
radius = 100

a = image_width/2.
b = image_height/2.

im = Image.new("RGB", (image_width,image_height))
pix = im.load()

for x in range(image_width):
    for y in range(image_height):
        pix[x,y] = (0,0,0)


for x in range(image_width):
    for y in range(image_height):
        r = sqrt(((x-a))**2 + ((y-b)/0.4)**2)
        if(int(r) == radius or int(r) == -radius):
            for i in range(3):
                pix[x+i,y] = (200,200,0)



im.rotate(0).save("cell1.png", "PNG")
