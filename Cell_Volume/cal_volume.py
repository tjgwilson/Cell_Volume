"""
Calculates volume of cell outline drawn on blank background.
T.Wilson 23/12/18

future updates: Speed up major axis finder - slow optimize routine
                Move functions to class - as originally indended
                Background function - allow for non blank Background
                Paralize loops
                Increase Accuracy

"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import math
from math import *

###############################################################
###############################################################
pixel_size = 1 #i.e. how many mm/m/km per pixel
integ_steps = 100 #number of steps to take along axis of rotation numbers
input_file = 'Cell4.png' #name of input file
output_file = 'output_'+input_file
###############################################################
###############################################################


#dynamic memory expansion
def expand_array(array):
    size = np.shape(array)
    array1 = np.zeros((2*size[0],size[1]))
    for i in range(size[0]):
        array1[i,0] = array[i,0]
        array1[i,1] = array[i,1]
    return array1

def cut_array(array,array_size):
    new_array = np.zeros((array_size,2))
    for i in range(array_size):
        new_array[i,0] = array[i,0]
        new_array[i,1] = array[i,1]
    return new_array

def find_pixels(pix,background):
    array_size = 100
    points = np.zeros((array_size,2))
    count = 0
    for x in range(image_width):
        for y in range(image_height):
            if(count >= array_size):
                points = expand_array(points)
                array_size *= 2

            if(pix[x,y][0] > background[0] or pix[x,y][1] >  background[1] or pix[x,y][2] > background[2]):
                points[count,0] = x
                points[count,1] = y
                count += 1
    points = cut_array(points,count)
    return points

#this takes the average background
#adjust this function if it does not correctly
#identitfy your cell boundary ie if your
#background is not black enough
def find_background_average(pix):
    sum = np.zeros((3))
    for x in range(image_width):
        for y in range(image_height):
            for i in range(3):
                sum[i] += pix[x,y][i]
    sum /= image_width*image_height
    print("Average background RGB value is: {}".format(sum))
    return sum

def find_major_axis(pixels,num_pixels):
    pos = np.zeros((2,2))
    temp = np.zeros((2))
    old_dist = 0.
    for j in range(num_pixels):
        temp[0] = pixels[j,0]
        temp[1] = pixels[j,1]
        for i in range(num_pixels):
            dist = sqrt((pixels[i,0] - temp[0])**2 + (pixels[i,1] - temp[1])**2)
            if(dist > old_dist):
                pos[0,0] = temp[0]
                pos[0,1] = temp[1]
                pos[1,0] = pixels[i,0]
                pos[1,1] = pixels[i,1]
                old_dist = dist
    return pos

def draw_axis(pix,axis):
    grad = (axis[1,1] - axis[0,1]) / (axis[1,0] - axis[0,0])
    c = axis[0,1] - (grad*axis[0,0])

    for x in range(image_width):
        for y in range(image_height):
            if( abs(y - ((grad*x) + c)) <= 0.2):
                pix[x,y] = (0,250,0)

def find_axis_angle(axis):
    theta = atan((axis[1,1] - axis[0,1]) / (axis[1,0] - axis[0,0]))
    print("angle of axis to horizontal is {0:.2f} deg".format(theta*180/pi))
    return theta

def find_widest_points_at_xy(x,y,pixels,num_pixels,pix):

    pix[x,y] = (250,0,0)
    grad = (axis[1,1] - axis[0,1]) / (axis[1,0] - axis[0,0])
    grad = - 1. / grad
    c = y - (grad*x)

    array_size = 2
    count = 0
    temp = np.zeros((array_size,2))
    dist = 0.
    prec = 1

    while(1):
        for i in range(num_pixels):
            if(abs(pixels[i,1] - (grad*pixels[i,0]) - c) <= prec):
                if(count >=array_size):
                    temp = expand_array(temp)
                    array_size *= 2

                temp[count,0] = pixels[i,0]
                temp[count,1] = pixels[i,1]
                count += 1

        greatest_sep = find_major_axis(cut_array(temp,count),count)
        distance = sqrt((greatest_sep[1,0] - greatest_sep[0,0])**2 + (greatest_sep[1,1] - greatest_sep[0,1])**2)

        if(count <= 10 or distance <= 2):
            prec += 0.5
        elif(count >= 1000):
            return np.inf
        else:
            break
    pix[greatest_sep[0,0],greatest_sep[0,1]] = (250,0,0)
    pix[greatest_sep[1,0],greatest_sep[1,1]] = (250,0,0)
    return distance

def integrate_area(pixels,num_pixels,axis,pix):
    volume = 0.0
    axis_length = sqrt((axis[1,0] - axis[0,0])**2 + (axis[1,1] - axis[0,1])**2)

    theta = find_axis_angle(axis)
    length = 0.
    dl = axis_length / integ_steps
    pos = axis[0,:]
    temp_diametre = 0.
    while(length <= axis_length):
        diametre = find_widest_points_at_xy(pos[0],pos[1],pixels,num_pixels,pix)
        if(diametre == "inf"):
            diametre = temp_diametre
            print("no, pixels found, using neares neighbour")
        else:
            volume += (0.25 * pi * (diametre*pixel_size)**2 *dl*pixel_size)
        pos[0] += dl * cos(theta)
        pos[1] += dl * sin(theta)
        length += dl
    return volume


###############################################################
###############################################################

im = Image.open(input_file)
pix = im.load()

image_width, image_height = im.size
print("Image dimensions {}x{}".format(image_width,image_height))

pixels = find_pixels(pix,find_background_average(pix))
num_pixels = np.shape(pixels)[0]
print("Number of coloured pixels found = {}".format(num_pixels))

#re-colouring to blue so you can check it correctly finds cell boundary
for i in range(num_pixels):
    pix[pixels[i,0],pixels[i,1]] = (0,0,250)


axis = find_major_axis(pixels,num_pixels)
draw_axis(pix,axis)
volume = integrate_area(pixels,num_pixels,axis,pix)

print(volume)

im.rotate(0).save(output_file,"PNG")
