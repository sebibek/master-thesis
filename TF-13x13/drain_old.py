import numpy as np
import math as m
import os

def multiply(X, Y): ## yields the matrix product X*Y
    result = [[sum(a * b for a, b in zip(X_row, Y_col)) for Y_col in zip(*Y)] for X_row in X]
    return result

def rotate(matrix, phi):
    xx = m.cos(phi) # cosine of phi in radian
    xy = -m.sin(phi)
    yx = m.sin(phi)
    yy = m.cos(phi)
    rotation = np.array([[xx, xy], [yx, yy]])

    return np.matmul(rotation, matrix) ## rotation applied from the left (as last step)

def mirrorX(): # mirroring at x-axis
    mirror = np.array([[-1, 0], [0, 1]])
    return mirror

def mirrorY(matrix): # mirroring at y-axis
    mirror = np.array([[1, 0], [0, -1]])
    return mirror

## begin with scaling and mirroring (in arbitrary order) to correctly interpret your results!, rotations on top...

def scale(xx,yy):
    scale = np.array([[xx, 0], [0, yy]])
    return scale

def normalize(matrix): ## normalize uniformly to max(det(tensor_field))*Pi for real data to capture relative change (instead of directional)
    det = np.linalg.det(matrix)
    c = 1/np.sqrt(det)#*np.pi)

    scale = np.array([[c, 0], [0, c]])
    return np.matmul(scale, matrix)

# --> create parameter to control intensity of light src until saturation - inverse mode (black/white) (1-I): "Negativ"

# calculater of txt file len
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

identity = np.array([[1, 0],[0, 1]]) ## row-major ordering: row-by-row
length = 117
matrixArray = np.ndarray(shape=(length,length), dtype=np.ndarray) # initialize ndarray w 0

# generate normalized tensors from linear transformations
with open('temp.txt', 'wb') as f:
    for j in range(length): # rows
        for i in range(length): # cols
            y = length - 1 - round(j-(length-1)/2) # use inverted y coordinates (mirroring at x-axis in array reference frame)
            x = round(i-(length-1)/2) # --> center around midpoint
            matrixArray[j][i] = normalize(np.array([[1, 0], [0, 1]])) # initialization w. normed identity
            rad = m.atan2(y,x) # use vector atan2 to determine angle and quadrant of vector [x|y]
            if x==0 and y==0: # if center (origin)..
                scaled = np.matmul(scale(1, 1), identity)  # use isotropic scaling
            else:
                scaled = np.matmul(scale(4,1),identity) # use anisotropic scaling
            rotated = rotate(scaled, rad)
            normalized = normalize(rotated) # ..->transforms
            matrixArray[j][i] = normalized # update matrixArray entry

            np.savetxt(f, matrixArray[j][i], fmt='%s', delimiter=' ', newline='\r\n') # add matrix to temp.txt

# reorder the output tensor field into a regular grid
i = 1
str1 = ""
str2 = ""
with open('temp.txt', 'r+') as txtfile:
    with open('tensor_field.txt', 'w+') as tensor_field:
        for line in txtfile:
            line = line.rstrip() # remove newline/whitespace chars from line ends (the right)
            if i % 2 == 0: # if even line
                str2 += " " + line # concatenate row2 (lower line)
            else:
                str1 += " " + line # concatenate row1 (upper line)

            if i % (2*length) == 0: # if one row (in matrixArray) passed.., write results
                tensor_field.write(str1) # write row 1 (upper line)
                tensor_field.write('\n' + str2 + '\n') # write row 2 (lower line)
                str1 = str2 = "" # ..reset line strings

            i += 1 #increment line index

os.remove('temp.txt')