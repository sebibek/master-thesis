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
length = 5

deg = -45
rad = deg*(np.pi/180.0)

matrixArray = np.ndarray(shape=(length,length), dtype=np.ndarray) # initialize ndarray w 0

# generate normalized tensors from linear transformations (polar form)
radrange = 2*m.pi
radres = m.pi/180 # 1deg/step
steps = radrange/radres
radarr = [float(i)*radrange/(steps-1) for i in range(round(steps))]
curDeltaX = [[1.0]*length]*length
curDeltaY = [[1.0]*length]*length
curDelta = np.ones((length,length))

for j in range(length): # rows
    for i in range(length): # cols
        matrixArray[j][i] = normalize(np.array([[1, 0], [0, 1]]))  # initialization w. normed identity

# create range for radius (radii)
radirange=range(m.ceil(length/2))
# walk through the radrange in radian steps
for i in radirange:
    for rad in radarr:
        if i == 0:
            r = i # set current radius
            x = round(r * m.cos(rad)) # get current cartesian x coordinate
            y = round(r * m.sin(rad)) # get current cartesian y coordinate
            rot = m.pi / 2  # use angle w. offset
            xIndex = round(x + (length - 1) / 2) # use shifted x index to obtain origin in center
            yIndex = length - 1 - round(y + (length - 1) / 2)  # use inverted y coordinates (mirroring at x-axis in array reference frame)
            scaled = np.matmul(scale(1, 1), identity)  # chronological transformation order: right->left
            rotated = rotate(scaled, rot)
            normalized = normalize(rotated) # ..->transforms
            matrixArray[yIndex][xIndex] = normalized  # update matrixArray entry
            continue

        r = i # set current radius
        x = round(r * m.cos(rad))  # get current cartesian x coordinate
        y = round(r * m.sin(rad))  # get current cartesian y coordinate
        rot = rad + m.pi/2 # use angle w. offset
        deltaX = abs(x - (r * m.cos(rad))) # calc errorX
        deltaY = abs(y - (r * m.sin(rad))) # calc errorX
        delta = m.sqrt(deltaX**2 + deltaY**2) # calc euclidean distance as absolute error --> minimize

        xIndex = round(x+(length-1)/2) # use shifted x index to obtain origin in center
        yIndex = length-1-round(y+(length-1)/2) # use inverted y coordinates (mirroring at x-axis in array reference frame)
        if xIndex < length and yIndex < length: # if indices inside range..
            if delta<curDelta[yIndex][xIndex]: # if absolute error < current error.. update
                scaled = np.matmul(scale(4, 1), identity)  # chronological transformation order: right->left
                rotated = rotate(scaled, rot)
                normalized = normalize(rotated) # ..->transforms
                matrixArray[yIndex][xIndex] = normalized  # update matrixArray entry
                curDelta[yIndex][xIndex] = delta # update current error


#iterate through matrixArray to add matrices in order
with open('temp.txt', 'wb') as f:
    for j in range(length): # rows
        for i in range(length): # cols
            np.savetxt(f,  matrixArray[j][i], fmt='%s', delimiter=' ', newline='\r\n')

print("filelen: " + str(file_len('temp.txt'))) # print filelen of temp: should equal length*length*matrixHeight
print("matrixarray: " + str(len(matrixArray))) # print len of matrixArray: should equal length

# reorder the output tensor field into a regular grid
i = 1
str1 = ""
str2 = ""
with open('temp.txt', 'r+') as txtfile:
    for line in txtfile:
        line = line.rstrip() # remove newline/whitespace chars from line ends (the right)
        if i % 2 == 0: # if even line
            str2 += " " + line # concatenate row2 (lower line)
        else:
            str1 += " " + line # concatenate row1 (upper line)

        if i%(2*length) == 0: # if one row (in matrixArray) passed.., write results
            txtfile.write('\n' + str1) # write row 1 (upper line)
            txtfile.write('\n' + str2) # write row 2 (lower line)
            str1 = str2 = "" # ..reset line strings

        i += 1 #increment line index

with open('temp.txt', 'r') as fin:
    data = fin.read().splitlines(True)
with open('tensor_field.txt', 'w+') as fout:
    fout.writelines(data[length*length*2+1:]) # in 2D...

os.remove('temp.txt')
