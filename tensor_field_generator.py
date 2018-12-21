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
    c = 1/np.sqrt(det*np.pi)

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

# generate normalized tensors from linear transformations (cartesian form)
# with open('temp.txt', 'wb') as f:
#     for j in range(length): # rows
#         for i in range(length): # cols
#             matrixArray[j][i] = normalize(np.array([[1, 0],[0, 1]])) # initialization w. normed identity
#             scaled = np.matmul(scale(4-i*(3/(length-1)),1),identity) # chronological transformation order: right->left
#             rotated = rotate(scaled, rad)
#             normalized = normalize(rotated)
#             matrixArray[j][i] = normalized   # create
#
#             np.savetxt(f, matrixArray[j][i], fmt='%s', delimiter=' ', newline='\r\n')

# generate normalized tensors from linear transformations (polar form)
radrange = 12*m.pi
radres = m.pi/180 # 2deg/step
radarr = [float(i)/57.296 for i in range(round(radrange/radres))]
curDeltaX = [[1.0]*length]*length
curDeltaY = [[1.0]*length]*length

for j in range(length): # rows
    for i in range(length): # cols
        matrixArray[j][i] = normalize(np.array([[1, 0], [0, 1]]))  # initialization w. normed identity

for rad in radarr:
    r = 0.1*rad
    x = round((length-1)/2+r*m.cos(rad))
    y = round((length-1)/2+r*m.sin(rad))
    rot = rad + 0.35
    deltaX = abs(x - ((length-1)/2+ r*m.cos(rad)))
    deltaY = abs(y - ((length-1)/2+ r*m.sin(rad)))
    if y < length and x < length:

        if deltaX<curDeltaX[y][x] and deltaY<curDeltaY[y][x]:
            scaled = np.matmul(scale(4, 1), identity)  # chronological transformation order: right->left
            rotated = rotate(scaled, rot)
            normalized = normalize(rotated)
            matrixArray[y][x] = normalized  # create

            curDeltaX[y][x] = deltaX
            curDeltaY[y][x] = deltaY
            xfit = x
            yfit = y


with open('temp.txt', 'wb') as f:
    for j in range(length): # rows
        for i in range(length): # cols
            np.savetxt(f,  matrixArray[j][i], fmt='%s', delimiter=' ', newline='\r\n')

print("filelen: " + str(file_len('temp.txt')))
print("matrixarray: " + str(len(matrixArray)))

# reorder the output tensor field into a regular grid
i = 1
str1 = ""
str2 = ""
with open('temp.txt', 'r+') as txtfile:
    for line in txtfile:
        line = line.rstrip()
        if i == 1:
            str1 += line
        elif i==1:
            str2 += line
        elif i % 2 == 0:
            str2 += " " + line
        else:
            str1 += " " + line

        if i%(2*length) == 0:
            txtfile.write('\n')
            txtfile.write(str1)
            txtfile.write('\n')
            txtfile.write(str2)
            str1 = ""
            str2 = ""

        i += 1


with open('temp.txt', 'r') as fin:
    data = fin.read().splitlines(True)
with open('tensor_field.txt', 'w+') as fout:
    fout.writelines(data[length*length*2+1:]) # in 2D...

#os.remove('temp.txt')
