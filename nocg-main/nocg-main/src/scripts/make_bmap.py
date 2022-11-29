"""

    make_bmap.py 

    Creates an image that can be used as a bump mapping texture.

    Mahesh Venkitachalam
    shader.in

"""

import numpy as np 
from PIL import Image
from math import sqrt

def main():

    NX, NY = 256, 256 

    nmap = np.zeros([NX, NY, 3], np.float32)

 
    r = 32.0
    rsq = r*r

    centers = [(64, 64), (192, 64), (64, 192), (192, 192)]

    for i in range(NX):
        for j in range(NY):
            inside = False
            for C in centers:
                x = (i-C[0])
                y = (j-C[1])
                if x*x + y*y < rsq :
                    nmap[i][j][0] = x / r
                    nmap[i][j][1] = y / r
                    nmap[i][j][2] = sqrt(rsq - (x*x + y*y))/ r
                    inside = True
            if not inside:
                nmap[i][j][0] = 0.0
                nmap[i][j][1] = 0.0
                nmap[i][j][2] = 1.0
    
    # [-1, 1] to [0, 255]
    nmap = 255.0*0.5*(nmap + 1.0)

    img = np.array(nmap, np.uint8)
    img = Image.fromarray(img)
    img.save("bmap.png")

# call main
if __name__ == '__main__':
    main()
