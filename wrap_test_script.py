# https://stackoverflow.com/questions/43048725/python-creating-video-from-images-using-opencv
import cv2
import numpy as np
import glob
from icecream import ic

img=[]
inPath = '/Users/dhueholt/Documents/GLENS_fig/20220106_snapshots/8_animate/1_6p/'
globsList = sorted(glob.glob(inPath+'*.png'))
for g in globsList:
    img.append(cv2.imread(g))

height,width,layers=img[1].shape

fourcc = cv2.VideoWriter_fourcc(*'mp4v')
video=cv2.VideoWriter(inPath+'animate_2yr_window.mp4',fourcc,5,(width,height))

for gc,gv in enumerate(globsList):
    video.write(img[gc])

cv2.destroyAllWindows()
video.release()
