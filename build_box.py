import numpy as np
import cv2
import pandas as pd
from tqdm import tqdm
from numba import njit
import glob
import os
import xml.etree.cElementTree as ET

def write_xml(folder,filename,path,width,height,names,xmins,ymins,xmaxs,ymaxs,xml_name):
    
    annotation = ET.Element("annotation")
    ET.SubElement(annotation, "folder").text = folder 
    ET.SubElement(annotation, "filename").text = filename
    ET.SubElement(annotation, "path").text = path
    source = ET.SubElement(annotation, "source")
    ET.SubElement(source, "database").text = "Unknown"
    size = ET.SubElement(annotation, "size")
    ET.SubElement(size, "width").text = width 
    ET.SubElement(size, "height").text = height 
    ET.SubElement(size, "depth").text = "1"
    ET.SubElement(annotation, "segmented").text = "0"
    for i in range(len(names)):
        object_i = ET.SubElement(annotation, "object")
        ET.SubElement(object_i, "name").text = names[i]
        ET.SubElement(object_i, "pose").text = "Unspecified"
        truncated = ET.SubElement(object_i, "truncated").text = "0"
        ET.SubElement(object_i, "difficult").text = "0"
        bndbox = ET.SubElement(object_i, "bndbox")
        ET.SubElement(bndbox, "xmin").text = xmins[i]
        ET.SubElement(bndbox, "ymin").text = ymins[i]
        ET.SubElement(bndbox, "xmax").text = xmaxs[i]
        ET.SubElement(bndbox, "ymax").text = ymaxs[i]

    tree = ET.ElementTree(annotation)
    tree.write(xml_name)
    return None





def sub_process(img,i,out):
    k=0
    j=0
    flag=False
    neighbor=7
    index_max=-2*np.ones((256,256,2))
    index_min=1000*np.ones((256,256,2))
    while(j<256 and flag==False):
        if img[i,j]>200:
            k=k+1
            if k>2:
                index_max[i,j,0]=i
                index_max[i,j,1]=j
                index_min[i,j,0]=i
                index_min[i,j,1]=j
                flag=True
            
        Mi=int(index_max[:,:,0].max())
        Mj=int(index_max[:,:,1].max())
        mi=int(index_min[:,:,0].min())
        mj=int(index_min[:,:,1].min())
        out[mi-neighbor:Mi+neighbor,mj-neighbor:Mj+neighbor]=255
        j=j+1  
    return out
def area(box):
    xmin = box[0]
    xmax = box[1]
    ymin = box[2]
    ymax = box[3]
    return (ymax-ymin)*(xmax-xmin) 
def process(img):
    result=np.empty((256,256,3))
    result[:,:,0]=img
    result[:,:,1]=img
    result[:,:,2]=img
    out = np.zeros((256,256))
    for i in range(256):
        out=sub_process(img,i,out)
    ymin =0
    ymax =0
    xmin =0
    xmax =0 
    if (out.sum()>255):
        ymin = int((np.where(out>250)[0]).min())
        ymax = int((np.where(out>250)[0]).max())
        xmin = int((np.where(out>250)[1]).min())
        xmax = int((np.where(out>250)[1]).max())
    print('done!')
    return xmin,xmax,ymin,ymax
for k in tqdm(range(10,100,1)):
    b = '{0:03d}'.format(k+1)
    for i in range(1):
        names=[]
        xmins=[]
        ymins=[]
        xmaxs=[]
        ymaxs=[]
        a = '{0:03d}'.format(i)
        file_name = 'surface_'+b+'_'+a+'.png'
        img = cv2.imread(file_name)[:,:,0]
        result = cv2.imread('surface_'+b+'_comp_'+a+'.png')
        folder=os.path.split(os.getcwd())[-1]
        filename='surface_'+b+'_comp_'+a+'.png'
        path=os.path.join(os.getcwd(),filename)
        width = str(img.shape[0])
        height = str(img.shape[0])
        names.append('external_defect')
        for p1 in range(2):
            for p2 in range(2):
                img_patch = np.zeros(img.shape)
                img_patch[p1*128:(p1+1)*128,p2*128:(p2+1)*128] = img[p1*128:(p1+1)*128,p2*128:(p2+1)*128]  
                box = process(img_patch)
                xmin = box[0]
                xmax = box[1]
                ymin = box[2]
                ymax = box[3]
                xmins.append(str(xmin))
                ymins.append(str(ymin))
                xmaxs.append(str(xmax))
                ymaxs.append(str(ymax))
                xml_name = 'surface_'+b+'_comp_'+a+'.xml'
                result = cv2.rectangle(result,(xmin,ymin),(xmax,ymax),(0,0,255),1)
        cv2.imwrite('../'+str(k)+'_'+str(i)+'.png',result) 
        write_xml(folder,filename,path,width,height,names,xmins,ymins,xmaxs,ymaxs,xml_name)
