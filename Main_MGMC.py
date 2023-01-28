from __future__ import print_function, division
import os
import matplotlib
#matplotlib.use('Agg')
import torch
import numpy as np
from pathlib import Path
import pathlib
import argparse
from PIL import Image
from MGMC_solver import MGMC, calculate_psnr,imsave
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
from scipy import io as scio
import time
import math
import copy
import cv2
import torch.optim as optim
import glob
from matplotlib import pyplot as plt
from torchvision.utils import save_image
import torch.nn as nn

import warnings

warnings.filterwarnings("ignore")
import random

use_gpu = True
Image.LOAD_TRUNCATED_IMAGES = True
IMG_FORMATS = ['bmp', 'jpg', 'jpeg', 'png', 'tif', 'tiff', 'dng', 'webp', 'mpo']  # acceptable image suffixes
VID_FORMATS = ['mov', 'avi', 'mp4', 'mpg', 'mpeg', 'm4v', 'wmv', 'mkv']  # acceptable video suffixes
FILE = Path(__file__).resolve()
ROOT = FILE.parents[0]

def imsave(img, path):
    im = Image.fromarray(img.cpu().numpy().astype(np.uint8).squeeze())
    im.save(path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir',type=str, default=  './test')
    opt = parser.parse_known_args()[0]
    opt.device = 0 # 0,1, cpu
    cpu = opt.device == 'cpu'
    if cpu:
        os.environ['CUDA_VISIBLE_DEVICES'] = '-1'  # force torch.cuda.is_available() = False
    elif opt.device:
        os.environ['CUDA_VISIBLE_DEVICES'] = opt.device  # set environment variable
        assert torch.cuda.is_available(), f'CUDA unavailable, invalid device {opt.device} requested'  # check availability
    cuda = not cpu and torch.cuda.is_available()  #
    device = torch.device("cuda:0" if cuda else "cpu")
    clean = scio.loadmat('./lena.mat')['im']
    noise = clean + 20 * np.random.randn(clean.shape[0], clean.shape[1])
    clean = torch.from_numpy(clean).to(device, non_blocking=True).float()
    noise = torch.from_numpy(noise).to(device, non_blocking=True).float()
    MGMC_result= MGMC().to(device)
    MGMC_result.cuda()
    denoise=MGMC_result.MMC_code(clean,noise,device)
    imsave(denoise,'./denoise.png')
