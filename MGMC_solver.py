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


class MGMC(nn.Module):
    def __init__(self):
        super().__init__()
        self.alpha = 15
        self.max_level = 4


    def MMC_code(self,clean,noise,device):
        self.device =device
        self.clean=clean
        self.noise=noise
        self.alpha = 15
        self.max_level = 4
        [nrow,ncol] = self.noise.shape
        u = torch.zeros([nrow,ncol]).to(self.device) 
        self.num_iter =10
        self.num_outer=100
        J0 = self.energy(u,self.noise,self.alpha)
        self.Energy=J0

        t1=0
        error,count=[],0
        Energy_out, error_out=[],[]
        for outer_iter in range(0,self.num_outer):
            u_out_old=copy.deepcopy(u)
            #print(f'outer iteration{outer_iter}' )
            level_opt =torch.arange(self.max_level,0,-1)
            for idx in range(0,len(level_opt)):
                level = level_opt[idx]
                #print(f'level {level}')
                self.partition_grid(self.noise, self.alpha, level)
                start = time.time()
                for ii in range(1,self.num_iter):
                    u_old =copy.deepcopy(u)
                    u = self.MGMC_solver(u, self.noise,  level)
                    error.append( torch.linalg.norm(u - u_old) / torch.linalg.norm(u))
                    count = count + 1
                    if error[count - 2] < 1e-2:
                        break
                stop = time.time()
                t1=t1+stop-start
            Energy_out.append(self.energy(u,self.noise,self.alpha))
            error_out.append(torch.linalg.norm(u - u_out_old) / torch.linalg.norm(u))
            if error_out[outer_iter-1] <1e-4:
                break
        print(t1,'s')
        print(f' psnr {calculate_psnr(u.cpu().numpy(), self.clean.cpu().numpy())}')
        #fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        #ax.imshow(u.cpu().numpy(),cmap='gray')
        return u 


    def MGMC_solver(self,w,noise,level):
        w = self.MMC_DDM(noise, w, self.row_odd_idx, self.col_odd_idx, self.row_odd, self.col_odd, self.bnd_idx, self.bnd_idy, level, self.alpha,
                          self.noise_ch1, self.nphi_1, self.vx_1, self.vy_1, self.dv1, self.s1)
        w = self.MMC_DDM(noise, w, self.row_odd_idx, self.col_even_idx, self.row_odd, self.col_even, self.bnd_idx, self.bnd_idy, level, self.alpha,
                          self.noise_ch2, self.nphi_2, self.vx_2, self.vy_2, self.dv2, self.s2)
        w = self.MMC_DDM(noise, w, self.row_even_idx, self.col_odd_idx, self.row_even, self.col_odd, self.bnd_idx, self.bnd_idy, level, self.alpha,
                          self.noise_ch3, self.nphi_3, self.vx_3, self.vy_3, self.dv3, self.s3)
        w = self.MMC_DDM(noise, w, self.row_even_idx, self.col_even_idx, self.row_even, self.col_even, self.bnd_idx, self.bnd_idy, level, self.alpha,
                          self.noise_ch4, self.nphi_4, self.vx_4, self.vy_4, self.dv4, self.s4)
        return w
    def MMC_DDM(self,noise,w,row_odd_idx,col_odd_idx,row_odd,col_odd,bnd_idx,bnd_idy,level,alpha,noise_ch1,nphi_1,vx_1,vy_1,dv1,s1):
        w1= w[ row_odd_idx.type(torch.int64), col_odd_idx.T.type(torch.int64)]
        dw1 = self.MMC_fix(w1,bnd_idx,bnd_idy,row_odd,col_odd,level,alpha,noise_ch1,nphi_1,vx_1,vy_1,dv1,s1)
        if row_odd_idx[0]== row_odd_idx[1]:
            dw1[0,:]=dw1[1,:]
        if row_odd_idx[-1] == row_odd_idx[-2]:
            dw1[-1,:]=dw1[-2,:]

        if col_odd_idx[0] == col_odd_idx[1]:
            dw1[:, 0]=dw1[:, 1]

        if col_odd_idx[-1] == col_odd_idx[-2]:
            dw1[:, -1]=dw1[:, -2 ]
   
     
        w[row_odd_idx.type(torch.int64), col_odd_idx.T.type(torch.int64)] = w1 + dw1
        return w
    def MMC_fix(self,w1,bnd_idx,bnd_idy,row_odd,col_even,level,alpha,noise_ch1,nphi_1,vx_1,vy_1,dv1,s1):
 
        [pr, pc] =  self.ker(level).shape
        U1ch = torch.zeros([pr * pc,  len(row_odd)  *  len(col_even) ]).to(self.device)
        nr =  len(row_odd)
        nc =  len(col_even)
        for ii in range(0,nr):
            U1ch[:, (ii) * nc : (ii+1) * nc] = w1[(ii) * pr: (ii + 1) * pr, :].T.reshape(-1, pr * pc).T
        noisech = noise_ch1 
        zbar = noisech - U1ch
        nphi = nphi_1 
        zstar = zbar * nphi
        zstar = torch.sum(zstar, axis=0) / s1.to(self.device)
        ux = U1ch[ bnd_idx[:, 0].type(torch.int64) ,:]-U1ch[ bnd_idx[:, 1].type(torch.int64),:]
        uy = U1ch[ bnd_idy[:, 0].type(torch.int64) , :] - U1ch[bnd_idy[:, 1].type(torch.int64), :]
        vx =vx_1.to(self.device)
        vy=vy_1.to(self.device)
        c=torch.zeros([1,ux.shape[1]]).to(self.device)
        A1 =vx**2+vy**2
        a=torch.ones([ux.shape[0],1]).to(self.device)
        C1 =ux*vx+uy*vy
        for i in range(0,6):
            tmpc = torch.mm(a, c)
            ucx=ux+tmpc*vx
            ucy=uy+tmpc*vy
            B1=ucx**2+ucy**2+1e-4
            B1=torch.sqrt(B1)
            A= A1/B1
            firs=torch.sum(A,axis=0)
            C=C1/B1
            sec=torch.sum(C,axis=0)
            c = (s1*zstar-alpha*sec)/(alpha*firs+s1)
            c = c.unsqueeze(dim=0)

   
        dc1 =torch.kron(c.T.reshape(nr,nc) ,torch.ones([pr,pc]).to(self.device) )
 
        du1=dc1*dv1
        return du1


    def partition_grid(self,noise,alpha,level):
        [nr,nc]=noise.shape
        [M,N] = self.ker(level).shape       
        pr,pc=M,N
        j=level
        npr = math.floor( (nr-1)/ 2**(j-1) +1)
        npc = math.floor((nc - 1) / 2 ** (j - 1) + 1)
        self.row_odd = torch.arange(0,npr,2).to(self.device) 
        self.row_even =torch.arange(1,npr,2).to(self.device) 
        self.col_odd = torch.arange(0,npc,2).to(self.device) 
        self.col_even =torch.arange(1,npc,2).to(self.device) 
        row_odd_ctrx = 1 + (self.row_odd - 1) * 2 ** (j - 1)
        row_even_ctrx = 1 + (self.row_even - 1) * 2 ** (j - 1)
        col_odd_ctrx = 1 + (self.col_odd - 1) * 2 ** (j - 1)
        col_even_ctrx = 1 + (self.col_even - 1) * 2 ** (j - 1)
        self.row_odd_idx = torch.zeros([pr * len(self.row_odd), 1]).to(self.device) 
        for ii in range(0,len(self.row_odd)) :
            self.row_odd_idx[(ii) * pr : (ii+1) * pr, 0]=torch.arange( row_odd_ctrx[ii] - (pr - 1) / 2,1+row_odd_ctrx[ii] + (pr - 1) / 2 ,1)
        self.row_odd_idx[self.row_odd_idx < 0] = 0
        self.row_odd_idx[self.row_odd_idx > nr-1] = nr-1

        self.row_even_idx = torch.zeros([pr * len(self.row_even), 1]).to(self.device) 
        for ii in range(0,len(self.row_even)):
            self.row_even_idx[(ii) * pr: (ii+1) * pr, 0] = torch.arange(row_even_ctrx[ii] - (pr - 1) / 2, 1+row_even_ctrx[ii] + (pr - 1) / 2, 1)
        self.row_even_idx[self.row_even_idx < 0] = 0
        self.row_even_idx[self.row_even_idx> nr-1] = nr-1

        self.col_odd_idx =torch.zeros([pr * len(self.col_odd), 1]).to(self.device) 
        for ii in range(0,len(self.col_odd)) :
            self.col_odd_idx[(ii) * pr: (ii+1) * pr, 0]=torch.arange( col_odd_ctrx[ii] - (pr - 1) / 2,1+col_odd_ctrx[ii] + (pr - 1) / 2 ,1)
        self.col_odd_idx[self.col_odd_idx < 0] = 0
        self.col_odd_idx[self.col_odd_idx > nc-1] = nc-1

        self.col_even_idx =torch.zeros([pr * len(self.col_even), 1]).to(self.device) 
        for ii in range(0,len(self.col_even)):
            self.col_even_idx[(ii) * pr: (ii+1) * pr, 0] = torch.arange(col_even_ctrx[ii] - (pr - 1) / 2, 1+col_even_ctrx[ii] + (pr - 1) / 2, 1)
        self.col_even_idx[self.col_even_idx < 0] = 0
        self.col_even_idx[self.col_even_idx> nc-1] = nc-1
        #----------------------------------------------
        pr,pc=pr-2,pc-2
        M,N=pr,pc
        bnd_idx = torch.zeros([4*M,2]).to(self.device) 
        bnd_idy = torch.zeros([4*M,2]).to(self.device) 
        if M==1 and N==1:
            bnd_idy[:,0]=torch.from_numpy(np.array( [5,8,6,9]).T).to(self.device)
            bnd_idx[:, 0]=bnd_idy[:, 0]
            bnd_idx[:, 1]=bnd_idx[:, 0]-(M + 2)
            bnd_idy[:, 1]=bnd_idy[:, 0]-1
        else:
            bnd_idy[0: M + 1, 0]=torch.arange(( M + 2) + 2,(M + 2) + 2 + (M + 2) * (N)+1,(M + 2)  ).to(self.device) 
            bnd_idy[M + 1: M + 2 + M, 0]=bnd_idy[0: M + 1, 0]+N
            bnd_idy[2 * M + 2: 2 * M + 1 + M, 0]=torch.arange((M + 2) + 3,(M + 2) * 2 - 1+1,1 ).to(self.device) 
            bnd_idy[3 * M + 1: 4 * M , 0]=bnd_idy[2 * M + 2: 2 * M + 1 + M, 0]+(N) * (M + 2)
            bnd_idy[:, 1]=bnd_idy[:, 0]-1
            bnd_idx[:, 0]=bnd_idy[:, 0]
            bnd_idx[:, 1]=bnd_idx[:, 0]-(M + 2)
        self.bnd_idx = bnd_idx-1
        self.bnd_idy = bnd_idy-1
        #-------------------------------------------------
        ntau = pr * pc
        pr = M + 2
        pc = N + 2
        [self.noise_ch1, self.nphi_1, self.vx_1, self.vy_1, self.dv1, self.s1] = self.fix_iteration_pro(level, noise, pr, pc,
                                                                                             self.row_odd, self.col_odd,
                                                                                             self.row_odd_idx, self.col_odd_idx,
                                                                                             self.bnd_idx, self.bnd_idy)
        [self.noise_ch2, self.nphi_2, self.vx_2, self.vy_2, self.dv2, self.s2] = self.fix_iteration_pro(level, noise, pr, pc,
                                                                                             self.row_odd, self.col_even,
                                                                                             self.row_odd_idx, self.col_even_idx,
                                                                                             self.bnd_idx, self.bnd_idy)
        [self.noise_ch3, self.nphi_3, self.vx_3, self.vy_3, self.dv3, self.s3] = self.fix_iteration_pro(level, noise, pr, pc,
                                                                                             self.row_even, self.col_odd,
                                                                                             self.row_even_idx, self.col_odd_idx,
                                                                                             self.bnd_idx, self.bnd_idy)
        [self.noise_ch4, self.nphi_4, self.vx_4, self.vy_4, self.dv4, self.s4] = self.fix_iteration_pro(level, noise, pr, pc,
                                                                                             self.row_even, self.col_even,
                                                                                             self.row_even_idx, self.col_even_idx,
                                                                                             self.bnd_idx, self.bnd_idy)



    def fix_iteration_pro(self,level,noise,M,N,row_odd,col_odd,row_odd_idx,col_odd_idx,bnd_idx,bnd_idy):
        [nr, nc] =  noise.shape
        z1 = noise[ row_odd_idx.type(torch.int64), col_odd_idx.T.type(torch.int64) ]
        noise_ch1 = torch.zeros([M * N, len(row_odd) * len(col_odd)] ).to(self.device) 
        nr = len(row_odd)
        nc = len(col_odd)
        for ii in range(0,nr):
           noise_ch1[:, (ii) * nc : (ii+1) * nc]= z1[(ii) * M : (ii+1) * M,:].T.reshape(-1,M*N).T
        phi = self.ker(level)
        phitmp = phi.reshape(-1,1)
        nphi_1=torch.mm(phitmp, torch.ones([1, noise_ch1.shape[1]]).to(self.device)  )
        s1 = torch.sum(nphi_1*nphi_1,axis=0)
        vx_1=nphi_1[bnd_idx[:,0].type(torch.int64),:]-nphi_1[ bnd_idx[:,1].type(torch.int64),:]
        vy_1=nphi_1[bnd_idy[:,0].type(torch.int64),:]-nphi_1[ bnd_idy[:,1].type(torch.int64),:]
        dv1 = torch.zeros([M * len(row_odd), N * len(col_odd)]).to(self.device)
        for ii in range(0, len(row_odd)):
            tmp = nphi_1[:, (ii) * len(col_odd) : (ii+1) * len(col_odd) ]
            dv1[(ii) * M : (ii+1) * N,:]= tmp.T.reshape(-1,M).T
        return noise_ch1,nphi_1,vx_1,vy_1,dv1,s1
    def ker(self, level):
        if level==1:
            k_a=torch.zeros([2**(level-1)+2, 2**(level-1)+2 ] ).to(self.device) 
        else:
            k_a =torch.zeros([2 ** (level - 1) + 2,2 ** (level - 1) + 2] ).to(self.device) 
            [m, n] =  k_a.shape
            k_a = torch.zeros([m + 1, n + 1]).to(self.device) 
        [m, n]=  k_a.shape
        k_a[1:m-1,1:n-1]=1
        return k_a

    def plane_product(self,pr):
        plane_part_a = torch.zeros([pr - 1, 2]).to(self.device) 
        plane_part_b = torch.zeros([pr - 1, 2]).to(self.device) 
        plane_part_a[:, 0]=torch.arange(1,pr,1).to(self.device) 
        plane_part_a[:, 1]=1 + pr * 2 + (pr - 2) * 2 - plane_part_a[:, 0]
        plane_part_b[0:  ((pr - 1) / 2).type(torch.int64), 0]=torch.arange(pr,  1+pr + 2 * (pr - 2 - 1) / 2,2).to(self.device) 
        if pr ==3:
            plane_part_b[-1, 0] = torch.arange(pr + 2 * (pr - 2 - 1) / 2 + 1, pr, -2).to(self.device) 
        else:
            plane_part_b[0 +  ((pr - 1) / 2).type(torch.int64): -1, 0] =torch.arange(pr + 2 * (pr - 2 - 1) / 2 + 1, pr, -2).to(self.device) 
        plane_part_b[:, 1]=1 + pr * 2 + (pr - 2) * 2 - plane_part_b[:, 0]
        plane_part_a = plane_part_a.T
        plane_part_b = plane_part_b.T

        plane_a_group =torch.row_stack( (plane_part_a,plane_part_b[0] ))
        plane_c_group = torch.row_stack( (plane_part_a,plane_part_b[1] ))
        plane_c_group = torch.row_stack((plane_c_group[1,:],plane_c_group[0,:],plane_c_group[2,:] ) )

        plane_b_group = torch.row_stack( ( plane_part_b,plane_part_a[0,:]))
        if pr==3:
            plane_b_group = torch.column_stack((plane_b_group[:, 0], plane_b_group[:, -1]))
        else:
            plane_b_group = torch.column_stack((plane_b_group[:, 0], plane_b_group[:, -1], plane_b_group[:, 1: -2]))

        plane_d_group =torch.row_stack((plane_part_b,plane_part_a[1,:]))
        if pr==3:
            plane_d_group = torch.column_stack((plane_d_group[:, 0], plane_d_group[:, -1]))
        else:
            plane_d_group = torch.column_stack((plane_d_group[:, 0], plane_d_group[:, -1], plane_d_group[:, 1: -1]))

        plane_d_group = torch.row_stack((plane_d_group[1,:],plane_d_group[0,:],plane_d_group[2,:]))

        plane = torch.int32(torch.column_stack(( plane_a_group, plane_b_group, plane_c_group, plane_d_group)))
        I = torch.argsort(plane)
        plane = plane[:,I]-1
        return plane

    def energy(self,U,I,alpha):
        [m, n] = self.noise.shape
        DxU = torch.column_stack( ( U[:, 0: n - 1] - U[:, 1: n], torch.zeros([m, 1]).to(self.device)  )  )
        DyU =torch.row_stack( ( U[0:m - 1,:] - U[1: m,:], torch.zeros([1, n]).to(self.device) ))
        J_U = sum(sum((U - I)** 2)) + alpha * sum(sum(torch.sqrt(DxU ** 2 + DyU ** 2)))
        return J_U




    def save_result(self, mod, epoch, title, outputs, x_qry, y_qry):
        # -------------------------------end finetunning and save some results----------------------------------------

        temp_out = outputs
        temp_input = x_qry
        temp_label = y_qry

        temp_out = temp_out.detach().cpu().numpy()
        temp_label = temp_label.detach().cpu().numpy()
        temp_input = temp_input.detach().cpu().numpy()
        # psnr1 = self.calculate_psnr(temp_out[0], temp_label[0])
        num_img = 5 if len(temp_out) > 6 else len(temp_out)
        fig, ax = plt.subplots(num_img, 3, figsize=(6, 6))
        if len(temp_out) == 1:
            [ax[0].imshow(temp_out[i].transpose((1, 2, 0))) for i, _ in enumerate(temp_out[0: num_img])]
            [ax[1].imshow(temp_input[i].transpose((1, 2, 0))) for i, _ in enumerate(temp_input[0: num_img])]
            [ax[2].imshow(temp_label[i].transpose((1, 2, 0))) for i, _ in enumerate(temp_label[0: num_img])]
        else:
            [ax[i][0].imshow(temp_out[i].transpose((1, 2, 0))) for i, _ in enumerate(temp_out[0: num_img])]
            [ax[i][1].imshow(temp_input[i].transpose((1, 2, 0))) for i, _ in enumerate(temp_input[0: num_img])]
            [ax[i][2].imshow(temp_label[i].transpose((1, 2, 0))) for i, _ in enumerate(temp_label[0: num_img])]

        f = self.train_result / f'{mod}_batch{epoch}_labels.png'
        plt.title(title)

        plt.savefig(f)
        plt.close()
    def save_val_result(self, mod, epoch, title, outputs, x_qry, y_qry):
        # -------------------------------end finetunning and save some results----------------------------------------

        temp_out = outputs
        temp_input = x_qry
        temp_label = y_qry

        fig, ax = plt.subplots(1, 3, figsize=(6, 6))
        ax[0].imshow(temp_out)
        ax[1].imshow(temp_input)
        ax[2].imshow(temp_label)
        f = self.train_result / f'{mod}_batch{epoch}_labels.png'
        plt.title(title)
        plt.savefig(f)
        plt.close()



def np2torch(x,device):
    x = torch.from_numpy(x).to(device, non_blocking=True).float()
    return x





def calculate_psnr(img1, img2, border=0):
    # img1 and img2 have range [0, 255]
    # img1 = img1.squeeze()
    # img2 = img2.squeeze()

    if not img1.shape == img2.shape:
        raise ValueError('Input images must have the same dimensions.')
    h, w = img1.shape[:2]
    img1 = img1[border:h - border, border:w - border]
    img2 = img2[border:h - border, border:w - border]

    img1 = (img1 ).astype(np.float64)
    img2 = (img2 ).astype(np.float64)
    mse = np.mean((img1 - img2) ** 2)
    if mse == 0:
        return float('inf')
    return 20 * math.log10(255 / math.sqrt(mse))


def imsave(img, path):
    im = Image.fromarray(img.cpu().numpy().astype(np.uint8).squeeze())
    im.save(path)




