#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from scipy import interpolate
import pdb

from flamesProps import*

class MarksteinNumber:

    def __init__(self,Ka,Sl_norm,localName = 'Unburbed gas'):
        self.Ka = Ka
        self.Sl_norm = Sl_norm
        self.localName = str(localName)

    def linear_fitting(self,Ka_threshold = 10):
        """ Linear Fitting  """
        print("\nMarkstein number: " + self.localName)
        self.x = []
        self.y = []
        for i in range(0,len(self.Ka)-1):
            if self.Ka[i] > Ka_threshold:
                break
            else:
                self.x.append(self.Ka[i])
                self.y.append(self.Sl_norm[i])        
        n = len(self.x)
        (ar, br) = np.polyfit(self.x, self.y, 1)
        xr = np.polyval([ar, br], self.x)
        # compute the mean square error
        err = np.sqrt(sum((xr - self.y)**2)/n)
        print('Linear regression:')
        print('regression: a=%.3f b=%.3f, ms error= %.3f' % ( ar, br, err))
        # Predict values
        self.x_fit = np.linspace(min(self.x),max(self.x),50)
        self.y_preditc = ar*self.x_fit + br
        print('Ma = {:.2f}\n'.format( -ar ) )
        
    def marks_line(self):
        return self.x_fit, self.y_preditc 


def main():

  data = pd.read_csv('./stretchResults/' + 'results.csv',index_col=None)
  
  phi = 1.
  phi = round(phi,1)
  D_th = flamesProp_CH4[phi]['D_th'] 
  Sl_o = flamesProp_CH4[phi]['Sl_o']
  deltaL = D_th/Sl_o
   
  # It limits to the linear stretch effect
  #threshold = 0.01
  threshold = 10 # For Sc
  
  # Burnt gas
  x = np.array(data['Kb']*deltaL/Sl_o)
  y = np.array(data['Su_b']/Sl_o)
  burntgas = MarksteinNumber(x,y,'Burnt gas') 
  burntgas.linear_fitting(Ka_threshold=threshold)

  # Unburnt gas
  x = np.array(data['Ku']*deltaL/Sl_o)
  y = np.array(data['Su']/Sl_o)
  #x = x[:-15]
  #y = y[:-15]
  unburntgas = MarksteinNumber(x,y,'Unburnt gas') 
  unburntgas.linear_fitting(Ka_threshold=threshold)

  # Consumption speed
  x = np.array(data['Kb']*deltaL/Sl_o)
  y = np.array(data['Sc']/Sl_o)
  consSpeed = MarksteinNumber(x,y,'Consumption speed') 
  consSpeed.linear_fitting(Ka_threshold=threshold)


  #--------------------------------------------------
  # PLOTS

  # Define fonts
  font = {'family': 'serif'}

  fontLegend = {'family': 'serif',
          'weight': 'normal',
          'size': 14,
          }

  fontLegend2 = {'family': 'serif',
          'weight': 'normal',
          'size': 14,
          }

  fontLegend3 = {'family': 'serif',
          'weight': 'normal',
          'size': 11,
          }
  #----------------------------------------------------------------------
  # Define fonts of plots [Math symbols and expressions]
  fonts1 = ["serif"]
  fonts2 = ["stix","stixsans","cm"]
  for font1,font2 in zip(fonts1,fonts2):
      plt.rcParams["font.family"] = font1
      plt.rcParams["mathtext.fontset"] = font2
      plt.rcParams["font.size"] = 11

  #----------------------------------------------------------

  fig = plt.figure(1,figsize=(6,5),facecolor='w')
  plt.scatter(data['Kb']*deltaL/Sl_o,data['Su_b']/Sl_o,s=6,color='k', lw =3)
  plt.plot(burntgas.marks_line()[0],burntgas.marks_line()[1],lw=2)
  plt.ylabel('$\\rm \\tilde{S}_d/S^0_l$',fontsize=17, fontdict=font)
  plt.xlabel('$\\rm Ka$',fontsize=17, fontdict=font)
  #plt.legend(loc='upper right',framealpha = 1,edgecolor='k',prop=fontLegend3)           
  plt.ylim([0.,1.2])
  plt.xlim([0,1.2])
  plt.tight_layout()
  plt.grid(color='gray', linestyle='--',lw=.5)

  fig = plt.figure(2,figsize=(6,5),facecolor='w')
  plt.scatter(data['Kb'],data['Su_b']/Sl_o,s=6,color='k', lw =3)
  plt.ylabel('$\\rm \\tilde{S}_d/S^0_l$',fontsize=17, fontdict=font)
  plt.xlabel('$\\rm Stretch\ rate\ (s^{-1})$',fontsize=17, fontdict=font)
  #plt.legend(loc='upper right',framealpha = 1,edgecolor='k',prop=fontLegend3)           
  plt.ylim([0.,1.2])
  #plt.xlim([0,350])
  plt.tight_layout()
  plt.grid(color='gray', linestyle='--',lw=.5)

  fig = plt.figure(3,figsize=(6,5),facecolor='w')
  plt.scatter(data['Ku']*deltaL/Sl_o,data['Su']/Sl_o,s=6,color='k', lw =3)
  plt.plot(unburntgas.marks_line()[0],unburntgas.marks_line()[1],lw=2)
  plt.ylabel('$\\rm {S}_u/S^0_l$',fontsize=17, fontdict=font)
  plt.xlabel('$\\rm Ka$',fontsize=17, fontdict=font)
  #plt.legend(loc='lower left',framealpha = 1,edgecolor='k',prop=fontLegend3)           
  plt.ylim([0.,1.2])
  plt.xlim([0,1.1])
  plt.tight_layout()
  plt.grid(color='gray', linestyle='--',lw=.5)

  fig = plt.figure(4,figsize=(6,5),facecolor='w')
  plt.scatter(data['Kb']*deltaL/Sl_o,data['Sc']/Sl_o,s=6,color='k', lw =3)
  plt.plot(consSpeed.marks_line()[0],consSpeed.marks_line()[1],lw=2)
  plt.ylabel('$\\rm {S}_c/S^0_l$',fontsize=17, fontdict=font)
  plt.xlabel('$\\rm Ka$',fontsize=17, fontdict=font)
  #plt.legend(loc='lower left',framealpha = 1,edgecolor='k',prop=fontLegend3)           
  plt.ylim([0.,1.2])
  plt.xlim([0,1.1])
  plt.tight_layout()
  plt.grid(color='gray', linestyle='--',lw=.5)
  plt.show()


if __name__=='__main__':
  print('\n\n*** Postprocessing Markstein number ***\n')
  main()
