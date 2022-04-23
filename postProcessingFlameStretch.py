#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import pdb
import time

from StretchRate import*

def main():

  phi = 1.
  path = './counterFlowResults/CH4/' + '{:.2f}'.format(phi) + '/'
  print('Data source: ' + path)
  pathToSave = './stretchResults/'
  if not os.path.isdir(pathToSave):
      os.makedirs(pathToSave)
      print("create folder : ", pathToSave)
  else:
      print(pathToSave, "folder already exists.")
  print('Path to Save: ' + pathToSave) 
  time.sleep(6) 
  
  uDatas = os.listdir(path)
  uDatas = sorted(uDatas, key=float)

  ulist = []
  Sulist = []
  Kulist = []
  Su_blist = []
  Kblist = []
  Sclist = []

  for uData in uDatas[1:-1]:  

    data = pd.read_csv( path + uData ,index_col=None)  

    # wdotO2 defines the point to take the flame speed
    refPoint = data['wdotO2'].argmin()

    Ku = np.gradient(data['velocity'],data['x'])
    maxStrLocation = abs(Ku).argmax()
    minVelocityPoint = data['velocity'][:maxStrLocation].argmin()

    # S_mass = m(x)/rho_u
    S_mass_u = data['rho']*data['velocity']/data['rho'][0]

    # Characteristic Strain Rate = K
    strainRatePoint = abs(Ku[:minVelocityPoint]).argmax()

    ph = PreHeat(data['wdotO2'],data['x'],data['velocity'],data['rho'])
    Su = ph.flame_speed_u()
    Ku_local =  ph.strain_rate_u()

    # Characteristic Flame Speed Su
    S_extrap = data['velocity'][strainRatePoint] - Ku_local*(data['x'][strainRatePoint:]-data['x'][strainRatePoint])

    print('** Pre-heat zone')
    print('Flame speed: {:.3f}'.format(Su))
    print('Straint rate: {:.2f}'.format(Ku_local))
    print('\n')

    # S_mass = m(x)/rho_b
    S_mass_b = data['rho']*data['velocity']/data['rho'][data['rho'].size-1]

    """ Kb = 1/rho*(grad(rho*u)) """
    grad = np.gradient(data['rho'][refPoint:-1]*data['velocity'][refPoint:-1],data['x'][refPoint:-1])
    Kb = np.array(-(1./(data['rho'][refPoint:-1]))*grad)

    rz = Reaction(data['wdotO2'],data['x'],data['velocity'],data['rho'])
    Sl_d = rz.flame_speed_b()
    Kb_local =  rz.strain_rate_b()

    print('** Reaction zone')
    print('Flame speed: {:.3f}'.format(Sl_d))
    print('Straint rate: {:.2f}\n'.format(Kb_local))

    fs = FlameSpeeds(data['x'],data['rho'],data['CH4'],data['wdotCH4'])
    Sc = fs.consumption_speed()
    
    print('** Global')
    print('Flame consumption speed: {:.3f}\n'.format(Sc))           

    ulist.append(uData)
    Sulist.append(Su)
    Kulist.append(Ku_local)
    Su_blist.append(Sl_d)
    Kblist.append(Kb_local)
    Sclist.append(Sc)


  df = pd.DataFrame()
  pd.options.display.float_format = '{:.6f}'.format  
  df['u'] = ulist
  df['Su'] = Sulist
  df['Ku'] = Kulist
  df['Su_b'] = Su_blist
  df['Kb'] = Kblist
  df['Sc'] = Sclist
  df.to_csv( pathToSave + 'results.csv', index = False)

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
  # PLOT STRETCH RATE AND FLAME SPEED AT PRE-HEAT ZONE
  fig = plt.figure(1,figsize=(7,5),facecolor='w')
  # Axial Velocity Plot
  L1 = plt.plot(data['x'],S_mass_u, 'magenta', lw=2, label=r'$m(x)/ \rho_u$')
  L2 = plt.plot(data['x'], data['velocity'], 'r', lw=2, label=r'$u$')
  L3 = plt.plot(data['x'][strainRatePoint:], S_extrap, 'k', linestyle = '--',lw=3, label=r'$S_u\ extrapolation$')
  plt.xlim(data['x'][0], data['x'][data['x'].size-1])

  # Identify the point where the strain rate is calculated
  plt.plot(data['x'][strainRatePoint], data['velocity'][strainRatePoint], 'gs')
  plt.annotate('Strain-Rate point',
               xy=(data['x'][strainRatePoint],
                   data['velocity'][strainRatePoint]),
               xytext=(0.001, 0.1),
               arrowprops={'arrowstyle': '->'})


  # Identify the point where the strain rate is calculated
  plt.plot(data['x'][refPoint],Su , 'gs')
  plt.annotate('Su point',
               xy=(data['x'][refPoint],Su),
               xytext=(0.01, 0.5),
               arrowprops={'arrowstyle': '->'})
  plt.xlabel('x axis (m)')
  plt.ylabel('Axial Velocity (m/s)')
  plt.twinx()
  L4 = plt.plot(data['x'], data['Qdot'], 'b', lw=2, label=r'$HRR$')
  plt.ylabel('HRR (Watt/m^3)')
  plt.axvline(x=data['x'][refPoint])
  plt.legend(L1+L2+L3+L4,[line.get_label() for line in L1+L2+L3+L4], \
             loc='upper right',framealpha = 1,edgecolor='k',prop=fontLegend3) 
  plt.tight_layout()

  fig = plt.figure(2,figsize=(7,5),facecolor='w')
  # Axial Velocity Plot
  L1 = plt.plot(data['x'], data['velocity'], 'r', lw=2,label=r'$u$')
  plt.xlim(data['x'][0], data['x'][data['x'].size-1])
  plt.xlabel('x axis (m)')
  plt.ylabel('Axial Velocity (m/s)')
  plt.twinx()
  L2 = plt.plot(data['x'], -Ku, 'b', lw=2,label=r'$K_u$')
  plt.ylabel('Stretch (s^-1)')
  plt.ylim([0,100])
  # Identify the point where the strain rate is calculated
  plt.plot(data['x'][strainRatePoint], Ku_local, 'gs')
  plt.annotate('Strain-Rate point',
               xy=(data['x'][strainRatePoint],Ku_local),
               xytext=(0.001, 60),
               arrowprops={'arrowstyle': '->'})
  plt.axvline(x=data['x'][refPoint])
  plt.legend(L1+L2,[line.get_label() for line in L1+L2], \
             loc='upper right',framealpha = 1,edgecolor='k',prop=fontLegend3) 
  plt.tight_layout()
  plt.show()

  # PLOT STRETCH RATE AND FLAME SPEED AT REACTION ZONE
  fig = plt.figure(1,figsize=(7,5),facecolor='w')
  # Axial Velocity Plot
  plt.title('Reaction zone')
  L1 = plt.plot(data['x'],S_mass_b, 'magenta', lw=2,label=r'$m(x)/\rho_b$')
  L2 = plt.plot(data['x'], data['velocity'], 'r', lw=2,label=r'$u$')
  L3 = plt.plot(data['x'][refPoint:],S_mass_b[refPoint:], 'k', lw = 1.5,label=r'$S_b\ extrapolation$')
  L4 = plt.plot(data['x'],data['rho']*data['velocity']/data['rho'][0], 'b',lw = 1.5,label=r'$S_u$')
  plt.xlim(data['x'][0], data['x'][data['x'].size-1])
  plt.xlabel('Distance (m)')
  plt.ylabel('Axial Velocity (m/s)')
  plt.twinx()
  plt.ylabel('HRR (Watt/m^3)')
  plt.legend(L1+L2+L3+L4,[line.get_label() for line in L1+L2+L3+L4], \
             loc='upper right',framealpha = 1,edgecolor='k',prop=fontLegend3) 
  plt.tight_layout()

  fig = plt.figure(2,figsize=(7,5),facecolor='w')
  plt.title('Reaction zone')
  L1 = plt.plot(data['x'], data['velocity'], 'r', lw=2,label=r'$u$')
  plt.xlim(data['x'][0], data['x'][data['x'].size-1])
  plt.xlabel('x axis (m)')
  plt.ylabel('Axial Velocity (m/s)')
  plt.twinx()
  L2 = plt.plot(data['x'][refPoint:-1], Kb, 'b', lw=2,label=r'$K_b$')
  plt.ylabel('Stretch (s^-1)')
  # Identify the point where the strain rate is calculated
  plt.plot(data['x'][refPoint], Kb_local, 'gs')
  plt.annotate('Strain-Rate point',
               xy=(data['x'][refPoint],
                   Kb_local),
               xytext=(0.001, 0.1),
               arrowprops={'arrowstyle': '->'})
  plt.axvline(x=data['x'][refPoint])
  plt.legend(L1+L2,[line.get_label() for line in L1+L2],\
             loc='upper right',framealpha = 1,edgecolor='k',prop=fontLegend3) 
  plt.tight_layout()
  plt.show()


if __name__=='__main__':
    print('\n*** Postprocessing the flame speed and stretch rata ***\n')
    main()
