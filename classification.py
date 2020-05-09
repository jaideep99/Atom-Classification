import pandas as pd
import numpy as np
import re

#dictionary with atom properties
pm = {'Ag': {'mp': 961, 'enthalphy': 11.3, 'en': 1.93, 'ar': 1.72}, 'Al': {'mp': 660.25, 'enthalphy': 10.7, 'en': 1.61, 'ar': 1.43}, 'Au': {'mp': 1064, 'enthalphy': 12.5, 'en': 2.54, 'ar': 1.66}, 'B': {'mp': 2076, 'enthalphy': 50.0, 'en': 2.04, 'ar': 0.85}, 'Be': {'mp': 1287, 'enthalphy': 7.95, 'en': 1.57, 'ar': 1.12}, 'C': {'mp': 3500, 'enthalphy': 105.0, 'en': 2.55, 'ar': 0.7}, 'Ca': {'mp': 839, 'enthalphy': 8.53, 'en': 1.0, 'ar': 2.31}, 'Ce': {'mp': 798, 'enthalphy': 5.5, 'en': 1.12, 'ar': 2.48}, 'Co': {'mp': 1495, 'enthalphy': 16.2, 'en': 1.88, 'ar': 2.0}, 'Cr': {'mp': 1857, 'enthalphy': 20.5, 'en': 1.66, 'ar': 1.28}, 'Cu': {'mp': 1084.6, 'enthalphy': 13.1, 'en': 1.9, 'ar': 1.4}, 'Dy': {'mp': 1412, 'enthalphy': 11.1, 'en': 1.22, 'ar': 2.35}, 'Er': {'mp': 1522, 'enthalphy': 19.9, 'en': 1.24, 'ar': 2.32}, 'Fe': {'mp': 1535, 'enthalphy': 13.8, 'en': 1.83, 'ar': 1.26}, 'Ga': {'mp': 29.76, 'enthalphy': 5.59, 'en': 1.81, 'ar': 1.87}, 'Gd': {'mp': 1312, 'enthalphy': 10.0, 'en': 1.2, 'ar': 1.8}, 'Hf': {'mp': 2227, 'enthalphy': 25.5, 'en': 1.3, 'ar': 2.25}, 'Ho': {'mp': 1470, 'enthalphy': 17.0, 'en': 1.23, 'ar': 2.33}, 'In': {'mp': 156.76, 'enthalphy': 3.26, 'en': 1.78, 'ar': 1.56}, 'La': {'mp': 920, 'enthalphy': 6.3, 'en': 1.1, 'ar': 2.5}, 'Mg': {'mp': 649, 'enthalphy': 8.7, 'en': 1.31, 'ar': 1.73}, 'Mn': {'mp': 1244, 'enthalphy': 13.2, 'en': 1.55, 'ar': 1.61}, 'Mo': {'mp': 2617, 'enthalphy': 36.0, 'en': 2.16, 'ar': 1.9}, 'Nb': {'mp': 2468, 'enthalphy': 26.8, 'en': 1.6, 'ar': 2.15}, 'Nd': {'mp': 1016, 'enthalphy': 7.1, 'en': 1.14, 'ar': 2.29}, 'Ni': {'mp': 1453, 'enthalphy': 17.2, 'en': 1.91, 'ar': 1.63}, 'P': {'mp': 44.3, 'enthalphy': 0.64, 'en': 2.19, 'ar': 1.95}, 'Pd': {'mp': 1552, 'enthalphy': 16.0, 'en': 2.2, 'ar': 1.63}, 'Pr': {'mp': 931, 'enthalphy': 6.9, 'en': 1.13, 'ar': 2.47}, 'Sc': {'mp': 1541, 'enthalphy': 16.0, 'en': 1.36, 'ar': 2.3}, 'Si': {'mp': 1414, 'enthalphy': 50.2, 'en': 1.9, 'ar': 1.11}, 'Sn': {'mp': 231.9, 'enthalphy': 7.0, 'en': 1.96, 'ar': 2.25}, 'Ta': {'mp': 3017, 'enthalphy': 36.0, 'en': 1.5, 'ar': 2.2}, 'Ti': {'mp': 1668, 'enthalphy': 18.7, 'en': 1.54, 'ar': 1.47}, 'Tm': {'mp': 1545, 'enthalphy': 16.8, 'en': 1.25, 'ar': 2.3}, 'Y': {'mp': 1526, 'enthalphy': 11.4, 'en': 1.22, 'ar': 2.4}, 'Zn': {'mp': 419.73, 'enthalphy': 7.35, 'en': 1.65, 'ar': 1.39}, 'Zr': {'mp': 1852, 'enthalphy': 21.0, 'en': 1.33, 'ar': 2.3}}

def diff(x):
  #Seperating atoms from composition
  s = re.sub(r'[^\w\s]','',x)
  s = re.sub('\d',' ',s)
  x = np.array([i for i in s.split(' ') if i in pm])

  print('Elements in BMG are : ', x)
  
  #making ranges for each atom
  ranges = {}
  for i in x:
    ranges[i] = 0.88*pm[i]['ar']
  
  #compiling scoring matrix
  score = {}
  
  for i in x:
    dct = {}
    for j in x:
      if pm[i]['ar']<ranges[j]:
        dct[j] = -1
      elif pm[i]['ar']>pm[j]['ar']:
        dct[j] = 1
      else :
        dct[j] = 0
          
    score[i] = dct
    
    
  dt = pd.DataFrame(score)

  print('\nScoring Matrix : ')
  print(dt)
  print()

  b = []
  sm = []
  
  #separating into big and small based on scoring matrix
  for i in score:
    sum = 0
    for j in score[i]:
      sum = sum+score[i][j]
      
    if(sum>0):
      b.append(i)
      
    else : 
      sm.append(i)
      
      
  return b,sm

bmg = input('Enter the BMG : ')
big,small = diff(bmg)

print('Big Atoms : ',big)
print('Small Atoms : ',small)
