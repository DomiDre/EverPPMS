import numpy as np
import lmfit, sys

class Evaluator:
  def __init__(self, filename, H_sat, H_step, verbose=False):
    self.log = ''
    self.H_sat = H_sat
    self.H_step = H_step
    self.H_a = None
    self.H_b = None
    self.M = None
    self.sM = None
    self.verbose = verbose
    self.load_dat(filename)

  def printlog(self, message):
    print(message)
    self.log = self.log + message + '\n'

  def load_dat(self, filename, marker='MEASURE STEP WISE'):
    self.printlog(f'Opening {filename}')

    H_a = []
    H_b = []
    M = []
    sM = []
    with open(filename, "r", errors='ignore') as f:
      is_data_line = False
      while not is_data_line:
        line = f.readline()
        if '[Data]' in line:
          is_data_line = True
      header = f.readline().rstrip().split(',')
      column_name_to_idx = {}
      for i, column_name in enumerate(header):
        column_name_to_idx[column_name] = i

      idx_T = column_name_to_idx['Temperature (K)']
      idx_B = column_name_to_idx['Magnetic Field (Oe)']
      idx_M = column_name_to_idx['Moment (emu)']
      idx_sM = column_name_to_idx['M. Std. Err. (emu)']

      is_forc_line = False # skip possible sweep measurement at the beginning

      forc_line_H_b = []
      forc_line_M = []
      forc_line_sM = []
      idx_forcs = 0
      for line in f:
        linestrip = line.strip()
        if linestrip.startswith("#") or linestrip == "":
          continue

        is_new_marker = marker in line
        if not is_forc_line: # still sweep measurement?
          is_forc_line = is_new_marker
          continue

        if is_new_marker:
          idx_forcs += 1
          H_a.append(self.H_sat - idx_forcs*self.H_step)
          H_b.append(forc_line_H_b)
          M.append(forc_line_M)
          sM.append(forc_line_sM)
          forc_line_H_b = []
          forc_line_M = []
          forc_line_sM = []
          continue

        split_line = linestrip.split(',')
        forc_line_H_b.append(np.round(float(split_line[idx_B]) / 1e4, 3))
        forc_line_M.append(float(split_line[idx_M]) * 1e3)
        forc_line_sM.append(float(split_line[idx_sM]) * 1e3)
    if len(forc_line_H_b) > len(H_b[-1]):
      idx_forcs += 1
      H_a.append(self.H_sat - idx_forcs*self.H_step)
      H_b.append(forc_line_H_b)
      M.append(forc_line_M)
      sM.append(forc_line_sM)

    self.H_a_raw = np.array(H_a)
    self.H_b_raw = np.array(H_b)
    self.M_raw = np.array(M)
    self.sM_raw = np.array(sM)
    
    Na = len(H_a)
    Nb = len(H_b[-1])
    self.H_a = np.array(H_a)
    self.H_b = np.array(H_b[-1])
    self.M = np.zeros((Na, Nb))
    self.sM = np.zeros((Na, Nb))

    for i in range(Na):
      for j in range(len(H_b[i])):
        self.M[i, -1 - j] = M[i][-1 - j]
        self.sM[i, -1 - j] = sM[i][-1 - j]

    self.H_a = self.H_a[::-1]
    self.M = self.M[::-1, :]
    self.sM = self.sM[::-1, :]

  def polynomialFORC(self, Ha, Hb, a1, a2, a3, a4, a5, a6):
    return a1 + a2*Ha + a3*Ha**2 + a4*Hb + a5*Hb**2 + a6*Ha*Hb

  def residuum(self, p, Ha, Hb, M, sM):
    return (
      M - self.polynomialFORC(
        Ha, Hb, p['a1'], p['a2'], p['a3'], p['a4'], p['a5'], p['a6'])
      )/sM

  def calcFORCdistribution(self, smoothing_factor=2):
    FORCdistribution = np.zeros(self.M.shape)
    for i in range(2*smoothing_factor, self.M.shape[0] - 2*smoothing_factor):
      fit_Ha = self.H_a[i- smoothing_factor:i+ smoothing_factor]
      print(f"Fitting {i}/{self.M.shape[0]- smoothing_factor} \r", end='')
      sys.stdout.flush()
      for j in range( i + 2*smoothing_factor, self.M.shape[1] - 2*smoothing_factor):
        fit_Hb = self.H_b[j - smoothing_factor:j + smoothing_factor]
        grid_Ha, grid_Hb = np.meshgrid(fit_Ha, fit_Hb, indexing='ij')
        p_init = lmfit.Parameters()
        p_init.add('a1',
          np.mean(
            self.M[i - smoothing_factor:i + smoothing_factor,
                   j - smoothing_factor:j + smoothing_factor])
        )
        p_init.add('a2', 0)
        p_init.add('a3', 0)
        p_init.add('a4', 0)
        p_init.add('a5', 0)
        p_init.add('a6', 0)

        fit_result = lmfit.minimize(
          self.residuum, p_init,
          args=(
            grid_Ha, grid_Hb,
            self.M[i - smoothing_factor: i + smoothing_factor,
                   j - smoothing_factor: j + smoothing_factor],
            self.sM[i - smoothing_factor: i + smoothing_factor,
                    j - smoothing_factor: j + smoothing_factor]
        ))
        if self.verbose:
          print(lmfit.fit_report(fit_result))
        FORCdistribution[i,j] = -fit_result.params['a6']
    self.FORCdistribution = FORCdistribution
    print('Done.')

  def rotateFORC(self):
    # Rotate coordinate system
    half_step = np.round(self.H_step/2., 3)
    H_c = np.arange(
      np.round((min(self.H_b) - max(self.H_a))/2, 3),
      np.round((max(self.H_b) - min(self.H_a))/2, 3) + half_step/2.,
      half_step
    )
    H_u = np.arange(
      np.round((min(self.H_b) + min(self.H_a))/2, 3),
      np.round((max(self.H_b) + max(self.H_a))/2, 3) + half_step/2.,
      half_step
    )

    FORC = np.zeros((len(H_c), len(H_u)))
    for i, valH_a in enumerate(self.H_a):
      for j, valH_b in enumerate(self.H_b):
        valH_c = np.round((valH_b - valH_a)/2., 4)
        valH_u = np.round((valH_b + valH_a)/2., 4)
        k = int(np.round((valH_c - H_c[0])/half_step, 0))
        l = int(np.round((valH_u - H_u[0])/half_step, 0))
        FORC[k, l] = self.FORCdistribution[i, j]

    # fill the spacings in between
    for k in range(1, FORC.shape[0]-1):
      for l in range(1, FORC.shape[1]-1):
        if k % 2 == 0:
          if l % 2 > 0 :
            FORC[k, l] = (FORC[k+1, l] + FORC[k, l+1] + FORC[k-1, l] + FORC[k, l-1])/4
        else:
          if l % 2 == 0 :
            FORC[k, l] = (FORC[k+1, l] + FORC[k, l+1] + FORC[k-1, l] + FORC[k, l-1])/4


    self.H_c = H_c
    self.H_u = H_u
    self.FORC = FORC