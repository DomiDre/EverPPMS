import numpy as np
import lmfit, sys

from . import closest_idx

class IRMDCD:
  def __init__(self, verbose=False):
    self.log = ''
    self.B_irm = None
    self.M_irm = None
    self.sM_irm = None
    self.B_dcd = None
    self.M_dcd = None
    self.sM_dcd = None

    self.verbose = verbose

    self._Ms = 0

    self.calculated_deltaM = False

  def printlog(self, message):
    print(message)
    self.log = self.log + message + '\n'

  def load_dat(self, filename, markerIRM = 'IRM MEASUREMENT',
  markerDCD = 'DCD MEASUREMENT'):
    self.printlog(f'Opening {filename}')
    self.printlog(f'IRM marker: {markerIRM}')
    self.printlog(f'DCD marker: {markerDCD}')

    B_irm = []
    M_irm = []
    sM_irm = []

    B_dcd = []
    M_dcd = []
    sM_dcd = []

    with open(filename, "r", errors='ignore') as f:

      line_counter = 0
      is_data_line = False
      # skip all comments line until [Data] is found in .DAT file
      while not is_data_line:
        line = f.readline()
        line_counter += 1
        if '[Data]' in line:
          is_data_line = True

      # after [Data] follows header
      # determine idx of relevant measurement parameters
      header = f.readline().rstrip().split(',')
      line_counter += 1
      column_name_to_idx = {}
      for i, column_name in enumerate(header):
        column_name_to_idx[column_name] = i

      idx_T = column_name_to_idx['Temperature (K)']
      idx_B = column_name_to_idx['Magnetic Field (Oe)']
      idx_M = column_name_to_idx['Moment (emu)']
      idx_sM = column_name_to_idx['M. Std. Err. (emu)']

      # create pointer for either IRM oder DCD measurement arrays
      helper_H, helper_M, helper_sM = None, None, None

      # read all following data lines
      for line in f:
        line_counter += 1
        # remove unnecessary whitespaces
        linestrip = line.strip()

        # check for any markers in the line
        if markerIRM in line:
          helper_H, helper_M, helper_sM = B_irm, M_irm, sM_irm
          continue
        if markerDCD in line:
          helper_H, helper_M, helper_sM = B_dcd, M_dcd, sM_dcd
          continue

        # ignore comments
        if linestrip.startswith("#") or linestrip == "":
          continue

        # if no marker was set yet, continue as you are before any measurement
        if helper_H is None:
          self.printlog(f'Line {line_counter} after header was skipped.')
          continue

        # a marker is set
        # check if the current line is valid as that it set the sample at a finite magnetic field


        #seperate line by commas into lists
        split_line = linestrip.split(',')
        B_line = np.round(float(split_line[idx_B]) / 1e4, 3)
        if np.abs(B_line) < 1e-3: # smaller than 1 mT?
          self.printlog(f"Set field in line {line_counter} is too small in magnitude. Skipping line.")
          continue

        # Sample was prepared at magnetic field. Check next line whether it is a remanent state
        next_line = f.readline().strip()
        line_counter += 1
        if next_line.startswith("#") or next_line == "":
          self.printlog(f'A supposed pair at line {line_counter} of IRM/DCD is a comment or empty! Skipping pair in total.')
          continue

        split_next_line = next_line.split(',')

        B_remanence = np.round(float(split_next_line[idx_B]) / 1e4, 3)

        if np.abs(B_remanence) > 1e-3: # remanent field should be smaller than 1 mT
          self.printlog(f'Remanence measurement in line {line_counter} after header is not zero. Skipping pair in total.')
          continue

        # found a valid pair. Store it
        helper_H.append(B_line)
        helper_M.append(float(split_next_line[idx_M]) * 1e3)
        helper_sM.append(float(split_next_line[idx_sM]) * 1e3)

    self.B_irm = np.array(B_irm)
    self.M_irm = np.array(M_irm)
    self.sM_irm = np.array(sM_irm)

    self.B_dcd = np.array(B_dcd)
    self.M_dcd = np.array(M_dcd)
    self.sM_dcd = np.array(sM_dcd)

    return self.B_irm, self.M_irm, self.sM_irm, self.B_dcd, self.M_dcd, self.sM_dcd

  @property
  def Ms(self):
    return self._Ms

  @Ms.setter
  def Ms(self, _Ms):
    self.printlog(f'Setting saturation magnetization to {_Ms}')
    self._Ms = _Ms

  def calcDeltaM(self):
    B_deltaM = []
    deltaM = []
    sDeltaM = []
    for idxIrm, B_irm in enumerate(self.B_irm):
      idxDcd = closest_idx(-self.B_dcd, B_irm)
      if abs(self.B_irm[idxIrm] + self.B_dcd[idxDcd]) < 1e-3:
        B_deltaM.append(self.B_irm[idxIrm])
        deltaM.append(self.M_dcd[idxDcd] - (self.Ms - 2*self.M_irm[idxIrm]))
        sDeltaM.append(np.sqrt(self.sM_dcd[idxDcd]**2 + 4*self.sM_irm[idxIrm]**2))

    self.B_deltaM = np.array(B_deltaM)
    self.deltaM = np.array(deltaM)
    self.sDeltaM = np.array(sDeltaM)
    self.calculated_deltaM = True
    return B_deltaM, deltaM, sDeltaM

  def export(self, filename):
    if len(filename.rsplit('.', 1)) >= 2:
      raw_filename = filename.rsplit('.', 1)[0] + '_raw.dat'
    else:
      raw_filename = filename + '_raw.dat'
      filename = filename + '.dat'

    with open(raw_filename, 'w') as f:
      f.write('#' + self.log.replace('\n','\n#') + '\n')

      f.write('\n#IRM\n')
      f.write('#B / T\tM / memu\tsM / memu\n')
      for idx_irm in range(len(self.B_irm)):
        f.write(f'{self.B_irm[idx_irm]}\t{self.M_irm[idx_irm]}\t{self.sM_irm[idx_irm]}\n')

      f.write('\n#DCD\n')
      f.write('#B / T\tM / memu\tsM / memu\n')
      for idx_dcd in range(len(self.B_dcd)):
        f.write(f'{self.B_dcd[idx_dcd]}\t{self.M_dcd[idx_dcd]}\t{self.sM_dcd[idx_dcd]}\n')

    if self.calculated_deltaM:
      with open(filename, 'w') as f:
        f.write('#' + self.log.replace('\n','\n#') + '\n')

        f.write('#deltaM\n')
        f.write('#B / T\tdeltaM / memu\tsDeltaM / memu\n')
        for idx_delta in range(len(self.B_deltaM)):
          f.write(f'{self.B_deltaM[idx_delta]}\t{self.deltaM[idx_delta]}\t{self.sDeltaM[idx_delta]}\n')
