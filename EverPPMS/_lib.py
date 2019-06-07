import matplotlib.colors as mcolors

def generate_FORC_sequence(
  sequence_filename,
  onsite_storage_folder,
  savefile_name,
  saturation_field,
  field_step,
  sweep_rate=25,
  averaging_time=1
  ):
  saturation_field = int(saturation_field)
  field_step = int(field_step)
  sweep_rate = int(sweep_rate)
  averaging_time = int(averaging_time)
  N = int(2*saturation_field/field_step)
  print(f'Generating sequence for {N} FORC measurements.')

  f = open(sequence_filename+'_0.seq', 'w')
  f.write('WAI WAITFOR 5 0 0 0 0 0\n')
  f.write(f'VSMDF "{onsite_storage_folder}\{savefile_name}" 0 0 ""\n')
  f.write(f'VSMLS 1 0 0 0 0 0\n')
  f.write('VSMCM "Saturate sample"\n')
  f.write(f'FLD FIELD {saturation_field} 100.0 0 0\n')
  f.write('WAI WAITFOR 5 0 1 0 0 0\n')
  count_lines = 8
  count_seq = 0
  for i in range(1, N + 1):
    if count_lines > 100:
      count_seq += 1
      next_seq_filename = sequence_filename+'_'+str(count_seq)+'.seq'
      f.write(f'CHN {onsite_storage_folder}\{next_seq_filename}')
      f.close()
      f = open(next_seq_filename, 'w')
      f.write(f'VSMDF "{onsite_storage_folder}\{savefile_name}" 0 1 ""\n')
      count_lines = 1
    min_field = saturation_field - i*field_step
    f.write(f'VSMCM "Measure step wise from {min_field} to {saturation_field}"\n')
    f.write(f'FLD FIELD {min_field} {sweep_rate} 0 1\n')
    f.write(f'WAI WAITFOR 1 0 1 0 0 0\n')
    f.write(f'VSMMH 1 36208647 0 1 0 2 40 {averaging_time} 0 2 0 1 2 -90000 {min_field} {saturation_field} {sweep_rate} 0 2 {field_step} 0 1 0 0 1 0 "A/C,0,10,10,0" "Q/M,0," \n')
    count_lines += 4
  f.write('FLD FIELD 0.0 100.0 0 0\n')
  f.close()

def generate_IRM_DCD_sequence(
  sequence_filename,
  onsite_storage_folder,
  savefile_name,
  saturation_field,
  field_step,
  sweep_rate=25,
  averaging_time=1
  ):
  saturation_field = int(saturation_field)
  field_step = int(field_step)
  N = int(saturation_field/field_step)

  sweep_rate = int(sweep_rate)
  averaging_time = int(averaging_time)
  print(f'Generating sequence for {N} IRM and DCD measurements.')

  f = open(sequence_filename+'_0.seq', 'w')
  f.write('WAI WAITFOR 5 0 0 0 0 0\n')
  f.write(f'VSMDF "{onsite_storage_folder}\{savefile_name}" 0 0 ""\n')
  f.write(f'VSMLS 1 0 0 0 0 0\n') # touchdown
  f.write(f'VSMCM "IRM Measurement"\n')
  count_lines = 4
  count_seq = 0

  # IRM loop
  for i in range(1, N + 1):
    if count_lines > 100:
      count_seq += 1
      next_seq_filename = sequence_filename+'_'+str(count_seq)+'.seq'
      f.write(f'CHN {onsite_storage_folder}\{next_seq_filename}')
      f.close()
      f = open(next_seq_filename, 'w')
      f.write(f'VSMDF "{onsite_storage_folder}\{savefile_name}" 0 1 ""\n')
      count_lines = 1
    next_field = i*field_step
    # go to next field and measure once
    f.write(f'FLD FIELD {next_field} {sweep_rate} 0 1\n')
    f.write(f'WAI WAITFOR 1 0 1 0 0 0\n')
    f.write(f'VSMCO 200 36208647 0 1 0 2 40 {averaging_time} 0 2 0 "A/C,0,10,10,0" "Q/M,0,"\n')
    # go to zero and measure once
    f.write(f'FLD FIELD 0.0 {sweep_rate} 0 1\n')
    f.write(f'WAI WAITFOR 1 0 1 0 0 0\n')
    f.write(f'VSMCO 200 36208647 0 1 0 2 40 {averaging_time} 0 2 0 "A/C,0,10,10,0" "Q/M,0,"\n')
    count_lines += 6

  f.write(f'VSMCM "DCD Measurement"\n')
  count_lines += 1
  # DCD loop
  for i in range(1, N + 1):
    if count_lines > 100:
      count_seq += 1
      next_seq_filename = sequence_filename+'_'+str(count_seq)+'.seq'
      f.write(f'CHN {onsite_storage_folder}\{next_seq_filename}')
      f.close()
      f = open(next_seq_filename, 'w')
      f.write(f'VSMDF "{onsite_storage_folder}\{savefile_name}" 0 1 ""\n')
      count_lines = 1
    next_field = -i*field_step
    # go to next field and measure once
    f.write(f'FLD FIELD {next_field} {sweep_rate} 0 1\n')
    f.write(f'WAI WAITFOR 1 0 1 0 0 0\n')
    f.write(f'VSMCO 200 36208647 0 1 0 2 40 {averaging_time} 0 2 0 "A/C,0,10,10,0" "Q/M,0,"\n')
    # go to zero and measure once
    f.write(f'FLD FIELD 0.0 {sweep_rate} 0 1\n')
    f.write(f'WAI WAITFOR 1 0 1 0 0 0\n')
    f.write(f'VSMCO 200 36208647 0 1 0 2 40 {averaging_time} 0 2 0 "A/C,0,10,10,0" "Q/M,0,"\n')
    count_lines += 6

  # go to zero field and finish
  f.write('FLD FIELD 0.0 100.0 0 0\n')
  f.close()

def get_cmap():
  c = mcolors.ColorConverter().to_rgb
  custom_colors = [(0, 0, 0, 0),\
            (0.18, 0.05, 0.05, 0.2),\
            (0.28, 0, 0, 1),\
            (0.4, 0.7, 0.85, 0.9),\
            (0.45, 0, 0.75, 0),\
            (0.6, 1, 1, 0),\
            (0.75, 1, 0, 0),\
            (0.92 , 0.6, 0.6, 0.6),\
            (1  , 0.95, 0.95, 0.95)]
  cdict = {'red': [], 'green': [], 'blue': []}
  for i, item in enumerate(custom_colors):
      pos, r, g, b = item
      cdict['red'].append([pos, r, r])
      cdict['green'].append([pos, g, g])
      cdict['blue'].append([pos, b, b])
  cmap =  mcolors.LinearSegmentedColormap('CustomMap', cdict)
  cmap.set_bad(color='black')
  return cmap