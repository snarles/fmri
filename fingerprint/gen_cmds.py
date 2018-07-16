for i in range(1,18):
  for j in range(i+1, 19):
    print('Rscript /data/MLcore/fmri/fingerprint/fingerprinting2batch_dc.R ' + str(i) + ' ' + str(j) + ' ' + str(8))
