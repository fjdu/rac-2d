def load_data(filename, comments='!'):
  import numpy as np

  data = np.loadtxt(filename, comments)

  f = open(filename, 'r')
  str_comment = f.readline()[1:].split()
  f.close()

  dic = {}
  for i in xrange(len(str_comment)):
    dic.update({str_comment[i]: data[:, i]})

  return dic
