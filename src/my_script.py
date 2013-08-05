def num2str4fig(num, scientific=True):
  import re
  s1 = '{0:16.2f}'.format(num).strip()
  s2 = [re.sub(r'^([\-+]*)(0+)([1-9]+)', r'\1\3', item) for item in \
        [item.replace('+','') for item in '{0:16.2e}'.format(num).strip().split('e')]]
  len1 = len(s1)
  len2 = 3
  for item in s2:
    len2 += len(item)
  choose = 1
  if (num != 0.0):
    if float(s1) == 0.0:
      choose = 2
  if (not scientific) or (len1 <= len2 and choose==1):
    return r'${' + s1 + r'}$'
  else:
    if float(s2[1]) == 0.0 or float(s2[1]) == 1.0:
      return r'${' + s1 + r'}$'
    else:
      return r'${' + s2[0] + r'}{\times}10^{' + s2[1] + r'}$'


def clip_fraction(num_str, th=0.01):
  num = float(num_str)
  num_round = round(num)
  if abs(num_round-num)/num < th:
    return num_str.split('.')[0]
  else:
    return num_str


def round2beauty(num):
  return num
