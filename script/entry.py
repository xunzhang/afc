#! /usr/bin/python
# Filename: entry.py

import os
import time
import math

dicts = {}
lst_nr = ['2000', '5000', '8000', '10000', '15000', '20000', '30000', '50000', '100000']
lst_r0 = [(i + 1) * 10 for i in xrange(10)]
lst_radio = [i + 1 for i in xrange(10)]

def sphere_vol(radius):
  return (4 * math.pi * radius ** 3) / 3

def script_entry(f):
  for nr in lst_nr:
    for r0 in lst_r0:
      volumn = sphere_vol(r0)
      for radio in lst_radio:
        indx = (nr, r0, radio)
      
        command = './a.out ' + nr + ' ' + str(r0) + ' ' + str(radio)

        # value = os.system(command)
        value_lst = os.popen(command).readlines()
        value = float(value_lst[0].strip(' '))

        diff = value - volumn
        r_err = math.fabs(diff) / volumn
        print math.fabs(diff), ' and ', volumn, 'and', math.fabs(diff)/volumn

        dicts[indx] = (value, r_err)
      
        para = float(r0) / float(radio)
        print '------------------------------------------------------------------------'
        print 'para:', para
        print 'std:', volumn
        print command
        print indx, '-->', dicts[indx] 
        print '------------------------------------------------------------------------'
	content = str(para) + '\t' + str(volumn) + '\t' + nr + '\t' + str(r0) + '\t' + str(radio) + '\t' + str(value) + '\t' + str(r_err)
	f.write(content+'\n')
        time.sleep(20)
  return
  
if __name__ == '__main__':
  f = file('data.txt', 'w+')
  script_entry(f)
  #print dicts
  f.close()
