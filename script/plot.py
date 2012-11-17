#! /usr/bin/python
# Filename: plot.py

import matplotlib
import matplotlib.pyplot as plt  
import time 

data_map = {}
lst_nr = [2000, 5000, 8000, 10000, 15000, 20000, 30000, 50000, 100000]
lst_r0 = [(i + 1) * 10 for i in xrange(10)]
lst_radio = [i + 1 for i in xrange(10)]

def gen_fig_single(x_lst, y_lst, title, x_label, y_label):
  fig = plt.figure()
  plt.xlabel(x_label)
  plt.ylabel(y_label)
  plt.title(title)
  plt.plot(x_lst, y_lst)
  plt.legend()
  return fig
 
def it_gen_fig_single(dem1, dem2, dem3, lst1, lst2, lst3, path):
  x_label = dem3
  y_label = 'relative_error'
  for var1 in lst1:
    for var2 in lst2:
      x_lst = lst3
      y_lst = []
      for var3 in lst3:
        if dem3 == 'radio':   
          indx = (var1, var2, var3)
        if dem3 == 'r0':
          indx = (var1, var3, var2)
        if dem3 == 'nr':
          indx = (var3, var1, var2)
        y_lst.append(data_map[indx][0]) 
      title = dem1 + ': ' + str(var1) + ' ' + dem2 + ': ' + str(var2)
      fig = gen_fig_single(x_lst, y_lst, title, x_label, y_label)
      filename = 'single_plot_' + dem1 + '_' + dem2 + '_' + str(var1) + '_' + str(var2)
      filename = path + filename
      fig.savefig(filename)
  return 

def it_gen_fig_multi(dem1, dem2, dem3, lst1, lst2, lst3, path):
  x_label = dem3
  y_label = 'relative_error'
  for var1 in lst1:
    title = dem1 + ': ' + str(var1)
    fig = plt.figure()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    for var2 in lst2:
      x_lst = lst3
      y_lst = []
      for var3 in lst3:
        if dem3 == 'radio' and dem2 == 'r0':   
          indx = (var1, var2, var3)
        if dem3 == 'radio' and dem2 == 'nr':
          indx = (var2, var1, var3)
        if dem3 == 'nr' and dem2 == 'radio':
          indx = (var3, var1, var2)
        if dem3 == 'nr' and dem2 == 'r0':   
          indx = (var3, var2, var1)
        if dem3 == 'r0' and dem2 == 'nr':   
          indx = (var2, var3, var1)
        if dem3 == 'r0' and dem2 == 'radio':   
          indx = (var1, var3, var2)
        y_lst.append(data_map[indx][0])
      curve_label = dem2 + ': ' + str(var2)
      fig.hold(True)
      plt.plot(x_lst, y_lst, label = curve_label)
      fig.hold(True)
      plt.legend()
      fig.hold(True)
    filename = 'multi_plot_' + dem1 + '_' + str(var1)
    filename = path + filename
    fig.savefig(filename)
  return 


if __name__ == '__main__':

  f = file('../result.txt', 'r')
  for line in f:
    line = line.strip('\n')
    lst = line.split('\t')
    key = (float(lst[2]), float(lst[3]), float(lst[4]))
    value = (float(lst[6]), float(lst[5]), float(lst[1]), float(lst[0]))
    data_map[key] = value
  
  # plot radio_para <--> r_err 
  #it_gen_fig_single('nr', 'r0', 'radio', lst_nr, lst_r0, lst_radio, './radio/')
  it_gen_fig_multi('nr', 'r0', 'radio', lst_nr, lst_r0, lst_radio, './multi_nr_r0/')
  it_gen_fig_multi('nr', 'radio', 'r0', lst_nr, lst_radio, lst_r0, './multi_nr_radio/')
  
  time.sleep(10)
  
  # plot radio_para <--> r_err 
  #it_gen_fig_single('nr', 'radio', 'r0', lst_nr, lst_radio, lst_r0, './r0/')
  it_gen_fig_multi('r0', 'radio', 'nr', lst_r0, lst_radio, lst_nr, './multi_r0_radio/')
  it_gen_fig_multi('r0', 'nr', 'radio', lst_r0, lst_nr, lst_radio, './multi_r0_nr/')
  
  time.sleep(10)

  # plot nr <--> r_err 
  #it_gen_fig_single('r0', 'radio', 'nr', lst_r0, lst_radio, lst_nr, './nr/')
  it_gen_fig_multi('radio', 'nr', 'r0', lst_radio, lst_nr, lst_r0, './multi_radio_nr/')
  it_gen_fig_multi('radio', 'r0', 'nr', lst_radio, lst_r0, lst_nr, './multi_radio_r0/')

  time.sleep(10)

  # the end of the script!
