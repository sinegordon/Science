#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess
from multiprocessing import Process
import os
import numpy as np

def f (program, in_filename, out_filename):
    subprocess.run([program, in_filename, out_filename], shell=True)

with open("./inTEMPLATE.txt") as f:
    tmp = f.read()

var_mas = {}

var_mas['tmin'] = 0.0 #По времени нижняя граница
var_mas['tmax'] = 100.0 #По времени верхняя граница
var_mas['xmin'] = -50.0 #По пространству нижняя граница
var_mas['xmax'] = 50.0 #По пространству верхняя граница
var_mas['nx'] = 20000 #По пространству число шагов
var_mas['threads'] = 1 #Количестов расчетных потоков
var_mas['b'] = 2.0 #b
var_mas['masnt'] = 100 #Число точек по оси времени в массиве
var_mas['divx'] = 20 #Делитель количества точек по оси координаты для сохранения в файл
var_mas['divt'] = 4000 #Делитель количества точек по оси времени для сохранения в файл
var_mas['v1'] = 0.5 #Скорость первого импульса
var_mas['x10'] = 0.0 #Координата центра первого начального импульса
var_mas['v2'] = 0.0 #Скорость второго импульса
var_mas['x20'] = 0.0 #Координата центра второго начального импульса
var_mas['xb'] = 0.0 #Координата дельта-барьера
var_mas['mu'] = 0.0 #Мощность дельта-барьера
var_mas['a'] = 0.0 #Амплитуда ВЧ поля
var_mas['beta'] = 0.02 #Коэффициент трения
var_mas['w'] = 10.0 #Частота поля
var_mas['tau'] = 5.0 #Характерное время включения поля

if(not os.path.exists("./res/")):
    os.makedirs("./res/")

count = 1
for a in [0.0, 1.0, 2.0, 3.0]:
    s = tmp
    var_mas['a'] = a
    for k, v in var_mas.items():
        s = s.replace("{" + k + "}", "{val}".format(val = v))
    dir = os.getcwd()
    in_filename = "{dir}/res/in_{a}_{b}.txt".format(dir = dir, a = var_mas['a'], b = var_mas['b'])
    out_filename = "{dir}/res/out_{a}_{b}.txt".format(dir = dir, a = var_mas['a'], b = var_mas['b'])
    exe = "{dir}/a.out {in_filename} {out_filename}".format(dir = dir, in_filename = in_filename, out_filename = out_filename)
    sh = "{dir}/res/run{n}.sh".format(dir = dir, n = count)
    with open(sh, "w") as f:
        f.write(exe)
    with open(in_filename, "w") as f:
        f.write(s)
    proc = subprocess.Popen(["chmod", "777", sh])
    proc.wait()
    #proc = subprocess.Popen(["sh", sh])
    #proc = subprocess.Popen(["open", "-a", "Terminal.app", sh])
    proc = subprocess.Popen(['xterm', '-e', sh])
    #proc.wait()
    #proc = subprocess.Popen(["rm", "-f", sh])
    #proc.wait()
    count += 1
    #subprocess.call(['xterm', '-e', 'python bb.py'])