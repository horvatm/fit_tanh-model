#!/usr/bin/env python
# -*- coding: utf8 -*-

import scipy.optimize 
import numpy as np
import math

def _sech2(x):
  return 1./math.cosh(x)**2 if x < 355 and x > -355 else 0.
  
sech2 = np.vectorize(_sech2)

'''
  The Tanh-model function
''' 
def model(x, p, norm = False):
  if norm:  # normalized model
    t = np.tanh(p[1])
    return 1 + p[2]*(np.tanh(p[0]*x + p[1]) - t) 
  else:     # un-normalized model
    return p[2]*np.tanh(p[0]*x + p[1]) + p[3]

'''
  Derivative of Tanh-model function w.r.t. parameters
'''
def dmodel_dpars(x, p, norm = False):
  s = sech2(p[0]*x + p[1])
  t = np.tanh(p[0]*x + p[1])
  
  if norm:  # normalized model
    s1 = sech2(p[1])
    t1 = np.tanh(p[1])
    return np.array([p[2]*x*s, p[2]*(s - s1), t - t1])
  else:     # un-normalized model
    return np.array([p[2]*x*s, p[2]*s, t, np.ones(len(x), dtype=float)])

'''
  Fitting model to the weather forecast uncentainty data
  
  Data:
  
    d = [[t_0 E_0], [t_1, E_1], ..., [t_{n-1}, E_{n-1}]]     
    
  Fitting model:
  
    E(t) = A tanh(a t  + b ) + B            for norm=False
    E(t) = 1 + c[tanh(at + b) − tanh(b)]    for norm=True
    
  Input:
    times:1D numpy array of floats
    E: 1D numpy array of floats
    norm: boolean
  
  Return:
    [a, b, c]         for norm=True
    [a, b, A, B]      for norm=False
'''

def fit_model(t, E, w = None, norm = False):
  
  x = np.array(t)
  y = np.array(E)
  n = len(x)
  
  # sorting data
  a = np.array([x,y]).T
  x, y = a[a[:,0].argsort()].T
  
  # rescale to normalized data
  scale = y[0]
  y /= scale
  
  ymax = np.amax(y)  
    
  if w is None:
    w = np.ones(n)
  
  def f(p):
    return np.sum(((y - model(x, p, norm))/w)**2)     
    
  def df(p):
    return np.dot(dmodel_dpars(x, p, norm), 2*(model(x, p, norm) - y)/w**2)
  
  min_kwargs = {"method": "BFGS", "jac": df, "tol": 1e-6}
  
  # init parameter: 
  #   assuming E ~ A tanh(a t) + B 
  #   x_0 ~ 0, y[0] = 1 
  p0 = np.array([(y[1] - 1)/(ymax - 1)/(x[1] - x[0]), 0., ymax - 1, 1.])
  T = f(p0)
  
  res = scipy.optimize.basinhopping(f, p0, T = T, minimizer_kwargs = min_kwargs, niter=200)
  
  p = res["x"]
  
  return {"err": (scale**2)*res["fun"], "p": np.concatenate((p[:2], scale*p[2:]))}
 
