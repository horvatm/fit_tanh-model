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
  Fiting data:
  
    d = [[t_0 E_0], [t_1, E_1], ..., [t_{n-1}, E_{n-1}]]     
  
  to model
  
    E(t) = A tanh(a t  + b ) + B            for norm=false
    E(t) = 1 + c[tanh(at + b) − tanh(b)]    for norm=true
    
  Input:
    times:1D numpy array of floats
    E: 1D numpy array of floats
    norm: boolean
  
  Return:
    [a, b, c]         for norm=true
    [a, b, A, B]      for norm=false
'''

def fit_model(t, E, w = None, norm = False):
  
  x = np.array(t)
  y = np.array(E)
  
  scale = y[0]
  y /= scale
  
  n = len(x)
  
  if w is None:
    w = np.ones(n)
  
  def f(p):
    return np.sum(((y - model(x, p, norm))/w)**2)     
    
  def df(p):
    return np.dot(dmodel_dpars(x, p, norm), 2*(model(x, p, norm) - y)/w**2)
  
  min_kwargs = {"method": "BFGS", "jac": df, "tol": 1e-6}
  
  p0 = np.ones(4)
  
  T = f(p0)
  
  res = scipy.optimize.basinhopping(f, p0, T = T, minimizer_kwargs = min_kwargs, niter=200)
  
  p = res["x"]
  
  return {"err": (scale**2)*res["fun"], "p": np.concatenate((p[:2], scale*p[2:]))}
 
