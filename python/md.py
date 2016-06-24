"""
A simple molecular dynamics code to be used for nuclei ground
state calculation.
"""

import numpy as np
import itertools as it

N = 40
Z = 20
nsteps = 10000
m = 938
q0 = 6.00
p0 = 0.313
v0 = 34.32
hbarra = 1
#Here we use


# q and p because we won't work with canonical momentum equal to m*v

# f and g are the gradients in q and p respectively

def phi(ndq, ndp, t1, t2):
  """
  Calculate 1/q*dH/dq

  Parameters
  ----------

  ndq : float
      Distance between the particles

  ndp : float
      Relative momentum between the particles

  t1 : integer
      Type of particle 1

  t2 : integer
      Type of particle 2

  Returns
  -------

  phi : float
      1/q*dH/dq
  """
  if t1 == t2:
    base = - v0/q0**2 * (hbarra/(q0 * p0))**3
    exponente = -ndp * ndp / (2 * p0**2) - ndq * ndq / (2 * q0**2)
    return base * np.exp(exponente)
  else:
    return 0

def psi(ndq, ndp, t1, t2):
  """
  Calculate 1/p*dH/dp

  Parameters
  ----------

  ndq : float
      Distance between the particles

  ndp : float
      Relative momentum between the particles

  t1 : integer
      Type of particle 1

  t2 : integer
      Type of particle 2

  Returns
  -------

  psi : float
      1/pd*dH/dp
  """

  if t1 == t2:
    base = - v0/p0**2 * (hbarra/(q0 * p0))**3
    exponente = -ndp * ndp / (2 * p0**2) - ndq * ndq / (2 * q0**2)
    return base * np.exp(exponente)
  else:
    return 0


def first_step(q, p, f, g):
  """
  First step of a symplectic integration with a potential that
  depends on both p and q. Modifies q and p in place.

  Parameters
  ----------

  q : numpy array
      Positions of the particles

  p : numpy array
      Momenta of the particles

  f : numpy array
      Forces (in q) exerted on the particles

  g : numpy array
      Forces (in p) [Gorces?] exerted on the particles
  """
  pass

def forces(q, p, t):
  """
  Calculate forces and gorces (?) for a q and p system

  Parameters
  ----------

  q : numpy array
      Positions of the particles

  p : numpy array
      Momenta of the particles

  t : numpy array
      Types of the particles

  Returns
  -------

  f, g : numpy arrays
  """

  natoms = np.shape(q)[0]
  f = np.zeros_like(q)
  g = np.zeros_like(q)
  for i in xrange(natoms):
    for j in xrange(i+1, natoms):
      dq = q[i, :] - q[j, :]
      ndq = np.linalg.norm(dq)
      dp = p[i, :] - p[j, :]
      ndp = np.linalg.norm(dp)
      t1 = t[i]
      t2 = t[j]
      f[i, :] += phi(ndq, ndp, t1, t2) * dq
      g[i, :] += psi(ndq, ndp, t1, t2) * dp
      f[j, :] -= phi(ndq, ndp, t1, t2) * dq
      g[j, :] -= psi(ndq, ndp, t1, t2) * dp
    g[i, :] += p[i, :] / m
  return f, g


def final_step(q, p, f, g):
  """
  Final step of a symplectic integration with a potential that
  depends on both p and q. Modifies q and p in place.

  Parameters
  ----------

  q : numpy array
      Positions of the particles

  p : numpy array
      Momenta of the particles

  f : numpy array
      Forces (in q) exerted on the particles

  g : numpy array
      Forces (in p) [Gorces?] exerted on the particles
  """
  h = 0.05
  p -= h*f
  q += h*g

def init_particles(N, Z):
  """

  """
  L = int(N**(1.0/3.0)) + 1
  q = np.zeros((N, 3))
  p = np.zeros((N, 3))
  t = np.zeros(N)
  ipart = 0
  d = 6.0

  for i, j, k in it.product(range(L), range(L), range(L)):
    q[ipart, :] = d*np.array([i, j, k])
    if ipart < Z:
      t[ipart] = 2
    else:
      t[ipart] = 1
    ipart += 1
    if ipart == N:
      break
  return q, p, t

q, p, t = init_particles(N, Z)
f = np.zeros_like(q)
g = np.zeros_like(q)
for i in xrange(nsteps):
  first_step(q, p, f, g)
  f, g = forces(q, p, t)
  final_step(q, p, f, g)
  print i
