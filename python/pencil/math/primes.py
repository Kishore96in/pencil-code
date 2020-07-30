#!/usr/bin/env python
#coding: utf-8
# primes.py
# Written by Fred Gent (fred.gent.ncl@gmail.com)
"""
Find the prime factorization of the simulation grid, useful for checking 
permissible remesh parameters, parallelization options and fft compatibility.
"""
import numpy as np
from . import natural_sort
def prime_factors(n):
    i = 2
    factors = list()
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
            factors.append(i)
    if n > 1:
        factors.append(n)
    return factors

def divisors(factors, nmax, nmin=8):
    divisors = list()
    nn = len(factors)
    n = factors[0]
    divisors.append(n)
    for i in range(0,nn-1):
        n = factors[i-1] * n
        if n <= nmax/nmin and not n in divisors and np.mod(nmax,n)==0:
            divisors.append(n)
        for j in range(i,nn):
            m = factors[j]*factors[i]
            if m <= nmax/nmin and not m in divisors and np.mod(nmax,m)==0:
                divisors.append(m)
            for k in range(j+1,nn):
                p = factors[k]*factors[i-1]*factors[j]
                if p <= nmax/nmin and not p in divisors and np.mod(nmax,p)==0:
                    divisors.append(p)
    if 1 not in divisors:
        divisors.append(1)
    return natural_sort(divisors)

def common_factors(nx, ny, nz, nmin=8):
    factors = list()
    xfactors = prime_factors(nx)
    xfactors.append(1)
    yfactors = prime_factors(ny)
    yfactors.append(1)
    zfactors = prime_factors(nz)
    zfactors.append(1)
    xdivisors = divisors(xfactors,nx,nmin=nmin)
    ydivisors = divisors(yfactors,ny,nmin=nmin)
    zdivisors = divisors(zfactors,nz,nmin=nmin)
    print('divisors x:',xdivisors,'y:',ydivisors,'z:',zdivisors)
    return xdivisors, ydivisors, zdivisors

def cpu_optimal(nx,ny,nz,mvar=8,maux=0,par=dict(),nmin=32,MBmin=5.0,minghosts=7,
                quiet=True):
    xdiv, ydiv, zdiv = common_factors(nx,ny,nz,nmin=nmin)
    
    #fft_xyz_parallel removes need for x constraint
    #lpower = False
    #if len(par) > 0:
    #    if par['lpower_spectrum']:
    #        lpower = True
    #if lpower:
    #    xdiv = [1]
    nvar=mvar+maux
    #size of farray in MB
    nsize=8*nvar*nx*ny*nz/1024**2
    ncpu_max = int(nsize/MBmin)+1
    nprocyz_max = int(nsize/MBmin)+1
    print('ncpu_max',ncpu_max)
    cpu_list=list()
    for div in (xdiv,ydiv,zdiv):
        tmp = list()
        for nt in div:
            if nt <= ncpu_max:
                tmp.append(nt)
        cpu_list.append(tmp)
    xdiv, ydiv, zdiv = cpu_list
    ncpu_list=list()
    int_ncpu_max = 1
    ncpu_list.append(1)
    for div_list in (xdiv,ydiv,zdiv):
        for div in div_list:
            if not quiet:
                print(div, ncpu_list[-1])
            if div*int_ncpu_max < ncpu_max:
                int_ncpu_max = max(int_ncpu_max,div*int_ncpu_max)
            for cpu in ncpu_list:
                if div*cpu < ncpu_max:
                    int_ncpu_max = max(int_ncpu_max, div*cpu)
            ncpu_list.append(int_ncpu_max)
            if not quiet:
               print(div, ncpu_list[-1])
    ncpu_max = max(ncpu_list)
    #fft_xyz_parallel removes need for x constraint
    #if lpower:
    #    nprocyz_max = min(ncpu_max,max(ny,nz))
    #    ncpu_max = min(ncpu_max,max(ny,nz))
    print('ncpu_max divisible',ncpu_max)
    nprocx = max(xdiv)
    nprocy = max(ydiv)
    nprocz = 1
    nprocz_last = 1
    nprocy_last = 1
    for iprocx in xdiv:
        for iprocy in ydiv:
            for iprocz in zdiv:
                if np.mod(ncpu_max,iprocz*nprocy*nprocx) == 0:
                    nprocz = min(nz/iprocz,ncpu_max/nprocy/nprocx,nprocyz_max/nprocy)
                    if nprocz_last > nprocz:
                        nprocz = nprocz_last
                    nprocz_last = nprocz
                if not quiet:
                    print('nprocx {}, nprocy {}, nprocz {}, iprocz {}'.format(nprocx,nprocy,nprocz,iprocz))
            if np.mod(ncpu_max,nprocz*iprocy*nprocx) == 0:
                nprocy = min(ny/iprocy,ncpu_max/nprocz/nprocx)
                if nprocy_last > nprocy:
                    nprocy = nprocy_last
                nprocy_last = nprocy
            if not quiet:
                print('nprocx {}, nprocy {}, nprocz {}, iprocy {}'.format(nprocx,nprocy,nprocz,iprocy))
        if np.mod(ncpu_max,nprocz*nprocy*iprocx) == 0:
            nprocx = min(nx/iprocx,ncpu_max/nprocz/nprocy/iprocx,iprocx)
        if not quiet:
            print('nprocx {}, nprocy {}, nprocz {}, iprocx {}'.format(nprocx,nprocy,nprocz,iprocx))
    return [xdiv, ydiv, zdiv], [int(nprocx), int(nprocy), int(nprocz)]
