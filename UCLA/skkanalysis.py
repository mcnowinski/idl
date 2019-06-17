# -*- coding: utf-8 -*-
"""
Created on Thu Aug 28 17:23:37 2014

@author: Szilard Gyalay
-.-.-

Description:
This file consists of a set of functions that are called by the last function (skkanalysis).
These functions are:
idl_tabulate
reverseshkuratov
skk

Requires numpy and scipy to run.
-.-.-

To use skkanalysis, run this file. Then can use the syntax
skkanalysis(wavelengths, albedo, refwave, refreal, q, S, filename)

Results will then be printed to the file, with iteration #'s printed.
-.-.-

skkanalysis inputs:

Wavelengths: The wavelengths corresponding to the albedos. Because the cauchy integral
is from 0 to inf with respect to ln(lambda), and ln(1)=0, the units used for the
wavelengths should be smaller than the smallest wavelength. E.g. if the smallest
wavelength is .4 microns, use nm.

albedo: The albedos corresponding to the wavelengths.

refwave: The reference wavelength for which the real index is known. Typically
within the wavelength range of the known albedos. In same units as wavelengths.
refreal: The real index of refraction at refwave

q: The volume fraction filled by particles (1-porosity)

S: Particle size (in same units as wavelengths)

filename: a string with the name of the file where the results are printed.
Typically a .txt file.
"""


# Need to import some things.
# SciPy is a module full of scientific functions
import scipy
# NumPy expands the mathematical functions available to python.
import numpy as np
import glob
import re

import logging
# logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("skkanalysis.log"),
        logging.StreamHandler()
    ])
logger = logging.getLogger('skkanalysis')

# This allows us to integrate from a set of x and f(x) values


def idl_tabulate(x, f, p=5):
    """
    Let's us integrate like IDL's int_tabulated
    This program was found online.
    """
    def newton_cotes(x, f):
        if x.shape[0] < 2:
            return 0
        rn = (x.shape[0] - 1) * (x - x[0]) / (x[-1] - x[0])
        weights = scipy.integrate.newton_cotes(rn)[0]
        return (x[-1] - x[0]) / (x.shape[0] - 1) * np.dot(weights, f)
    ret = 0
    for idx in xrange(0, x.shape[0], p - 1):
        ret += newton_cotes(x[idx:idx + p], f[idx:idx + p])
    return ret


# To use the number pi
from math import pi
# To use log functions.
from math import log

# This allows us to get k from n (and other properties)


def reverseshkuratov(albedo, n, S, q, wavelengths):
    """
    Using reflectance data, we can determine complex index of refraction.
    Uses SHkuratov et al. 1999 (S99) eqs 13

    Albedo is the hemispheric albedo
    n is the real index of refraction
    S is the mean path length of light, basically particle size
    q is the volume fraction filled by particles--basically 1-porosity
    Wavelengths is a list of the corresponding wavelengths
    """

    # K: basically initializing a bunch of arrays
    # Here we pre-populate a shitload of lists.
    # r_0 the Fresnel coefficient at normal incidence
    r_0 = [0.0 for i in range(len(wavelengths))]
    # R_e and R_i the average external and internal reflection coefficients
    # Use S99 eqs 8 for n=1.4-1.7 (accurately, can use for broader range)
    R_e = [0.0 for i in range(len(wavelengths))]
    R_i = [0.0 for i in range(len(wavelengths))]
    # Average backward reflectance coefficient R_b S99 eq 8a
    R_b = [0.0 for i in range(len(wavelengths))]
    # Average transmittances using S99 eq 7
    T_e = [0.0 for i in range(len(wavelengths))]
    T_i = [0.0 for i in range(len(wavelengths))]
    # Definig variables from S99 eqs 13
    y = [0.0 for i in range(len(wavelengths))]
    a = [0.0 for i in range(len(wavelengths))]
    b = [0.0 for i in range(len(wavelengths))]
    c = [0.0 for i in range(len(wavelengths))]
    # The actual imaginary index of refraction S99 eq 13
    k = [0.0 for i in range(len(wavelengths))]
    # Python is not matrix/array oriented, so need a for loop to do the more complicated math on all elements of a list
    for i in range(len(wavelengths)):
        # Uses Shkuratov et al. 1999 equations listed above
        r_0[i] = (((n[i]-1.0)/(n[i]+1.0))**2)
        R_e[i] = (r_0[i]+0.05)
        R_i[i] = (1.04-1.0/(n[i]**2))
        R_b[i] = ((0.28*n[i]-0.2)*R_e[i])
        T_e[i] = (1.0-R_e[i])
        T_i[i] = (1.0-R_i[i])
        y[i] = (((1.0-albedo[i])**2)/(2.0*albedo[i]))
        a[i] = (T_e[i]*T_i[i]*(y[i]*R_i[i]+q*T_e[i]))
        b[i] = (y[i]*R_b[i]*R_i[i]+0.5*q*(T_e[i]**2)*(1.0+T_i[i]) -
                T_e[i]*(1.0-q*R_b[i]))
        c[i] = (2.0*y[i]*R_b[i]-2.0*T_e[i]*(1.0-q*R_b[i])+q*T_e[i]**2)
        # Python will throw errors if the argument of the square root or log is negative, so we make these NaN's.
        if ((b[i]/a[i])**2-(c[i]/a[i])) < 0:
            k[i] = float('nan')
        elif ((b[i]/a[i])+((b[i]/a[i])**2-(c[i]/a[i]))**0.5) < 0:
            k[i] = float('nan')
        else:
            k[i] = (-wavelengths[i]*log((b[i]/a[i]) +
                                        ((b[i]/a[i])**2-(c[i]/a[i]))**0.5)/(4.0*pi*S))
    return(k)


# so we can use euler's number
from math import e
# To test whether a number is not a number ot not
from math import isnan
# To use fanyc-schmancy integration.
from scipy import integrate

# This allows us to get n from k


def skk(wavesprime, imsprime, refwave, refreal):
    """
    For now, assuming we actually have values of k for every wavelength
    Can find n for every wavelength
    Uses formula n(l_i)=n_v+(2/pi)* P int_0^inf(l**2*k(l)/
        (l_v**2-l**2)*(l_i**2-l**2))d(ln(l))
    from Dalton+Pitman 2012, Roush 2010
    """
    # declaring some lists
    wavelengths = []
    imsafter = []
    for i in range(len(wavesprime)):
        # Checking if list elements are finite
        if isnan(imsprime[i]) == False:
            # Populating the lists with the finite results
            wavelengths.append(wavesprime[i])
            imsafter.append(imsprime[i])

    # Need to get rid of nan values from the imaginary index
    # Setting up an array to fill later
    realsafter = []
    # Boundary conditions for our range of k values
    k_a = imsafter[0]
    lambda_a = wavelengths[0]
    k_b = imsafter[len(imsafter)-1]
    lambda_b = wavelengths[len(wavelengths)-1]
    # Need to get logs ofwavelengths since integral is taken with respect to it
    lnwavelengths = []
    for i in range(0, len(wavelengths)):
        lnwavelengths.append(log(wavelengths[i]))
    # Now we can actually start the damn thing
    # Need to compute the integral for every wavelength we want n for.
    for i in range(0, len(wavelengths)):
        # This part is what the integral is multiplied by
        multi = (2/pi)*(refreal**2-wavelengths[i]**2)
        # When we get to the refwave, it is equal to refreal, so no need for integral
        if multi == 0:
            cauchy = 0
        else:
            # We can essentially split the integral into 7 different parts
            # Need to set up singularities: The wavelength we're calculating n for
            # and the reference wavelength
            if wavelengths[i] < refreal:
                # if current wavelength lest than ref, first sing is current
                lambda1 = wavelengths[i]
                lambda2 = refreal
            elif wavelengths[i] > refreal:
                # the other way
                lambda1 = refreal
                lambda2 = wavelengths[i]
            # Part 1: ln(l) from 0 to ln(lambda_a)
            cauchy1 = 0
            # Creates a function to integrate over.

            def func1(x): return ((e**x)**2)*k_a / \
                ((refreal**2-(e**x)**2)*(wavelengths[i]**2-(e**x)**2))
            # Does the integration.
            inte1 = integrate.quad(func1, log(1), log(lambda_a))
            # Integration gives a few results, the first is the answer.
            cauchy1 = inte1[0]
            # Part 2: from ln(lambda_a) to ln(lambda1)
            cauchy2 = 0
            # Going to do a tabluated integral like in IDL
            # Declare some vatiables for our range
            waveslt1 = []
            lnwaveslt1 = []
            imslt1 = []
            for j in range(0, len(wavelengths)):
                if wavelengths[j] < lambda1:
                    # Get the wavelengths from beginning of our range up to first singularity
                    waveslt1.append(wavelengths[j])
                    # get natural log of them for using the function
                    lnwaveslt1.append(lnwavelengths[j])
                    # Has the imaginary indices for that range
                    imslt1.append(imsafter[j])
            # Make the function for our range
            funclt1 = []
            for j in range(0, len(waveslt1)):
                funclt1.append(
                    (waveslt1[j]**2)*imslt1[j]/((lambda1**2-waveslt1[j]**2)*(lambda2**2-waveslt1[j]**2)))
            # Making a numpy array to use the idl_tabulate function
            lnwaveslt1np = np.array(lnwaveslt1)
            funclt1np = np.array(funclt1)
            # doing the actual integral
            cauchy2 = idl_tabulate(lnwaveslt1np, funclt1np)
            # Part 4: Between the singularities
            cauchy4 = 0
            # Tabulated integral again, like in Part 2
            # Gets wavelengths and function of those wavelengths between the singularities
            wavesbtw = []
            lnwavesbtw = []
            imsbtw = []
            for j in range(0, len(wavelengths)):
                if wavelengths[j] > lambda1:
                    if wavelengths[j] < lambda2:
                        wavesbtw.append(wavelengths[j])
                        lnwavesbtw.append(lnwavelengths[j])
                        imsbtw.append(imsafter[j])
            funcbtw = []
            for j in range(0, len(wavesbtw)):
                funcbtw.append(
                    (wavesbtw[j]**2)*imsbtw[j]/((lambda1**2-wavesbtw[j]**2)*(lambda2**2-wavesbtw[j]**2)))
            lnwavesbtwnp = np.array(lnwavesbtw)
            funcbtwnp = np.array(funcbtw)
            cauchy4 = idl_tabulate(lnwavesbtwnp, funcbtwnp)
            # Part 3: Around the singularity at ln(lambda1)
            cauchy3 = 0
            if len(waveslt1) > 0:
                # If there exists a wavelength before the singularity that we have
                # Then that wavelength is where the integral is calculated from
                lambda1_1 = waveslt1[len(waveslt1)-1]
                k1_1 = imslt1[len(imslt1)-1]
                if len(wavesbtw) > 0:
                    # If there exists wavelengths after the singularity
                    lambda1_2 = wavesbtw[0]
                    k1_2 = imsbtw[0]
                    # Makes the imaginary value between between the imaginary values corresponding to those wavelengths
                    ksing1 = k1_1+(k1_2-k1_1) * \
                        (lambda1-lambda1_1)/(lambda1_2-lambda1_1)
                elif len(wavesbtw) == 0:
                    # If there are no wavelengths between our singularities, then other limit of integration between them
                    lambda1_2 = lambda1+0.5*(lambda2-lambda1)
                    k1_2 = imsafter[i]
                    ksing1 = k1_1+(k1_2-k1_1) * \
                        (lambda1-lambda1_1)/(lambda1_2-lambda1_1)
            elif len(waveslt1) == 0:
                # If there are no wavelengths right before the first singularity,
                # means that the first wavelength we have is our singularity
                # So we integrate from slightly before the singularity
                lambda1_1 = wavelengths[0]-0.001
                # im index same as that for first wavelength
                k1_1 = imsafter[0]
                # upper limit same as before
                if len(wavesbtw) > 0:
                    lambda1_2 = wavesbtw[0]
                    k1_2 = imsbtw[0]
                    ksing1 = k1_1
                elif len(wavesbtw) == 0:
                    lambda1_2 = lambda1+0.5*(lambda2-lambda1)
                    k1_2 = imsafter[i]
                    ksing1 = k1_1
            # Now making the function to integrate over

            def func3(x): return ((e**x)**2)*ksing1 / \
                ((lambda1**2-(e**x)**2)*(lambda2**2-(e**x)**2))
            # Declaring the bounds we figued out earlier
            bound1_1 = log(lambda1_1)
            bound1_2 = log(lambda1_2)
            # Doing the actual integration
            inte3 = integrate.quad(
                func3, bound1_1, bound1_2, weight='cauchy', wvar=log(lambda1))
            # Getting the answer from integrating
            cauchy3 = inte3[0]
            # Part 6: From ln(lambda2) to ln(lambda_b)
            cauchy6 = 0
            # Tabulated integral again, like in Part 2, 4
            # This time stuff after the second singularity
            wavesgt2 = []
            lnwavesgt2 = []
            imsgt2 = []
            for j in range(0, len(wavelengths)):
                if wavelengths[j] > lambda2:
                    wavesgt2.append(wavelengths[j])
                    lnwavesgt2.append(lnwavelengths[j])
                    imsgt2.append(imsafter[j])
            funcgt2 = []
            for j in range(0, len(wavesgt2)):
                funcgt2.append(
                    (wavesgt2[j]**2)*imsgt2[j]/((lambda1**2-wavesgt2[j]**2)*(lambda2**2-wavesgt2[j]**2)))
            lnwavesgt2np = np.array(lnwavesgt2)
            funcgt2np = np.array(funcgt2)
            cauchy6 = idl_tabulate(lnwavesgt2np, funcgt2np)
            # Part 5: Around singularity at ln(lambda2)
            cauchy5 = 0
            if len(wavesbtw) > 0:
                # If stuff between singularity, grab last wavelength before
                lambda2_1 = wavesbtw[len(wavesbtw)-1]
                k2_1 = imsbtw[len(imsbtw)-1]
                if len(wavesgt2) > 0:
                    # If stuff after singularity, grab first wavelength after
                    lambda2_2 = wavesgt2[0]
                    k2_2 = imsgt2[0]
                    ksing2 = k2_1+(k2_2-k2_1) * \
                        (lambda2-lambda2_1)/(lambda2_2-lambda2_1)
                    # Cauchy weighting makes values ridiculously negative
                elif len(wavesgt2) == 0:
                    # If nothing after second singularity, grab last wavelength
                    # But wait, that's a bit silly. That would make the last wavelength the singularity
                    # No wonder integration doesn't work here. Have to fix.
                    lambda2_2 = wavelengths[len(wavelengths)-1]
                    k2_2 = imsafter[len(imsafter)-1]
                    ksing2 = k2_2
                    # Cauchy weighting does not want to work here
            elif len(wavesbtw) == 0:
                # If no wavelengths between, make one between to integrate from
                lambda2_1 = lambda1+0.5*(lambda2-lambda1)
                k2_1 = imsafter[i]
                if len(wavesgt2) > 0:
                    # If wavelengths after, grab the first one after
                    lambda2_2 = wavesgt2[0]
                    k2_2 = imsgt2[0]
                    ksing2 = k2_1+(k2_2-k2_1) * \
                        (lambda2-lambda2_1)/(lambda2_2-lambda2_1)
                    # Cauchy weighting here makes the first value ridiculously negative
                elif len(wavesgt2) == 0:
                    # If no wavelengths after, use the last wavelength
                    lambda2_2 = wavelengths[len(wavelengths)-1]
                    k2_2 = imsafter[len(imsafter)-1]
                    ksing2 = k2_2
                    # Same values for all whether this is weighted cauchy or just breakpoints labeled.

            def func5(x): return ((e**x)**2)*ksing2 / \
                ((lambda1**2-(e**x)**2)*(lambda2**2-(e**x)**2))
            # Bounds for the equation
            bound2_1 = log(lambda2_1)
            bound2_2 = log(lambda2_2)
            # Integration
            # inte5=integrate.quad(func5,bound2_1,bound2_2,weight='cauchy',wvar=log(lambda2))
            inte5 = integrate.quad(
                func5, bound2_1, bound2_2, points=[log(lambda2)])
            cauchy5 = inte5[0]
            # Part 7: From ln(lambda_b) to infinity
            cauchy7 = 0

            def func7(x): return k_b/((1/((e**x)**2)) *
                                      (lambda1**2-(e**x)**2)*(lambda2**2-(e**x)**2))
            inte7 = integrate.quad(func7, log(lambda_b), log(1e12))
            cauchy7 = inte7[0]
            # Adding it all together
            cauchy = (cauchy1+cauchy2+cauchy3+cauchy4+cauchy5+cauchy6+cauchy7)
        # Getting final results
        result = refreal+multi*cauchy
        # Adding the result to the list.
        realsafter.append(result)
    # Getting the results plus NaNs.
    actualrealsafter = []
    anumber = 0
    for i in range(len(wavesprime)):
        if isnan(imsprime[i]):
            actualrealsafter.append(float('nan'))
            anumber = anumber+1
        else:
            actualrealsafter.append(realsafter[i-anumber])
    return(actualrealsafter)

# Now for the actual iterator


def skkanalysis(wavelengths, albedo, refwave, refreal, q, S, filename):
    """
    Here we do some hard-core iterating.
    Recommended that you use nm. Integrals go from 0-inf of ln(wave)
    """
    itr = 0
    # Will iterate 50 times if not broken
    while itr < 50:
        print('iteration '+str(itr+1))
        # Initializing
        if itr == 0:
            realsbefore = []
            for i in range(0, len(wavelengths)):
                realsbefore.append(refreal)
        elif itr > 0:
            realsbefore = realsafter
        # Ensuring this won't loop forever
        itr = itr+1
        # Getting k from n
        imsafter = reverseshkuratov(albedo, realsbefore, S, q, wavelengths)
        # Getting n from k
        realsafter = skk(wavelengths, imsafter, refwave, refreal)
        # Testing to see if n has converged
        # i.e. n differs <1% between iterations
        testing = []
        for i in range(0, len(realsafter)):
            if ((realsbefore[i]-realsafter[i])/(realsbefore[i])) > 0.01:
                testing.append(1)
        if len(testing) == 0:
            # Stops the loop if it has converged
            break
    imsafter = reverseshkuratov(albedo, realsbefore, S, q, wavelengths)
    results = []
    for i in range(0, len(wavelengths)):
        # Puts results into an array
        results.append([wavelengths[i], realsafter[i], imsafter[i]])
    # Writes the results to file
    file = open(filename, 'w')
    file.write(
        'Optical Constants derived by Shkuratov et al. Theory and SK-K Analysis\n')
    file.write('Wavelength    real index    im index\n')
    for i in range(len(results)):
        file.write(str(results[i][0])+'    ' +
                   str(results[i][1])+'    '+str(results[i][2])+'\n')
    file.close()
    print('file '+filename+' printed after '+str(itr)+' iterations.\n')


def main():
    logger.info('Starting skkanalysis...')

    input_path = ''

    skk_files = glob.glob(input_path+'*.skk')
    print 'Found %d input (.skk) files.' % len(skk_files)

    for skk_file in skk_files:
        wavelengths = []  # wavelengths for reflectance data
        albedos = []  # reflectance valued for specific wavelength
        n = None  # real part of index of refraction
        wavelength = None  # wavlength where r is valid; in same units as wavelengths
        q = None  # porosity (0-1)
        S = None  # particle size in same units as wavelengths
        try:
            logger.info('Opening file (%s)...' % skk_file)
            f = open('%s' % skk_file)
        except Exception as e:
            logger.error('Failed to open file (%s).' % skk_file)
            continue
        # get line containing r, wavelength, q, and S
        for line in f:
            if re.search('^#', line.strip()):
                continue
            else:
                try:
                    n, wavelength, q, S = map(float, line.split(','))
                    logger.info('Index of refraction (real) is %.3f.' % n)
                    logger.info('Reference wavelength is %.3f.' % wavelength)
                    logger.info('Porosity is %.3f.' % q)
                    logger.info('Particle size is %.3f.' % S)
                except Exception as e:
                    logger.error('Invalid file format (%s).' % skk_file)
                    raise
                break
        # get wavelength, albedo array
        for line in f:
            if re.search('^#', line.strip()):
                continue
            (x, y) = map(float, line.split(','))
            wavelengths.append(x)
            albedos.append(y)
        # logger.debug(wavelengths)
        # logger.debug(albedos)
        logger.info('Read %d albedo vs. wavelength data points.' %
                    len(wavelengths))
        f.close()

        output_path = input_path + skk_file + \
            '.n%.3f@%.3f.q%.3f.S%0.3f.out' % (n, wavelength, q, S)
        logger.debug(output_path)
        # perform the skk calculations
        skkanalysis(wavelengths, albedos, wavelength, n, q, S, output_path)


if __name__ == "__main__":
    main()
