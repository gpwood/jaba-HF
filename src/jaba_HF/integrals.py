import numpy as np
import math
import rdkit
from rdkit import Chem
import re
from pathlib import Path
from dataclasses import dataclass
import decimal

@dataclass
class pt_table:
    """Class  pt table  """

    symbols  = [
     'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl',
     'Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As',
     'Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',
     'Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb',
     'Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl',
     'Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk',
     'Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs']

@dataclass
class gauss_contraction:
    """ class to store a gaussian contraction """
    shell: str
    exp: list[float]
    coeff: list[float]
    scale: float

def read_gaussian_basis(file: Path) -> dict:
    """ read a basis set file of gaussians in gaussian format and return as dictionary of 
    atoms with coefficients and exponents"""

    basis_set: dict = {}
    basis_element = ""
    gauss_types: list = []
    gauss_numbers: list = []
    gauss_scales: list = []

    exponents: list = []
    coefficents: list = []

    gaussians: dict = {}

    with open(file, 'r') as f:
        data = f.readlines()

    for line in data:

        line = line.lstrip().strip()
        element_re = re.search("^([A-Z][a-z]?)\s+0$",line)
        gauss1_re = re.search("^(\d+.\d+D[+,-]\d+)\s+([-]?\d+.\d+D[+,-]\d+)$",line)
        gauss2_re = re.search("^(\d+.\d+D[+,-]\d+)\s+([-]?\d+.\d+D[+,-]\d+)\s+([-]?\d+.\d+D[+,-]\d+)$",line)
        gausstype_re = re.search("(S|SP|D|F)\s+(\d+)\s+(\d+.\d+)$",line)

        if re.search('^!',line):
            # this is a comment
            pass
        elif line == "":
            # this is a blank line
            pass
        elif element_re:
            pt = pt_table()
            if element_re.group(1) in pt.symbols:
                basis_element = element_re.group(1)
            else:
                raise Exception(f"Unrecognized element in basis set file {element_re.group(1)}")
        elif gausstype_re:
            gauss_types.append(gausstype_re.group(1))
            gauss_numbers.append(int(gausstype_re.group(2)))
            gauss_scales.append(float(gausstype_re.group(3)))
        elif gauss1_re:
            exponents.append(decimal.Decimal(gauss1_re.group(1).replace('D','E')))
            coefficents.append(decimal.Decimal(gauss1_re.group(2).replace('D','E')))
        elif gauss2_re:
            exponents.append(decimal.Decimal(gauss2_re.group(1).replace('D','E')))
            exponents.append(decimal.Decimal(gauss2_re.group(1).replace('D','E')))
            coefficents.append(decimal.Decimal(gauss2_re.group(2).replace('D','E')))
            coefficents.append(decimal.Decimal(gauss2_re.group(3).replace('D','E')))
        elif line == "****":
            gaussians: list = []
            counter = 0
            for i,number in enumerate(gauss_numbers):
                gaussian = {}
                gaussian2 = {}
                exps = []
                coeffs = []
                exps2 = []
                coeffs2 = []
                gauss_type = gauss_types[i]
                gauss_scale = gauss_scales[i]
                if gauss_type == 'SP':
                    for j in range(number):
                        exps.append(exponents[counter])
                        coeffs.append(coefficents[counter])
                        counter += 1
                        exps2.append(exponents[counter])
                        coeffs2.append(coefficents[counter])
                        counter += 1
                    ctg=gauss_contraction('S',exps,coeffs,gauss_scale)
                    gaussians.append(ctg)
                    ctg=gauss_contraction('P',exps2,coeffs2,gauss_scale)
                    gaussians.append(ctg)
                else:
                    for j in range(number):
                        exps.append(exponents[counter])
                        coeffs.append(coefficents[counter])
                        counter += 1
                    ctg=gauss_contraction(gauss_type,exps,coeffs,gauss_scale)
                    gaussians.append(ctg)

            basis_set[basis_element] = gaussians
            gauss_types = []
            gauss_numbers = []
            gauss_scales = []
            exponents = []
            coefficents = []
            gaussians = []
        else:
            pass

    return basis_set

def overlap(ctg1,ctg2,r):
    pass


#def integrals(iop,N,R,zeta1,zeta2,Za,Zb):
#
#    ### N is STO-NG:
#H     0
#S    3   1.00
#      0.3425250914D+01       0.1543289673D+00
#      0.6239137298D+00       0.5353281423D+00
##      0.1688554040D+00       0.4446345422D+00
#****
#He     0
#S    3   1.00
#      0.6362421394D+01       0.1543289673D+00
#      0.1158922999D+01       0.5353281423D+00
#      0.3136497915D+00       0.4446345422D+00
#****
#
#
#    ### the integrals to calculate
#    ## numpy.longdouble?
#    s12,t11,t12,t22,v11a,v12a,v22a,v11b,v12b,v22b,v1111,v2111,v2121,v2211,v2221,v2222
    ###  basis set
#    coeffcients = np.zeros((3,3)) ## the basis set coeffs read in
##    exponents = np.zeros((3,3)) ## the exponents read in
#    d1 = np.zeros(3) ## normalize the coefficients
#    a1 = np.zeros(3) ## scale the exponents by the slater coefficient
#    d2 = np.zeros(3)
#    a2 = np.zeros(3)
#
#    # CGF-1S(zeta=1.0,STO-1G) = 1.0*GF(0.270905)
####
