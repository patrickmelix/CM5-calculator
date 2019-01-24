#Calculate CM5 Charges
#
#
# Based on the work of Marenich et al.
# http://dx.doi.org/10.1021/ct200866d
# Please cite them when using this.
#

import numpy as np

#data from paper for element 1-118
_radii = np.array([0.32, 0.37, 1.30, 0.99, 0.84, 0.75,
          0.71, 0.64, 0.60, 0.62, 1.60, 1.40,
          1.24, 1.14, 1.09, 1.04, 1.00, 1.01,
          2.00, 1.74, 1.59, 1.48, 1.44, 1.30,
          1.29, 1.24, 1.18, 1.17, 1.22, 1.20,
          1.23, 1.20, 1.20, 1.18, 1.17, 1.16,
          2.15, 1.90, 1.76, 1.64, 1.56, 1.46,
          1.38, 1.36, 1.34, 1.30, 1.36, 1.40,
          1.42, 1.40, 1.40, 1.37, 1.36, 1.36,
          2.38, 2.06, 1.94, 1.84, 1.90, 1.88,
          1.86, 1.85, 1.83, 1.82, 1.81, 1.80,
          1.79, 1.77, 1.77, 1.78, 1.74, 1.64,
          1.58, 1.50, 1.41, 1.36, 1.32, 1.30,
          1.30, 1.32, 1.44, 1.45, 1.50, 1.42,
          1.48, 1.46, 2.42, 2.11, 2.01, 1.90,
          1.84, 1.83, 1.80, 1.80, 1.73, 1.68,
          1.68, 1.68, 1.65, 1.67, 1.73, 1.76,
          1.61, 1.57, 1.49, 1.43, 1.41, 1.34,
          1.29, 1.28, 1.21, 1.22, 1.36, 1.43,
          1.62, 1.75, 1.65, 1.57])


_Dz = np.array([0.0056, -0.1543, 0.0000, 0.0333, -0.1030, -0.0446,
      -0.1072, -0.0802, -0.0629, -0.1088, 0.0184, 0.0000,
      -0.0726, -0.0790, -0.0756, -0.0565, -0.0444, -0.0767,
       0.0130, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
       0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
      -0.0512, -0.0557, -0.0533, -0.0399, -0.0313, -0.0541,
       0.0092, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
       0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
      -0.0361, -0.0393, -0.0376, -0.0281, -0.0220, -0.0381,
       0.0065, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
       0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
       0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
       0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
       0.0000, 0.0000, -0.0255, -0.0277, -0.0265, -0.0198,
      -0.0155, -0.0269, 0.0046, 0.0000, 0.0000, 0.0000,
       0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
       0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
       0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
       0.0000, 0.0000, 0.0000, 0.0000, -0.0179, -0.0195,
      -0.0187, -0.0140, -0.0110, -0.0189])

_alpha = 2.474
_C = 0.705
_DHC = 0.0502
_DHN = 0.1747
_DHO = 0.1671
_DCN = 0.0556
_DCO = 0.0234
_DNO = -0.0346

def calc_cm5(atoms, hirschfeld):
    #all matrices have the naming scheme matrix[k,k'] according to the paper
    distances = atoms.get_all_distances(mic=True)
    atomicNumbers = np.array(atoms.numbers)
    Rz = _radii[atomicNumbers-1]
    RzSum = np.tile(Rz,(len(Rz),1))
    RzSum = np.add(RzSum, np.transpose(RzSum))
    Bkk = np.exp(-_alpha * (np.subtract(distances,RzSum)), out=np.zeros_like(distances), where=distances!=0)
    assert (np.diagonal(Bkk) == 0).all()

    Dz = _Dz[atomicNumbers]
#    Tkk = np.tile(Dz,(len(Dz),1))
#    Tkk = np.subtract(Tkk, np.transpose(Tkk))
    Tkk = np.zeros(shape=Bkk.shape)
    shape = Tkk.shape
    for i in range(shape[0]):
        for j in range(shape[1]):
            numbers = [atomicNumbers[i], atomicNumbers[j]]
            if numbers[0] == numbers[1]:
                continue
            if set(numbers) == set([1,6]):
                Tkk[i,j] = _DHC
                if numbers == [6,1]:
                    Tkk[i,j] *= -1.0
            elif set(numbers) == set([1,7]):
                Tkk[i,j] = _DHN
                if numbers == [7,1]:
                    Tkk[i,j] *= -1.0
            elif set(numbers) == set([1,8]):
                Tkk[i,j] = _DHO
                if numbers == [8,1]:
                    Tkk[i,j] *= -1.0
            elif set(numbers) == set([6,7]):
                Tkk[i,j] = _DCN
                if numbers == [7,6]:
                    Tkk[i,j] *= -1.0
            elif set(numbers) == set([6,8]):
                Tkk[i,j] = _DCO
                if numbers == [8,6]:
                    Tkk[i,j] *= -1.0
            elif set(numbers) == set([7,8]):
                Tkk[i,j] = _DNO
                if numbers == [8,7]:
                    Tkk[i,j] *= -1.0
            else:
                Tkk[i,j] = _Dz[numbers[0]-1] - _Dz[numbers[1]-1]
    assert (np.diagonal(Tkk) == 0).all()
    product = np.multiply(Tkk, Bkk)
    assert (np.diagonal(product) == 0).all()
    result = np.sum(product,axis=1)
    return np.array(hirschfeld) + result
