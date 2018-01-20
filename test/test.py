import numpy as np
import couplingmatrix as cp
import CP
import timeit

normalizedFreq = np.arange(-2.5, 2.5, 0.5)
M = np.array([[0.00000,0.44140,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000],
[0.44140,0.75180,0.16820,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000],
[0.00000,0.16820,0.78050,0.12880,0.00000,0.00000,0.00000,0.00000,0.00000],
[0.00000,0.00000,0.12880,0.78280,0.12080,-0.00420,-0.00840,0.00000,0.00000],
[0.00000,0.00000,0.00000,0.12080,0.79140,0.12890,0.00000,0.00000,0.00000],
[0.00000,0.00000,0.00000,-0.00420,0.12890,0.78290,0.12890,0.00000,0.00000],
[0.00000,0.00000,0.00000,-0.00840,0.00000,0.12890,0.78340,0.18000,0.00000],
[0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.18000,0.78280,0.46540],
[0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.46540,0.00000]])

S11_CP, S21_CP = CP.CM2S(M, normalizedFreq)
S11_cp, S21_cp = cp.CM2S(M, normalizedFreq)

np.testing.assert_array_equal(S11_CP, S11_cp)
np.testing.assert_array_equal(S21_CP, S21_cp)
np.allclose(S11_CP, S11_cp)
np.allclose(S21_CP, S21_cp)

setup_str = "import numpy as np; \
import CP; \
import couplingmatrix as cp; \
normalizedFreq = np.arange(-2.5, 2.5, 0.01); \
M = np.array([[0.00000,0.44140,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000], \
[0.44140,0.75180,0.16820,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000], \
[0.00000,0.16820,0.78050,0.12880,0.00000,0.00000,0.00000,0.00000,0.00000], \
[0.00000,0.00000,0.12880,0.78280,0.12080,-0.00420,-0.00840,0.00000,0.00000], \
[0.00000,0.00000,0.00000,0.12080,0.79140,0.12890,0.00000,0.00000,0.00000], \
[0.00000,0.00000,0.00000,-0.00420,0.12890,0.78290,0.12890,0.00000,0.00000], \
[0.00000,0.00000,0.00000,-0.00840,0.00000,0.12890,0.78340,0.18000,0.00000], \
[0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.18000,0.78280,0.46540], \
[0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.46540,0.00000]])"

num = 100
print(timeit.timeit('S11, S21 = CP.CM2S(M, normalizedFreq)', setup=setup_str, number=num))
print(timeit.timeit('S11, S21 = cp.CM2S(M, normalizedFreq)', setup=setup_str, number=num))