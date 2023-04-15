from libraries import *
from constants import *
from potential import *


class Calculation:

    def __init__(self, n):

        self.n = n
        self.points = [self.Point.generate() for i in tqdm(range(n))]
        self.sigma = 0
        self.S = 0

    def sum_points(self):
        for i in range(self.n):
            self.S += self.points[i].fg * self.points[i].wg
        self.S = self.S / self.n
        for i in range(self.n):
            self.sigma += (self.points[i].fg * self.points[i].wg - self.S) ** 2
        self.S = self.S / (2 * np.pi * h) ** 5 * (1.408456 * 10 ** (-16)) ** 5
        self.sigma = (self.sigma / (self.n * (self.n - 1)))**0.5 / (2 * np.pi * h)**5 * (1.408456 * 10**(-16)) ** 5

    class Point:

        def __init__(self, Pr=0, cos_th=0, R=3.87, PR=0, wg=0):
            self.Pr = Pr
            self.cos_th = cos_th
            self.R = R
            self.PR = PR
            self.wg = wg
            self.fg = 0

        @classmethod
        def generate(cls):
            R = cls.sample_dist()
            if R >= 3.857452:
                cos_min, cos_max = [-1, 1]
            else:
                cos_max = cls.dihotomy_th(R)  # подбор пределов для th
                cos_min = -cos_max
            cos_th = cls.sample_angle(cos_min, cos_max)
            if V(R, cos_th) > 0:
                print(R, cos_th, cos_max)
            Pr_max = (2 * mu_O2 * (-V(R, cos_th))) ** 0.5

            Pr = cls.sample_dist(0, Pr_max)
            PR_max = (2 * mu_xeO2 * (-V(R, cos_th) - Pr**2 / (2 * mu_O2))) ** 0.5
            PR = cls.sample_dist(0, PR_max)
            wg = cls.Wg(Pr_max=Pr_max, PR_max=PR_max, cos_max=cos_max, cos_min=cos_min)
            return cls(Pr, cos_th, R, PR, wg)

        @staticmethod
        def sample_dist(dist_min=3.43523014, dist_max=8):
            return (dist_min ** 3 + np.random.random() * (dist_max ** 3 - dist_min ** 3)) ** (1 / 3)

        @staticmethod
        def sample_angle(cos_min, cos_max):
            return cos_min + np.random.random() * (cos_max - cos_min)

        @staticmethod
        def dihotomy_th(r, a=0, c=1, N=300):
            n = 0
            while n != N:
                n = n + 1
                while True:
                    try:
                        if V(r, a) * V(r, (a + c) / 2) <= 0:
                            c = (a + c) / 2
                        else:
                            a = (a + c) / 2
                        break
                    except (TypeError, RuntimeWarning):
                        a = a + 10 ** (-4)
            return (a + c) / 2

        @staticmethod
        def Wg(Pr_max= 0, R_min=3.435229845, R_max=8, PR_max=0, cos_max=1, cos_min=0):
            return 4 * np.pi * rO2 ** 2 * \
                   2 * np.pi / 3 * (R_max ** 3 - R_min ** 3) * (cos_max - cos_min) * \
                   np.pi * Pr_max ** 2 * \
                   4 / 3 * np.pi * PR_max ** 3

        def H(self):
            return self.Pr ** 2 / (2 * mu_O2) + self.PR ** 2 / (2 * mu_xeO2) + V(self.R, self.cos_th)

        def Fg(self):
            self.fg = np.exp(-self.H() / T)
