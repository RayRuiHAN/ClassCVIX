# -*- coding: UTF-8 -*-

import numpy as np
from ClassPrep import IOPrep


class CVIX(object):

    def __init__(self, input_file, input_path='.'):
        """
        rnd.info Date
        rnd.mu, rnd.sigma, rnd.skew, rnd.kurt
        """

        self.Data = IOPrep(input_path + '/' + input_file)
        self.info = self.Data.info
        self.Date = self.Data.TradingDate
        data1 = self.Data.data1
        data2 = self.Data.data2
        data1['order'] = range(len(data1))
        data2['order'] = range(len(data2))
        k1_list = data1.StrikePrice.unique()
        k2_list = data2.StrikePrice.unique()
        t1 = self.Data.t1
        t2 = self.Data.t2

        if self.Data.IfValid:

            def delta_k(x, i):
                k_list = i
                if x.StrikePrice == k_list[0]:
                    return k_list[1] - k_list[0]
                if x.StrikePrice == k_list[-1]:
                    return k_list[-1] - k_list[-2]
                else:
                    a = int(x.order + 1)
                    b = int(x.order - 1)
                    return 0.5 * (k_list[a] - k_list[b])

            def contrib(x, i):
                return x.delta_k / x.StrikePrice ** 2 * np.exp(self.Data.r * i) * x.ClosePrice_Merge

            data1['delta_k'] = data1.apply(lambda x: delta_k(x, k1_list), axis=1)
            data2['delta_k'] = data2.apply(lambda x: delta_k(x, k2_list), axis=1)
            data1['contrib'] = data1.apply(lambda x: contrib(x, self.Data.t1), axis=1)
            data2['contrib'] = data2.apply(lambda x: contrib(x, self.Data.t2), axis=1)

            var1 = (2 * data1.contrib.sum() - (self.Data.f1 / self.Data.k0_1 - 1) ** 2) / t1
            var2 = (2 * data2.contrib.sum() - (self.Data.f2 / self.Data.k0_2 - 1) ** 2) / t2
            vix0 = 100 * np.sqrt((0.5*var1/t1 + 0.5*var2/t2) / 12)
            m = 30
            w1 = 365 / m * t1 * (t2 - m / 365) / (t2 - t1)
            w2 = 365 / m * t2 * (m / 365 - t1) / (t2 - t1)
            self.cvix = w1 * var1 + w2 * var2
            #print(data1, var1, var2, vix, vix0/100**2)
            #print(self.Data.t1, self.Data.t2, 1/12)

            # w1 = (self.t2 - m) / (self.t2 - self.t1)
            # w2 = -(self.t1 - m) / (self.t2 - self.t1)
            # std_iv = np.sqrt((self.t1 * w1 * x.IV1 ** 2 + self.t2 * w2 * x.IV2 ** 2) * 12)
            # 0.0256/0.1644 + 0.0275/0.2411 0.106/0.0833

        else:
            nan = ''
            self.mu = nan


if __name__ == '__main__':
    v = CVIX("20200101.csv", input_path='./IO19-22')
    print(v.Date, v.cvix)
