import numpy as np
import matplotlib.pyplot as plt
import multiprocessing

import unidec
from unidec_modules import unidectools as ud

__author__ = 'Michael.Marty'


# Test for error function:
# 1. Peak Width
# 2. Noise
# 3. Number of Charge States
# 4. Baseline
# 5. Ambiguity
# Invariant towards:
# 1. Data density
# 2. Peak Shape
# 3. Different peak widths
# 4. Different Intensities
# 5. MZsig

class UniDecTest(unidec.UniDec):
    def test_spectra(self, mlist=None, ilist=None, res=1000, noise=0.0, pad=2000, window=None,
                     plot=False, restricted=True, massbins=10, psfun=0, mzsig=None, **kwargs):
        if mlist is None:
            mlist = [200000, 210000]
        if ilist is None:
            ilist = [100, 100]

        m1 = mlist[0]
        m2 = mlist[1]
        if mzsig is None:
            mzsig = m1 / float(res) / float(ud.predict_charge(m1)) / 2.

        self.open_test_spectrum(masslist=mlist, intlist=ilist, noise=noise, resolution=res, psfun=psfun, silent=True,
                                **kwargs)
        self.config.massbins = massbins
        self.config.mzsig = mzsig
        self.config.psfun = psfun
        self.config.linflag = 2
        self.config.mzbins = 0
        self.config.peaknorm = 2
        if window is None:
            self.config.peakwindow = pad
        else:
            self.config.peakwindow = window
        if restricted:
            self.config.masslb = m1 - pad
            self.config.massub = m2 + pad
            self.config.endz = np.amax(self.data.ztab)
            self.config.startz = np.amin(self.data.ztab)
        self.process_data(silent=True)
        self.run_unidec(silent=True)
        self.config.peakthresh = 0.1
        self.pick_peaks()
        self.autointegrate()
        self.get_errors(**kwargs)

        # check
        results = np.empty((len(mlist), 4))
        vals = np.array(
            [self.errorgrid[:, 0, 1], self.errorgrid[:, 1, 1], self.errorgrid[:, 0, 3], self.errorgrid[:, 1, 3]])
        peaks = []
        for i in xrange(0, len(mlist)):
            truemass = mlist[i]
            trueint = ilist[i]

            for j, e in enumerate(self.errorgrid):
                fitmass = self.errorgrid[j, 0, 0]
                corrmass = self.errorgrid[j, 1, 0]
                if np.abs(truemass - corrmass) > self.config.peakwindow:
                    pass
                else:
                    fitmasserr = self.errorgrid[j, 0, 1]
                    # TODO: Aysymmetric CIs
                    if np.abs(fitmass - truemass) < fitmasserr:
                        results[i, 0] = True
                    else:
                        results[i, 0] = False

                    corrmasserr = self.errorgrid[j, 1, 1]
                    if np.abs(corrmass - truemass) < corrmasserr:
                        results[i, 1] = True
                    else:
                        results[i, 1] = False

                    fitarr = self.errorgrid[j, 0, 2]
                    fitarrerr = self.errorgrid[j, 0, 3]
                    # TODO: Aysymmetric CIs
                    if np.abs(fitarr - trueint) < fitarrerr:
                        results[i, 2] = True
                    else:
                        results[i, 2] = False

                    corrint = self.errorgrid[j, 1, 2]
                    correrr = self.errorgrid[j, 1, 3]

                    if np.abs(corrint - trueint) < correrr:
                        results[i, 3] = True
                    else:
                        results[i, 3] = False
                    peaks.append(self.pks.peaks[j])
        totres = np.all(results, axis=0)

        if False:
            print "Fit Mass:", totres[0]
            print "Peak Mass:", totres[1]
            print "Fit Area:", totres[2]
            print "Corr Area:", totres[3]
            print vals
        if plot:
            self.make_plot(massrange=[m1 - pad, m2 + pad])
            # self.config.print_config()
        return totres, vals, peaks

        pass

    def simpleplotter(self, m1, m2, e1, e2, i):
        plt.subplot(2, 4, i + 1)
        plt.plot(m1)
        plt.plot(m2)
        plt.subplot(2, 4, i + 5)
        plt.plot(e1)
        plt.plot(e2)
        n = len(m1)
        tvalue = ud.get_tvalue(n - 1)
        zvalue = ud.get_zvalue()
        #print tvalue,zvalue
        #plt.hlines(np.std(m1,ddof=1) * zvalue, 0, n, color="r")
        #plt.hlines(np.std(m2,ddof=1) * zvalue, 0, n, color="r")
        plt.hlines(np.std(m1,ddof=1)*tvalue, 0, n, color="m")
        plt.hlines(np.std(m2,ddof=1)*tvalue, 0, n, color="m")
        plt.hlines(np.mean(e1), 0, n, color="orange")
        plt.hlines(np.mean(e2), 0, n, color="orange")

    def plot_repeats(self, peaks):

        plt.figure(figsize=(12, 12))

        m1 = [p.fitmassavg for p in np.array(peaks)[:, 0]]
        m2 = [p.fitmassavg for p in np.array(peaks)[:, 1]]
        e1 = [p.fitmasserr for p in np.array(peaks)[:, 0]]
        e2 = [p.fitmasserr for p in np.array(peaks)[:, 1]]
        self.simpleplotter(m1, m2, e1, e2, 0)

        m1 = [p.massavg for p in np.array(peaks)[:, 0]]
        m2 = [p.massavg for p in np.array(peaks)[:, 1]]
        e1 = [p.masserr for p in np.array(peaks)[:, 0]]
        e2 = [p.masserr for p in np.array(peaks)[:, 1]]
        self.simpleplotter(m1, m2, e1, e2, 1)

        m1 = [p.fitarea for p in np.array(peaks)[:, 0]]
        m2 = [p.fitarea for p in np.array(peaks)[:, 1]]
        e1 = [p.score for p in np.array(peaks)[:, 0]]
        e2 = [p.score for p in np.array(peaks)[:, 1]]
        print np.mean(e1)
        print np.mean(e2)
        self.simpleplotter(m1, m2, e1, e2, 2)

        m1 = [p.corrint for p in np.array(peaks)[:, 0]]
        m2 = [p.corrint for p in np.array(peaks)[:, 1]]
        e1 = [p.correrr for p in np.array(peaks)[:, 0]]
        e2 = [p.correrr for p in np.array(peaks)[:, 1]]
        self.simpleplotter(m1, m2, e1, e2, 3)

        plt.show()

    def plot_test(self, tests, vals, testdat, testmdat, peaks):
        import matplotlib.pyplot as plt
        print np.array(tests)
        np.set_printoptions(3, formatter={'float': '{: 0.2f}'.format})
        print np.array(vals)

        plt.figure(figsize=(14,6))
        offset = 100
        for i, d in enumerate(testdat):
            plt.subplot(121)
            plt.plot(d[:, 0], d[:, 1] * 100 - i * offset)
            plt.xlabel("m/z")
            plt.subplot(122)
            plt.xlabel("Mass")

            m = testmdat[i]
            plt.plot(m[:, 0], m[:, 1] - i * offset)

            for p in peaks[i]:
                plt.errorbar(p.fitmassavg, p.fitarea - i * offset, yerr=p.fitareaerr, xerr=p.fitmasserr, label="Fit",
                             linestyle="", color="r")
                plt.errorbar(p.massavg, p.corrint - i * offset, xerr=p.masserr, yerr=p.correrr, label="Corr",
                             linestyle="", color="b")
                plt.plot(p.mass, p.integral - i * offset, color="g", marker='o')
        plt.ylim(100 - (len(testdat) * 100), 100)
        plt.show()

    def test_width(self, **kwargs):
        reslist = [300, 500, 1000, 2000, 5000]
        tests = []
        vals = []
        testdat = []
        testmdat = []
        peaks = []
        for res in reslist:
            test, val, peak = self.test_spectra(res=res, noise=0.01, pad=2000, window=2000, plot=False, **kwargs)
            testdat.append(self.data.data2)
            testmdat.append(self.data.massdat)
            peaks.append(peak)
            tests.append(test)
            vals.append(np.average(val, axis=1))
        self.plot_test(tests, vals, testdat, testmdat, peaks)

    def test_noise(self, **kwargs):
        noiselist = [0, 0.001, 0.01, 0.1, 0.3]
        tests = []
        vals = []
        testdat = []
        testmdat = []
        peaks = []
        for n in noiselist:
            test, val, peak = self.test_spectra(res=1000, noise=n, pad=20000, window=2000, plot=False, **kwargs)
            testdat.append(self.data.data2)
            testmdat.append(self.data.massdat)
            peaks.append(peak)
            tests.append(test)
            vals.append(np.average(val, axis=1))
        self.plot_test(tests, vals, testdat, testmdat, peaks)

    def test_zwidth(self, **kwargs):
        zwlist = [0.5, 1, 3]
        tests = []
        vals = []
        testdat = []
        testmdat = []
        peaks = []
        for zw in zwlist:
            test, val, peak = self.test_spectra(res=1000, noise=0.01, pad=20000, window=2000, plot=False, zwidth=zw,
                                                **kwargs)
            testdat.append(self.data.data2)
            testmdat.append(self.data.massdat)
            peaks.append(peak)
            tests.append(test)
            vals.append(np.average(val, axis=1))
        self.plot_test(tests, vals, testdat, testmdat, peaks)

    def test_baseline(self, **kwargs):
        blist = [0, 0.01, 0.1, 0.3, 0.6, -0.6, -0.3, -0.1, -0.01]
        blist = [0,0.05,0.1]
        tests = []
        vals = []
        testdat = []
        testmdat = []
        peaks = []
        for b in blist:
            test, val, peak = self.test_spectra(res=1000, noise=0.01, pad=20000, window=2000, plot=False, baseline=b,
                                                **kwargs)
            testdat.append(self.data.data2)
            testmdat.append(self.data.massdat)
            peaks.append(peak)
            tests.append(test)
            vals.append(np.average(val, axis=1))
        self.plot_test(tests, vals, testdat, testmdat, peaks)

    def test_ambiguity(self, **kwargs):
        m1 = 200000
        m2list = [202000, 210000, 205000, 215000, 220000]
        tests = []
        vals = []
        testdat = []
        testmdat = []
        peaks = []
        for m2 in m2list:
            mlist = [m1, m2]
            test, val, peak = self.test_spectra(mlist=mlist, res=300, noise=0.01, pad=20000, window=1000, plot=False,
                                                **kwargs)
            testdat.append(self.data.data2)
            testmdat.append(self.data.massdat)
            peaks.append(peak)
            tests.append(test)
            vals.append(np.average(val, axis=1))
        self.plot_test(tests, vals, testdat, testmdat, peaks)

    def test_density(self, **kwargs):
        massbinlist = [3, 5, 10]
        tests = []
        vals = []
        testdat = []
        testmdat = []
        peaks = []
        num = 10
        for b in massbinlist:
            for n in xrange(num):
                test, val, peak = self.test_spectra(res=1000, noise=0.01, pad=20000, window=2000, plot=False,
                                                    massbins=b,
                                                    **kwargs)
                testdat.append(self.data.data2)
                testmdat.append(self.data.massdat)
                peaks.append(peak)
                tests.append(test)
                vals.append(np.average(val, axis=1))
                print n
        self.plot_test(tests, vals, testdat, testmdat, peaks)

    def test_psfun(self, **kwargs):
        psfuns = [0, 1, 2]
        tests = []
        vals = []
        testdat = []
        testmdat = []
        peaks = []
        num = 10
        for p in psfuns:
            for n in xrange(num):
                test, val, peak = self.test_spectra(res=1000, noise=0.01, pad=20000, window=2000, plot=False, psfun=p,
                                                    **kwargs)
                testdat.append(self.data.data2)
                testmdat.append(self.data.massdat)
                peaks.append(peak)
                tests.append(test)
                vals.append(np.average(val, axis=1))
        self.plot_test(tests, vals, testdat, testmdat, peaks)

    def test_widths(self, **kwargs):
        r1 = 1000
        r2list = [300]#, 500, 800, 1200, 3000]
        tests = []
        vals = []
        testdat = []
        testmdat = []
        peaks = []
        mlist = [200000, 202000]
        num = 1
        for r2 in r2list:
            rlist = [r1, r2]
            for n in xrange(num):
                test, val, peak = self.test_spectra(mlist=mlist, rlist=rlist, res=1000, noise=0.00, pad=2000,
                                                    window=1000,
                                                    plot=False,
                                                    psfun=0, plot_corr=False, plot_fits=False, **kwargs)
                testdat.append(self.data.data2)
                testmdat.append(self.data.massdat)
                peaks.append(peak)
                tests.append(test)
                vals.append(np.average(val, axis=1))
                print [p.score for p in self.pks.peaks]
        self.plot_test(tests, vals, testdat, testmdat, peaks)

    def test_sigs(self, **kwargs):

        rlist = [1000, 500]
        tests = []
        vals = []
        testdat = []
        testmdat = []
        peaks = []
        mlist = [200000, 202000]
        sigs = [1, 3, 5]
        num = 1
        for sig in sigs:
            for n in xrange(0, num):
                test, val, peak = self.test_spectra(mlist=mlist, rlist=rlist, res=1000, noise=0.01, pad=2000,
                                                    window=1000,
                                                    plot=False, mzsig=sig,
                                                    psfun=0, plot_corr=False, plot_fits=False, **kwargs)
                testdat.append(self.data.data2)
                testmdat.append(self.data.massdat)
                peaks.append(peak)
                tests.append(test)
                vals.append(np.average(val, axis=1))
                print [p.score for p in self.pks.peaks]
        self.plot_test(tests, vals, testdat, testmdat, peaks)

    def test_repeat(self, **kwargs):
        rlist = [1000, 2000]
        tests = []
        vals = []
        testdat = []
        testmdat = []
        peaks = []
        mlist = [200000, 210000]
        sigs = [1]
        num = 5
        for sig in sigs:
            for n in xrange(0, num):
                test, val, peak = self.test_spectra(mlist=mlist, rlist=rlist, res=1000, noise=0.01, pad=2000,
                                                    window=1000,
                                                    plot=False, mzsig=sig,
                                                    psfun=0, plot_corr=False, plot_fits=False, **kwargs)
                testdat.append(self.data.data2)
                testmdat.append(self.data.massdat)
                peaks.append(peak)
                tests.append(test)
                vals.append(np.average(val, axis=1))
                print [p.score for p in self.pks.peaks]
        self.plot_repeats(peaks)
        self.plot_test(tests, vals, testdat, testmdat, peaks)

    def test_chopped_spectrum(self,**kwargs):
        test, val, peak = self.test_spectra(mlist=[200000,202000],  res=1000,plot=True, **kwargs)


if __name__ == "__main__":
    multiprocessing.freeze_support()
    eng = UniDecTest()

    #eng.test_ambiguity()
    #eng.test_baseline()
    #eng.test_zwidth()
    # eng.test_density()
    # eng.test_psfun()
    #eng.test_widths()
    eng.test_chopped_spectrum()
    #eng.test_width()
    #eng.test_sigs()
    #eng.test_repeat()
    #eng.test_noise()
