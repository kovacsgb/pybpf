import unittest
import os

try:
    from .. import tcdata as tcdata
    from .. import calcion as calcion
    from ..bpfreader import bpfDataRead
except (SystemError,ImportError,ValueError):
    #Testing purposes, later someone cleverer may remove this part
    import sys
    sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)),os.pardir))
    import tcdata
    import calcion
    from bpfreader import bpfDataRead
    





class TestBaseData(unittest.TestCase):
    
    def test_Constructor(self):
        data = {'label1' : 1, 'label2' : 2}
        NewBaseData = tcdata.BaseData(data,tuple(data.keys()))
        self.assertIsInstance(NewBaseData,tcdata.BaseData)
    
    def test_GetArg(self):
        data = {'label1' : 1, 'label2' : 2}
        NewBaseData = tcdata.BaseData(data,tuple(data.keys()))
        self.assertEqual(NewBaseData.data('label1'),1)
        self.assertEqual(NewBaseData.label2,2)
        self.assertRaises(KeyError,NewBaseData.data,'nonexistent')
        self.assertRaises(KeyError,NewBaseData.__getattr__,'nonexistent2')

    def test_InsertColumn(self):
        data = {'label1' : 1, 'label2' : 2}
        NewBaseData = tcdata.BaseData(data,tuple(data.keys()))
        NewBaseData = NewBaseData.insertColumn({'label3':3},tuple('label3'))
        self.assertEqual(NewBaseData.data('label3'),3)
        NewBaseData = NewBaseData.insertColumn([4],tuple(['label4']))
        self.assertEqual(NewBaseData.data('label4'),4)

class TestDataPoint(unittest.TestCase):
    def test_DataPointRead(self):
        a_point = tcdata.DataPoint("45\t136\t12\t50\t5.6\t2e5\t4.13e+29\t125.68\t1.25e40\t1.2e38\t12000\t1e29\t12e5\t1.5e4\t1e2\t0.98\n")
        self.assertEqual(a_point.dm,4.13e+29)

class TestRawProfiles(unittest.TestCase):
    def test_RawProfileInit(self):
        TEST_DIR=os.path.join(os.path.dirname(os.path.abspath(__file__)),'testfiles')
        TEST_FILE=TEST_DIR + '/fort.19'
        testraw = tcdata.RawProfiles(TEST_FILE)
        self.assertIsInstance(testraw,tcdata.RawProfiles)

class TestBPFReader(unittest.TestCase):
    TEST_DIR=os.path.join(os.path.dirname(os.path.abspath(__file__)),'testfiles')
    def test_default(self):
        mod,his,lim = bpfDataRead(self.TEST_DIR)
        self.assertIsInstance(mod,tcdata.Model)
        self.assertIsInstance(his,tcdata.History)
        self.assertIsInstance(lim,tcdata.LimitCycle)
    def test_with_rawprofile(self):
        mod,his,lim,raw =bpfDataRead(self.TEST_DIR,return_rawprofile = True)
        self.assertIsInstance(raw,tcdata.RawProfiles)
    def test_with_ionization(self):
        mod,his,lim = bpfDataRead(self.TEST_DIR,do_ionization=True,X=0.75,Y=0.246)
        self.assertIsInstance(mod,tcdata.Model)
        self.assertIsInstance(his,tcdata.History)
        self.assertIsInstance(lim,tcdata.LimitCycle)
    def test_ionization_error(self):
        self.assertRaises(RuntimeError,bpfDataRead,self.TEST_DIR,do_ionization=True)




if (__name__ == '__main__'):
    unittest.main(verbosity=2,exit=False)

    from line_profiler import LineProfiler
    print('======Profiling=======')
    lp = LineProfiler()
    lp.add_function(bpfDataRead)
    lp.add_function(tcdata.Model.__init__)
    lp.add_function(tcdata.History.__init__)
    lp.add_function(tcdata.LimitCycle.__init__)
    lp.add_function(tcdata.RawProfiles.CalcSpecVol)
    lp.enable_by_count()
    mod,his,lim = bpfDataRead(os.path.join(os.path.dirname(os.path.abspath(__file__)),'testfiles'))
    lp.print_stats()