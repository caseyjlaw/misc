import multiprocessing as mp
import numpy as n
import time

def scheme1():
    """ Uses two threads optionally called in while-1 loop.
    Works!
    """

    def gendata(d, l):
        """ Data generating stage of parallel data function.
        data0 is where data is held to be processed. 
        if reading gets ahead, data1 is used for next chunk.
        """

        name = mp.current_process().name
        n.random.seed()
        ts = n.random.random()
        print 'gendata %s needs %.1f seconds' % (name, ts)
        time.sleep(ts)
        with l:
            kk = d.keys()
            if 'data0' not in kk:   # if data0 empty, fill it in one of two ways
                if 'data1' not in kk:
                    d['data0'] = n.random.random(10)
                    print 'filling new data0'
                else:
                    d['data0'] = d.pop('data1')
                    print 'moving to data0'
            elif 'data1' not in kk:  # once data0 is ok, then consider filling data1
                d['data1'] = n.random.random(10)
                print 'filling new data1'

    def processdata(d, l):
        """ Processing stage of parallel data function. 
        Only processes from data0. Uses a "first in, first out" model, where data0 is next in line.
        """


        name = mp.current_process().name
        n.random.seed()
        ts = n.random.random()
        print 'processdata %s needs %.1f seconds' % (name, ts)
        time.sleep(ts)
        kk = d.keys()
        with l:
            if 'data0' in kk:
                d['mean'] = d.pop('data0').mean()

    def start():
        """ Threading for parallel data reading and processing.
        Either side can be faster than the other, since data are held for processing in shared buffer.
        """

        mgr = mp.Manager()
        d = mgr.dict()
        l = mp.Lock()
        pgenlist = []; pproclist = []
        try:
            while 1:
                kk = d.keys()
#                print kk, len(pgenlist), len(pproclist)

                if ((('data0' not in kk) or ('data1' not in kk)) and (len(pgenlist) <= 1)):
                    pgen = mp.Process(target=gendata, args=(d,l))
                    pgen.start()
                    pgenlist.append(pgen)
                if (('data0' in kk) and (len(pproclist) == 0)):
                    pproc = mp.Process(target=processdata, args=(d,l))
                    pproc.start()
                    pproclist.append(pproc)

                # process cleanup
                    for proc in pgenlist:
#                        print 'checking proc', proc, proc.is_alive()
                        if not proc.is_alive():
#                            print 'joining proc', proc
                            proc.join()
                            pgenlist.pop(pgenlist.index(proc))
                    for proc in pproclist:
#                        print 'checking proc', proc, proc.is_alive()
                        if not proc.is_alive():
#                            print 'joining proc', proc
                            proc.join()
                            pproclist.pop(pproclist.index(proc))
                            print 'mean = ', d['mean']

        except KeyboardInterrupt:
            print 'Ctrl-C received. Shutting down threads...'
            for pgen in pgenlist:
                pgen.terminate()
                pgen.join()
            for pproc in pproclist:
                pproc.terminate()
                pproc.join()


def scheme2():
    """ Tries to mimic reading data with an iterator triggered from outside, a la msiter.
    """

    def gendata(d, l, e):
        """ Data generating stage of parallel data function.
        data0 is where data is held to be processed. 
        if reading gets ahead, data1 is used for next chunk.
        reading triggered by event
        """

        data = iter(range(10))
        while 1:
            e.wait()
            try:
                name = mp.current_process().name
                n.random.seed()
                ts = n.random.random()
                ts = 0.2
                print 'gendata %s needs %.1f seconds' % (name, ts)
                time.sleep(ts)
                with l:
                    kk = d.keys()
                    if 'data0' not in kk:   # if data1 empty, fill it
                        if 'data1' not in kk:
                            d['data0'] = data.next()
                            print 'filling new data0'
                        else:
                            d['data0'] = d.pop('data1')
                            print 'moving to data0'
                    elif ('data1' not in kk):  # if data0 is ok, fill from data1, if filled
                        d['data1'] = data.next()
                        print 'filling new data1'
            except StopIteration:
                print 'End of iteration'
                break
            e.clear()

    def processdata(d, l, e):
        """ Processing stage of parallel data function. 
        Only processes from data0. Uses a "first in, first out" model, where data0 is next in line.
        Needs event signal to process.
        """

        while 1:
            e.wait()
            name = mp.current_process().name
            n.random.seed()
            ts = n.random.random()
            ts = 0.1
            print 'processdata %s needs %.1f seconds' % (name, ts)
            time.sleep(ts)
            kk = d.keys()
            with l:
                if 'data0' in kk:
                    d['mean'] = float(d.pop('data0'))
            e.clear()
            print 'mean = ', d['mean']

    def start():
        """ Threading for parallel data reading and processing.
        Either side can be faster than the other, since data are held for processing in shared buffer.
        """

        mgr = mp.Manager()
        d = mgr.dict()
        l = mp.Lock()
        egen = mp.Event()
        eproc = mp.Event()

        try:
            # start processes
            pgen = mp.Process(target=gendata, args=(d,l,egen))
            pgen.start()
            pproc = mp.Process(target=processdata, args=(d,l,eproc))
            pproc.start()

            # loop monitors data areas to determine whether to trigger more reading or processing
            while 1:
                kk = d.keys()
                if ('data0' not in kk) or ('data1' not in kk):
                    egen.set()
                if 'data0' in kk:
                    eproc.set()

        except KeyboardInterrupt:
            print 'Ctrl-C received.'

        finally:
            print 'Shutting down threads...'
            pgen.terminate()
            pgen.join()
            pproc.terminate()
            pproc.join()

    start()
