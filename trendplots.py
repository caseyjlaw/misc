#!/usr/bin/env python

""" Script to use pyGTrends v 0.81.
Extracts trend info into csv format, parses, then plots relative popularities.
claw, 4dec10
"""

import sys, pylab, random, pyGTrends


def csvarray(pygt):
    """Takes csv output and creates numpy array.
    Assumes six terms input.
    """

    csv1 = (pygt.csv()).split('\n')
    date = []; c1 = []; c2 = []; c3 = []; c4 = []; c5 = []
    c6 = []; c7 = []; c8 = []; c9 = []; c10 = []; c11 = []; c12 = []; c13 = []

    for i in range(1,len(csv1)):
        fields = (csv1[i]).split(',')
        len(fields)
        try:
            c1.append(float(fields[1]))
            date.append(fields[0][:3] + ' \'' + fields[0][-2:])  # just take the month and year
            c2.append(float(fields[3]))
            c3.append(float(fields[5]))
            c4.append(float(fields[7]))
            c5.append(float(fields[9]))
            c6.append(float(fields[11]))
            c7.append(float(fields[13]))
            c8.append(float(fields[15]))
            c9.append(float(fields[17]))
            c10.append(float(fields[19]))
            c11.append(float(fields[21]))
        except ValueError:
            continue
        except IndexError:
            continue

    return date, (c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11)


def figs(username, password):
    """Generates figures of trends for pairs of trends.
    Plots one from list "goods" with one from list "celebs".  Each pair creates plot.
    """

    goods = ('martin luther king', 'muhammad yunus', 'mother teresa', 'the madonna', 'martti ahtisaari', 'gandhi', 'nelson mandela')
    celebs = ('kanye west', 'paris hilton', 'madonna', 'kim kardashian', 'jessica simpson', 'hannah montana', 'brad pitt')

    pygt = pyGTrends.pyGTrends(username, password)
    pygt.download_report(goods)
    date, rep1ar = csvarray(pygt)
#    goods = pygt.keywords
    xaxis1 = range(len(date))
    pygt = pyGTrends.pyGTrends(username, password)
    pygt.download_report(celebs)
#    celebs = pygt.keywords
    date, rep2ar = csvarray(pygt)
    xaxis2 = range(len(date))
    skip = 30
    
    for i in range(len(rep1ar)):
        if len(rep1ar[i]) == 0:  continue
        for j in range(len(rep2ar)):                             
            if len(rep2ar[j]) == 0:  continue
            pylab.plot(xaxis1, rep1ar[i], label=goods[i])
            pylab.plot(xaxis2, rep2ar[j], label=celebs[j])
            pylab.xticks(range(0, len(date), skip), (date[::skip]), rotation=60, size='x-small')
            pylab.ylabel('Popularity')
            pylab.legend(loc=2)
            pylab.title('Relative popularity plot # %i' % (random.randint(1,1000)))
            pylab.savefig('compare%i.png' % (j + i*len(goods)))
            pylab.clf()

if __name__ == "__main__":
    figs(sys.argv[1], sys.argv[2])
