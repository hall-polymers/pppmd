import kshen.ion_dynamics
import pandas as pd

def main():
    intrvl = 38
    nBlock = 152
    blockSize = 39
    timescale = 0.005
    firstshell = 1.4
    secondshell = 2.4

    r, ir, timestep, box_bounds, id2type, id2index = special_read('loglog.lammpstrj', types=[3,4], num_frames=((nBlock-1)*intrvl + blockSize))
    r, boxsize, bound_lo = wrap(r, box_bounds)
    bin2id, id2bin, bins = binbox(r, boxsize, bound_lo, dist_range=firstshell)
    id2neighbors_1stshell = buildnlist(r, bin2id, id2bin, bins, boxsize, id2type, dist_range=firstshell, nearest=False)    
    bin2id, id2bin, bins = binbox(r, boxsize, bound_lo, dist_range=secondshell)
    id2neighbors_2ndshell = buildnlist(r, bin2id, id2bin, bins, boxsize, id2type, dist_range=secondshell, nearest=False)

    autocorr = np.zeros([nBlock, 4, blockSize], np.float)
    avg_num = np.zeros([nBlock, blockSize], np.float)
    nfreeions = np.zeros([nBlock, blockSize], np.float)
    nSSIP = np.zeros([nBlock, blockSize], np.float)
    nCIP = np.zeros([nBlock, blockSize], np.float)
    nAGG = np.zeros([nBlock, blockSize], np.float)
    for n in range(nBlock):
        lo = n*intrvl
        hi = n*intrvl + blockSize
        print 'Calculating cluster autocorrelation of block' , n, 'with timesteps between:', timestep[lo], timestep[hi-1]
        autocorr[n], avg_num[n], nfreeions[n], nSSIP[n], nCIP[n], nAGG[n] = buildclusterlist(r[lo:hi], id2neighbors_1stshell[lo:hi], id2neighbors_2ndshell[lo:hi])
    pd.DataFrame(zip(timestep[:blockSize]*timescale, autocorr.mean(axis=0)[0], autocorr.mean(axis=0)[1], autocorr.mean(axis=0)[2], autocorr.mean(axis=0)[3], autocorr.std(axis=0)[0], autocorr.std(axis=0)[1], autocorr.std(axis=0)[2], autocorr.std(axis=0)[3])).to_excel('ionpair_lifetime.xlsx'.format(nBlock), header=['Time', 'Avg_FI', 'Avg_SSIP', 'Avg_CIP', 'Avg_AGG', 'Std_FI', 'Std_SSIP', 'Std_CIP', 'Std_AGG'], index=False)

    print '============Cluster stats============'
    print 'Average cluster size =', avg_num.mean(), avg_num.std()
    print 'Averages (Freeions, SSIP, CIP, AGG):'
    print 'Standard deviation (Freeions, SSIP, CIP, AGG):'
    print nfreeions.mean(), nSSIP.mean(), nCIP.mean(), nAGG.mean()
    print nfreeions.std(), nSSIP.std(), nCIP.std(), nAGG.std()

if __name__ == '__main__':    
    sys.exit(main())