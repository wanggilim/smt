#! /usr/bin/env python
import argparse
import logging
import sys

# setup logging
log = logging.getLogger(__file__)

def itime3(noddwell, repeats, rewinds=0, dithers=0, loops=1,
           chopeff=1/1.43, c2nc2=False, nxcac=False,
           nodtime=5, dithertime=5, filtchange=30, lost=25,
           *args,**kwargs):
    '''Does not take into account possibly different LOS rewind time for dither vs.
    non-dither
    If no rewind number needed, rewind should be set to 0
    If no dithers, dithers should be set to 0
    Duration per rewind only shown if rewind not set equal to 0
    If nxcac is indicated, the grism keyword is ignored
    '''

    repeats = int(repeats)
    rewinds = int(rewinds)
    if rewinds < 0:
        rewinds = 0
    dithers = int(dithers)
    loops = int(loops)

    if (dithers == 0):
        if c2nc2:
            #print('C2NC2 mode must have dithers')
            #return
            msg = 'C2NC2 mode must have dithers'
            log.error(msg)
            raise ValueError(msg)
        if rewinds == 0:
            log.info('No dithers and no rewinds during AOR')
            dur = repeats * ((noddwell*2)+nodtime) + nodtime * (repeats % 2)
            durperrew = 0  # dur?
            fdur = dur + filtchange
        else:
            log.info('No dithers')
            if (repeats % rewinds) > 0:
                dur1 = rewinds * ((noddwell*2)+nodtime) + nodtime*(rewinds % 2) + lost
                durperrew = dur1
                numlos = int(repeats/rewinds)
                fdur1 = numlos*dur1
                remain = repeats - (rewinds*numlos)
                fdur2 = remain*((noddwell*2) + nodtime) + nodtime*(rewinds % 2) + lost
                fdur = fdur1 + fdur2 + filtchange
            else:
                dur= rewinds*((noddwell*2.0)+nodtime)  + nodtime*(rewinds % 2)
                durperrew = dur + lost
                numlos = int(repeats/rewinds)
                fdur = numlos*durperrew + filtchange

    else:
        if c2nc2:
            log.info('C2NC2 mode')
            if (repeats != 0) or (rewinds != 0):
                #print('Repeats and rewinds must be zero in C2NC2 mode')
                #return
                msg = 'Repeats and rewinds must be zero in C2NC2 mode'
                log.error(msg)
                raise ValueError(msg)
            
            num_aba = int(int(dithers)/2)
            num_ab = dithers % 2
            dur = (num_aba*noddwell*3.0) + (num_aba*nodtime*2.0) + \
                  (num_ab *noddwell*2.0) + (num_ab *nodtime*2.0)
            durperrew = dur + lost
            if loops == 1:
                log.info('%i loop'%loops)
            else:
                log.info('%i loops'%loops)
                
            fdur = (durperrew*loops) + filtchange
            
        else:
            if repeats == 0:
                #print('C2N and NXCAC modes requires at least 1 repeat.')
                #return
                msg = 'C2N and NXCAC modes requires at least 1 repeat.'
                log.error(msg)
                raise ValueError(msg)
            else:
                if repeats < rewinds:
                    if (rewinds % repeats) > 0:
                        #print('Rewinds must be an integer multiple of repeats')
                        #return
                        msg = 'Rewinds must be an integer multiple of repeats'
                        log.error(msg)
                        raise ValueError(msg)
                    else:
                        if rewinds == 0:
                            numlos = 1
                        else:
                            numlos = int(repeats*dithers/rewinds)
                    dur = repeats*((noddwell*2.0)+nodtime)  + nodtime*(repeats % 2)
                    if (((dithers*repeats) % rewinds) != 0.0):
                        fdur = dithers*(dur+dithertime) + (numlos+1)*lost + filtchange
                    else:
                        fdur=dithers*(dur+dithertime) + numlos*lost + filtchange
                    if repeats == rewinds:
                        log.info('Rewind when moving to next dither position')
                        durperrew = dur + dithertime + lost
                    else:
                        log.info('Dither before rewind')
                        durperrew = (rewinds/repeats)*(dur+dithertime)

                else:
                    if rewinds == 0:
                        log.info('No rewinds during AOR')
                        fdur = dithers*(repeats*((noddwell*2)+nodtime) + \
                                        nodtime*(repeats % 2) + dithertime) + filtchange
                        dur = fdur
                        #print(dur)
                        durperrew = 0
                    else:
                        if (repeats % rewinds) > 0:
                            #print('Repeats must be an integer multiple of rewinds')
                            #return
                            msg = 'Repeats must be an integer multiple of rewinds'
                            log.error(msg)
                            raise ValueError(msg)
                        else:
                            log.info('Rewind before dither')
                            dur = rewinds*((noddwell*2.0)+nodtime)  + nodtime*(rewinds % 2) + lost
                            durperdith = (repeats/rewinds)*dur + dithertime
                            fdur = durperdith*dithers +filtchange
                            durperrew = dur




    # Print durations and int times to screen
    log.info('Duration per rewind:\t%is or %.1fm' % (int(durperrew), durperrew/60))
    log.info('Total AOR duration:\t%is or %.1fm' % (int(fdur), fdur/60))

    if dithers == 0:
        dithers = 1
    if rewinds == 0:
        rewinds = 1

    '''
    if grism or nxcac:
        eff = 3
    else:
        eff = 1.43
    '''

    eff = 1/chopeff

    if c2nc2:
        inttime= ((noddwell/2.)*dithers*loops)/eff
        log.info('AOR effective on-source integration time: %is'%inttime)
        log.info('Observation efficiency: %.1f%%' % ((inttime/fdur)*100))

    else:
        if nxcac:
            inttime=((noddwell/2)*dithers*repeats)/eff
            log.info('AOR effective on-source integration time: %is'%inttime)
            log.info('Observation efficiency: %.1f%%' % ((inttime/fdur)*100))
        else:
            inttime=(noddwell*2*dithers*repeats)/eff
            log.info('AOR effective on-source integration time: %is'%inttime)
            log.info('Observation efficiency: %.1f%%' % ((inttime/fdur)*100))


    # (duration per rewind, total aor duration, on-source integration, efficiency
    return durperrew,fdur,inttime,inttime/fdur



def main():
    parser = argparse.ArgumentParser(description='Integration time efficiency calculator')
    parser.add_argument('noddwell',type=float,help='Nod dwell time in sec')
    parser.add_argument('repeats',type=int,help='Number of repeats')
    parser.add_argument('rewinds',type=int,default=0,nargs='?',help='Number of repeats per rewind e.g. rewind cadence (default=0; no rewinds during AOR)')
    parser.add_argument('dithers',type=int,default=0,nargs='?',help='Number of dithers (default=0)')
    parser.add_argument('-loops',type=int,default=1,help='Number of loops for C2NC2 mode (default=1)')
    parser.add_argument('-chopeff',type=float,default=0.7,help='Chop efficiency (default=0.70)')
    parser.add_argument('--c2nc2',action='store_true',help='C2NC2 mode')
    parser.add_argument('--nxcac',action='store_true',help='NXCAC mode')
    parser.add_argument('-nodtime',type=float,default=5,help='Nod overhead in sec (default=5)')
    parser.add_argument('-dithertime',type=float,default=5,help='Dither overhead in sec (default=5)')
    parser.add_argument('-filtchange',type=float,default=30,help='Filter change overhead in sec (default=30)')
    parser.add_argument('-lost',type=float,default=25,help='Additional overhead in sec (default=25)')

    args = parser.parse_args()

    logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(message)s')

    durperrew,fdur,inttime,eff = itime3(args.noddwell,args.repeats,
                                        args.rewinds,args.dithers,
                                        loops=args.loops,
                                        chopeff=args.chopeff,c2nc2=args.c2nc2,
                                        nxcac=args.nxcac,nodtime=args.nodtime,
                                        dithertime=args.dithertime,
                                        filtchange=args.filtchange,
                                        lost=args.lost)

    

if __name__ == '__main__':
    main()
