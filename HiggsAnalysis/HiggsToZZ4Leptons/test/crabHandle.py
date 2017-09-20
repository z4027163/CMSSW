#!/usr/bin/env python

import sys
import os

try:
    from CRABAPI.RawCommand import crabCommand as cc
    import CRABClient.UserUtilities
    import CRABClient.ClientUtilities
#    CRABClient.UserUtilities.setConsoleLogLevel(CRABClient.ClientUtilities.LOGLEVEL_MUTE)
except ImportError:
    print 'Must source CRAB setup script first!'
    raise

try:
    pwd=sys.argv[1]
    os.chdir(pwd)
except IndexError:
    pwd=os.getcwd()

try:
    dirs=sys.argv[2:]
    dirs=[[y for y in x.split('/') if 'crab' in y][-1] for x in dirs]
except IndexError:
    dirs=[direct for direct in os.listdir(pwd) if (os.path.isdir(direct))]# and os.path.isfile(direct+'_ND'))]

if not dirs:
    print 'All crab jobs finished!'
    sys.exit(0)

for directory in dirs:
    try:
        print 'Running crab status on directory %s' % directory
        results=cc('status',dir=directory)
        statuses=results['jobsPerStatus']
        if 'failed' in statuses: #Want to resubmit if failed:
            if 'finished' in statuses or 'transferring' in statuses: #But only if some of them have NOT failed
                print 'Resubmitting...'
                cc('resubmit',dir=directory)
        if results['status']=='COMPLETED':
#            print "Fetching output..."
            sys.exit(0)
#            if statuses['completed'] > 500:
#                jmin=1
#                jmax=500
#                while jmin<statuses['completed']:
#                    cc('getoutput',dir=directory,jobids='%i-%i'%(jmin,jmax))
#                    jmin+=500
#                    jmax+=500
#                sys.exit(0)
#            else:
#                go=cc('getoutput',dir=directory)
#                print go
#                os.remove(directory+'_ND')
    except Exception:
        print "Exception caught and ignored in directory %s" % directory
        raise

sys.exit(1)
