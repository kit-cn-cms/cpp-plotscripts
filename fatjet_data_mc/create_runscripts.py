import glob
import os
import sys
import stat
import ROOT

program='PATH_TO_PROGRAM'

outpath='storage_folder'
scriptpath='...'
cmsswpath='...'


samples=[('ttbar',['/nfs/dust/cms/user/hmildner/trees1027/ttbar/ttbar*nominal*.root']),
         ('data',['...'])]

#events_per_job=1000000
files_per_job=20

if not os.path.exists(scriptpath):
    os.makedirs(scriptpath)

if not os.path.exists(outpath):
    os.makedirs(outpath)

if not os.path.exists(cmsswpath):
    print 'WRONG CMSSW PATH!'
    print cmsswpath
    sys.exit()

def getEventsInFiles(files):
    l=[]
    for f in files:
        tf=ROOT.TFile(f, 'readonly')
        tree=tf.Get('MVATree')
        l.append(tree.GetEntries())
        tf.Close()
    return l

       
def createScript(scriptname,programpath,processname,filenames,outfilename,maxevents,skipevents,suffix):
    script="#!/bin/bash \n"
    script+="export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch \n"
    script+="source $VO_CMS_SW_DIR/cmsset_default.sh \n"
    script+='cd '+cmsswpath+'/src\neval `scram runtime -sh`\n'
    script+='cd - \n'
    script+='export PROCESSNAME="'+processname+'"\n'
    script+='export FILENAMES="'+filenames+'"\n'
    script+='export OUTFILENAME="'+outfilename+'"\n'
    script+='export MAXEVENTS="'+str(maxevents)+'"\n'
    script+='export SKIPEVENTS="'+str(skipevents)+'"\n'
    script+='export SUFFIX="'+suffix+'"\n'
    script+=programpath+'\n'
    f=open(scriptname,'w')
    f.write(script)
    f.close()
    st = os.stat(scriptname)
    os.chmod(scriptname, st.st_mode | stat.S_IEXEC)

    
def createScriptsForSamples(samples,suffix,programpath):
    samples_=[]
    for s in samples:
        allfiles=[]
        neventsperfile=[]
        for f in s[1]:
            allfiles+=glob.glob(f)
            samples_.append((s[0],allfiles,neventsperfile))
    samples=samples_

    for sample in samples:
        process=sample[0]
        files=sample[1]
        ijob=1
        nfiles=0
        files_in_job=[]
        for f in files:
            files_in_job.append(f)
            nfiles+=1
            if nfiles>=files_per_job or f == files[-1]:
                scriptname=scriptpath+'/'+process+suffix+'_'+str(ijob)+'.sh'
                filenames=' '.join(files_in_job)
                outfilename=outpath+'/'+process+suffix+'_'+str(ijob)+'.root'
                maxevents=9999999999
                skipevents=0
                createScript(scriptname,programpath,process,filenames,outfilename,maxevents,skipevents,suffix)
                ijob+=1
                nfiles=0
                files_in_job=[]

createScriptsForSamples(samples,'',programpath)
