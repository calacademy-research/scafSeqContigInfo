#!/usr/bin/env python

# looks through all scaffolds in scafSeq file (and potentially singletons is -C option used)
# and reports one line for each summarizing scaffold len, contig lens, num contigs, largest contig, N50, %Ns
# the output is sorted by the scaffold length largest to smalllest.
# By default only scaffold lengths >= 1K are shown, this can be changed to as low as 201 by -s option
#
# 10May15 JBH put gapThreshold logic in getContigInfo()
#             add -1 .. -99 option to allow 1 to 99 N's in a subsequence for it still to be considered a contig

# usage: scafSeqContigInfo.py <filename.scafSeq> 

import sys, time, re
from collections import Counter, defaultdict
from math import floor

aryScafContigInf = [] # each entry holds a tuple of:
    # scaf_len, scaf_name, contig_len_sum, num_contigs, largest_contig, N50, N50_contig_num, %Ns, numNs, smallContigBPs
    
cutOffLen = 1000 # only show info for scaffolds this length and greater
showScaffoldSummaryLines = True
showContigLines = True # show Longest and shortest contig lines
showTopContigScaffNames = False; minShowNameSize = 10000
gapThreshold = 0 # 10May15 number of N's allowable in a row to still be considered a contig

cntrContigLens = Counter()
mapLensToScaffs = defaultdict(list) # if showTopContigScaffNames we keep track of the scaffold that have lengths >= 10K
SMALL_CONTIG = 100
SMALL_SCAFFOLD = 200 # used to count small scaffolds and most importantly singleton Contigs

def scafProcess(fname):   
    fh = open(fname)
    
    current_milli_time = lambda: int(round(time.time() * 1000))
    lastTime = current_milli_time(); lastScaf = ""; newMsg = ""; msgLen = 0

    tstMax = 0 # set to 0 for no max
    totBPs = 0; totNs = 0; totContigs = 0;
    smallContigBPs = 0; scafsChecked = 0; scafsAdded = 0;
    smallScaffolds = 0; smallScaffoldBP = 0;
    
    nxthdr = fh.readline().rstrip()
    while nxthdr and nxthdr[0] == '>':
        nxthdr, info = getContigInfo(fh, nxthdr)
        if info[0] >= cutOffLen:
            aryScafContigInf.append(info)
            scafsAdded += 1
            totBPs += info[2]; totNs += info[8];
            totContigs += info[4]; smallContigBPs += info[9]
            tstMax = tstMax - 1
            if tstMax == 0:
                break
                
        if info[0] <= SMALL_SCAFFOLD:
            smallScaffolds += 1
            smallScaffoldBP += info[2]
                
        scafsChecked += 1
        if (current_milli_time() - lastTime) >= 1000: # update every second or so
            if lastScaf != info[0]:
                newMsg = str(scafsChecked) + " scaffolds processed ("+str(scafsAdded)+" ge cutoff), "
                newMsg += info[1] + " has " + str(info[2]) + " contig bases and N50 " + str(info[5])
                if len(newMsg) < msgLen:
                    newMsg = newMsg + ' '*(msgLen-len(newMsg))
                sys.stderr.write('\b'*msgLen + newMsg)
                msgLen = max(msgLen, len(newMsg)); lastTime = current_milli_time()
            lastScaf = info[0]
        
    fh.close()
    
    sys.stderr.write('\b'*msgLen + ' '*msgLen + '\n')
    displayResults(totBPs, totNs, totContigs, smallContigBPs, smallScaffolds, smallScaffoldBP)
   
def displayResults(totBPs, totNs, totContigs, smallContigBPs, smallScaffolds, smallScaffoldBP):
    aryScafContigInf.sort(reverse=True)
    items = [(k, v) for k, v in cntrContigLens.items()]
    contigCountByLength = sorted(items, reverse=True)
    
    midPoint = totBPs / 2
    
    #compute the contig N50 from the contigs in the processed scaffolds
    BPsofar = 0; N50 = 0; contigsSoFar = 0;
    for lngthAndcount in contigCountByLength: # looping from biggest contig lengths to smallest, see when we cross 50%
        BPinContigsThisSize = lngthAndcount[0] * lngthAndcount[1] # for example contigs of len 978 occured 14 time
        BPsofar += BPinContigsThisSize
        if BPsofar >= midPoint:
            N50 = lngthAndcount[0] # this is the length of the contig that sent us over the halfway point
            amountOver = BPsofar - midPoint # see how many of contigs this length were needed to go over midpoint
            numContigsNeeded = (amountOver / N50) + 1 # e.g. with 300bp contig we are 450 over, so 2 contigs needed
            contigsSoFar += numContigsNeeded
            break
        contigsSoFar += lngthAndcount[1]
    numScafs = len(aryScafContigInf)
    
    print "#", N50, "Contig N50 in", numScafs, "scaffolds >=", str(cutOffLen)+",",
    print "contig", contigsSoFar, "of", totContigs, "reached N50"

    #compute and display the N50 for the scaffolds
    totSoFar = 0; scaffsSoFar = 0; N50 = 0
    for inf in aryScafContigInf:
        totSoFar += inf[2] # add the whole scaffold contig length and see if we go over
        scaffsSoFar +=  1
        if totSoFar >= midPoint:
            N50 = inf[4] # take largest contig as the one that put us over (not a true N50)
            break
    
    pctNs = percentStr(totNs/((totBPs+totNs)*1.0)) # kind of a fluffyness metric
    contigsPerMb = totContigs / (totBPs / 1000000.0) # kind of a fragmentation metric
    print "#", "{:,}".format(totBPs), "contig bases, midpoint in", inf[1], "(largest contig", str(N50)+"),",
    print "scaffold #"+str(scaffsSoFar), "of", numScafs
    print "#", pctNs + " Ns in", numScafs, "scaffolds,", "{0:.2f}".format(contigsPerMb), "contigs per Mbases.",

    if showContigLines:
        if True: # 10May15 show gapThreshold info instead of small scaffold line
            print gapThreshold+1, "or more consecutive Ns end a contig."
        else:
            if smallScaffolds > 0:
                print "Scaffolds/Singletons len <=", str(SMALL_SCAFFOLD)+":", smallScaffolds, "summing to", smallScaffoldBP, "bases.",
            print
        showShortestContigs(contigCountByLength, 20)
        showLargestContigs(contigCountByLength, 20)
        showContigAvgs(contigCountByLength)
    else:
        print # finish up line with #N's
    
    if showScaffoldSummaryLines:
        #output a line for each scaffold
        for inf in aryScafContigInf:
            print inf[1], "{:,}".format(inf[0]), "len,", inf[2], "bases in", inf[3], "contigs,", 
            print "largest contig", str(inf[4])+",", "N50", inf[5], "in", inf[6], "contigs,", inf[7], "Ns"
        
def showLargestContigs(aryLengthAndCounts, numToShow, prefix="#"):
    maxShow = min(numToShow, len(aryLengthAndCounts))
    print prefix, maxShow, "longest contigs:",
    for ix in xrange(maxShow):
        lngthAndcount = aryLengthAndCounts[ix]
        if showTopContigScaffNames:
            print str(lngthAndcount[0]) + "(" + ' '.join(mapLensToScaffs[lngthAndcount[0]]) +")",
        else:
            if lngthAndcount[1] == 1:
                print lngthAndcount[0],
            else:
                print str(lngthAndcount[0]) + "(" + str(lngthAndcount[1]) + ")",
    print

def showContigAvgs(aryLengthAndCounts):
    contigList = longContigList(aryLengthAndCounts, 20);
    print "# ", len(contigList), "top contigs, Avg len:", getAvg(contigList), " Hmean:", int(getHmean(contigList))
    contigList = longContigList(aryLengthAndCounts, 100);
    print "#", len(contigList), "top contigs, Avg len: ", getAvg(contigList), " Hmean: ", int(getHmean(contigList))
  
def getAvg(numList):
    if len(numList) < 1:
        return 0
    return sum(n for n in numList) / len(numList)
    
def getHmean(numList):
    if len(numList) < 1:
        return 0
    return len(numList)/sum(1.0/n for n in numList)
        
def longContigList(aryLengthAndCounts, numContigs):
    maxContigs = min(numContigs, len(aryLengthAndCounts))
    numContigs = 0
    numList = []
    for ix in xrange(maxContigs):
        lngthAndcount = aryLengthAndCounts[ix]
        if lngthAndcount[1] == 1:
            numList.append( lngthAndcount[0] )
            numContigs += 1
        else:
            for c in xrange(lngthAndcount[1]):
                numList.append( lngthAndcount[0] )
                numContigs += 1
                if numContigs >= maxContigs:
                    break
        if numContigs >= maxContigs:
            break
    
    return numList    
        
def showShortestContigs(aryLengthAndCounts, numToShow, prefix="# "):
    toShow = min(numToShow, len(aryLengthAndCounts))
    print prefix + str(toShow) + " shortest contigs:",
    
    ix = len(aryLengthAndCounts)
    while toShow > 0:
        ix -= 1; toShow -= 1
        lngthAndcount = aryLengthAndCounts[ix]
        if lngthAndcount[1] == 1:
            print lngthAndcount[0],
        else:
            print str(lngthAndcount[0]) + "(" + str(lngthAndcount[1]) + ")",
    print
        
def getContigInfo(fh, hdr):
    ixSpace = hdr.find(" ") # this is what works for SOAP 2 scafSeq files
    if ixSpace != -1:
        scaffID = hdr[1:ixSpace]
    else:
        ixSpace = hdr.find("\t") # this is what works for SOAP 1 scafSeq files
        if ixSpace != -1:
            scaffID = hdr[1:ixSpace]
        else:
            scaffID = hdr[1:]

    ln = fh.readline().strip()
    if not ln:
        return ln, ()
            
    aryBasesRuns = []
    aryNruns = []
        
    actgRun = 0; Nrun = 0; NrunCount = 0; chTotal = 0; Ntotal = 0; smallContigBPs = 0; Nseq = 0
    
    while ln and ln[0] != '>':
        lnlen = len(ln) # 10May15 have to change from for stmt to while, since gapThreshold skips internally
        ix = 0
        while ix < lnlen:   # was: for ch in ln:
            ch = ln[ix]; ix += 1
            chTotal += 1
            if actgRun > 0:
                if ch != "N":
                    actgRun += 1
                else: # end of run or we might skip a few N's if gapThreshold not met
                    if gapThreshold > 0:
                        Nseq += 1; ixN = ix
                        while ixN < lnlen and ln[ixN] == "N" and Nseq <= gapThreshold:
                            Nseq += 1
                            ixN += 1
                        if Nseq <= gapThreshold:
                            actgRun += Nseq
                            ix = ixN
                            Nseq = 0
                            continue
                        Nseq = 0 # reset gapThreshold counter
                        
                    # end of this run report it
                    aryBasesRuns.append(actgRun)
                    cntrContigLens[actgRun] += 1 # count how many contigs of length actgRun we have in total set of scaffolds' contigs
                    if actgRun <= SMALL_CONTIG:
                        smallContigBPs += actgRun
                    elif actgRun >= minShowNameSize and showTopContigScaffNames:
                        mapLensToScaffs[actgRun].append(scaffID)
                    actgRun = 0
                    Nrun = 1
                    NrunCount += 1
            else:
                if ch == "N":
                    Nrun += 1
                else: # end of this run of N's report it
                    aryNruns.append(Nrun)
                    Ntotal += Nrun
                    Nrun = 0
                    actgRun = 1
            
        ln = fh.readline().strip()
    
    if actgRun > 0:
        aryBasesRuns.append(actgRun)
        cntrContigLens[actgRun] += 1
        if actgRun <= SMALL_CONTIG:
            smallContigBPs += actgRun
        elif actgRun >= minShowNameSize and showTopContigScaffNames:
            mapLensToScaffs[actgRun].append(scaffID)
        
    aryBasesRuns.sort(reverse=True)
    pct = Ntotal/ (chTotal*1.0)
    contigLen = chTotal-Ntotal
    
    N50, numContigs = contigN50(contigLen, aryBasesRuns)
    
    return ln, (chTotal, scaffID, contigLen, len(aryBasesRuns), aryBasesRuns[0], N50, numContigs, percentStr(pct), Ntotal, smallContigBPs)

def percentStr(val, digits=2):
    val *= 10 ** (digits + 2)
    return '{1:.{0}f}%'.format(digits, floor(val) / 10 ** digits) 

def contigN50(contigTotalLen, aryDescContigLens):
    midPoint = contigTotalLen / 2
    totSoFar = 0
    contigsSoFar = 0
    N50contig = 0
    for contig_len in aryDescContigLens:
        totSoFar += contig_len
        contigsSoFar +=  1
        if totSoFar >= midPoint:
            N50contig = contig_len
            break
            
    return N50contig, contigsSoFar

def toInt(intStr, default=1000): # allow 1M, 10K, 100K for numbers as well as any integer
    if intStr.isdigit():
        return int(intStr)
    elif intStr[0].isdigit():
        intStr = intStr.replace("M", "000000") # so 1M becomes 1000000
        intStr = intStr.replace("K", "000") # so 10K becomes 10000
    if intStr.isdigit():
        return int(intStr)
    return default # return the default cutOffLen
    
def usage():        
    print("\nUsage: scafSeqContigInfo.py <filename.scafSeq> [ [-x] [-s <smallest_scaff_to_include>] | -All ] [-1 .. -99]\n\n"
          "       Looks through scaffolds in the scafSeq file (by default those 1K or larger) and displays overview info,\n"
          "       By default reports one line for scaffold, summarizing scaffold len, contig lens, num contigs, largest contig, N50, %Ns\n"
          "       the output is sorted by the scaffold length largest to smallest (this length is the contig bases + N gaps).\n\n"
          "       -s option: by default only scaffold lengths >= 1K are shown, this can be changed to as low as 201 by -s option\n"
          "       -S option: same as -s and list of Longest Contigs will include Scaffold names\n"
          "       -x option: this excludes the individual info lines for each scaffold and only shows the overview info lines\n"
          "       -nc option: do not show longest and shortest contig info lines\n"
          "       -All option: shows overview info for all scaffolds and singleton contigs, individual info lines are not shown\n"
          "       -1 .. -99, e.g. -4, sets gap threshold to number of N's allowed in contig before splitting it (0 by default)\n"
          )
    sys.exit(0)
    
def main(): 
    global cutOffLen, showScaffoldSummaryLines, showContigLines, showTopContigScaffNames, gapThreshold
    if len(sys.argv) >= 2:
        filename = ""
        maxArgs = len(sys.argv)-1
        ixArg = 0
        while ixArg < maxArgs:
            ixArg += 1
            arg = sys.argv[ixArg]
            if arg[0] != '-' and filename == "":
                filename = arg
            elif arg == "-x": # exclude the individual summary lines for each scaffold, just show overview lines
                showScaffoldSummaryLines = False
            elif arg == "-All": # look at all scaffolds and singletons for overview  info but don't output individual info lines
                cutOffLen = 0
                showTopContigScaffNames = True
                showScaffoldSummaryLines = False    
            elif arg == "-nc": # no Contig longest and shortest lines
                showContigLines = False
            elif arg == "-s" or arg == "-S":
                if arg == "-S":
                    showTopContigScaffNames = True
                ixArg += 1
                arg = sys.argv[ixArg]
                cutOffLen = toInt(arg)
                if cutOffLen < 201:
                    cutOffLen = 201
            elif arg[0] == '-' and arg[1:].isdigit(): #10May15 if -1 .. -99 set gap threshold, so small number of N's won't split a contig
                gapThreshold = int(arg[1:])

        if filename != "":
            scafProcess(filename)
        else:
            usage()
    else:
        usage()

main()