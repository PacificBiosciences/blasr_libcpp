#ifndef _BLASR_FIND_MAX_INTERVAL_IMPL_HPP_
#define _BLASR_FIND_MAX_INTERVAL_IMPL_HPP_

template<typename T_Sequence, typename T_AnchorList>
float DefaultWeightFunction<T_Sequence, T_AnchorList>::operator() (
    T_Sequence &text, T_Sequence &read, T_AnchorList matchPosList) {
    int i;
    float weight = 0;
    for (i = 0; i < matchPosList.size(); i++) {
        weight += matchPosList[i].weight();
    }
    return weight;
}

template<typename T_Pos>
int MatchPosQueryOrderFunctor<T_Pos>::operator()(T_Pos &pos) {
    return pos.q;
}

template<typename T_MatchList>
void PrintLIS(T_MatchList &matchList, DNALength curPos, 
    DNALength curGenomePos, DNALength nextGenomePos, 
    DNALength clp, DNALength cle) {
    int i;
    cout << curPos << " " << curGenomePos << " " 
        << nextGenomePos << " " << clp << " " << cle << endl;
    for (i = 0; i < matchList.size(); i++) {
        cout.width(8);
        cout << matchList[i].l << " ";
    }
    cout << endl;
    for (i = 0; i < matchList.size(); i++) {
        cout.width(8);
        cout << matchList[i].q << " ";
    }  
    cout << endl;
    for (i = 0; i < matchList.size(); i++) {
        cout.width(8);
        cout << matchList[i].t << " ";
    }
    cout << endl << endl;
}


template<typename T_MatchList, typename T_SequenceDB>
void FilterMatchesAsLIMSTemplateSquare(T_MatchList &matches, 
    DNALength queryLength, DNALength limsTemplateLength, 
    T_SequenceDB &seqDB) {
    int seqIndex;
    //
    // Make sure there is sequence coordinate information.
    //
    if (seqDB.nSeqPos == 0) {
        return;
    }
    int matchIndex = 0;
    for (seqIndex = 1; seqIndex < seqDB.nSeqPos; seqIndex++) {
        DNALength refLength = seqDB.seqStartPos[seqIndex] - 
            seqDB.seqStartPos[seqIndex - 1];
        // account for indel error.
        refLength = queryLength * 1.15 + limsTemplateLength; 

        //
        // Flag matches that are beyond the (rough) square with the length
        // of the query for removal.
        //
        while (matchIndex < matches.size() and 
                matches[matchIndex].t < seqDB.seqStartPos[seqIndex]) {
            if (matches[matchIndex].t > seqDB.seqStartPos[seqIndex-1] 
                + refLength) {
                matches[matchIndex].l = 0;
            }
            matchIndex++;
        }
        int curMatchIndex = 0;
        matchIndex = 0;
        for (matchIndex = 0; matchIndex < matches.size(); matchIndex++) {
            if (matches[matchIndex].l != 0) {
                matches[curMatchIndex] = matches[matchIndex];
                curMatchIndex++;
            }
        }
        matches.resize(curMatchIndex);
    }
}

template<typename T_MatchList, typename T_SequenceBoundaryDB>
void AdvanceIndexToPastInterval(T_MatchList &pos, DNALength nPos,
    DNALength intervalLength, DNALength contigLength,
    T_SequenceBoundaryDB &SeqBoundary,
    DNALength startIndex, DNALength startIntervalBoundary,
    DNALength &index, DNALength &indexIntervalBoundary) {
    (void)(contigLength);
    if (index >= pos.size()) {
        return;
    }
    indexIntervalBoundary = SeqBoundary(pos[index].t);
    DNALength boundaryIndex = SeqBoundary.GetIndex(pos[index].t);
    DNALength nextBoundary  = SeqBoundary.GetStartPos(boundaryIndex + 1);
    (void)(nextBoundary);
    while (// index is not past the end of the genome
            index < nPos and 
            //
            // Stop when the index goes too far ahead.
            //
            pos[index].t - pos[startIndex].t <= intervalLength and
            //
            // Still searching in the current contig.
            //
            indexIntervalBoundary == startIntervalBoundary) {
        index++;
        if (index < nPos) {
            indexIntervalBoundary = SeqBoundary(pos[index].t);
        }
    }
}


template<typename T_MatchList>  
int RemoveZeroLengthAnchors(T_MatchList &matchList) {       
    int origSize = matchList.size();
    int cur = 0, m;
    for (m = 0; m < matchList.size(); m++) {
        if (matchList[m].l > 0) {
            matchList[cur] = matchList[m];
            cur++;
        }
    }
    matchList.resize(cur);
    return origSize - cur;
}

template<typename T_MatchList>
int RemoveOverlappingAnchors(T_MatchList &matchList) {       
  int m;
  int n;
  for (m = int(matchList.size()); m > 0; m--) {
    n = m - 1;
    //
    // Skip past repeats in the query.
    while (n > 0 and matchList[n].t == matchList[m].t) {
      n--;
    }

    bool mergeFound = false;
    int ni = n;
    while (mergeFound == false and 
           n > 0 and 
           matchList[n].t == matchList[ni].t) {
      if (matchList[n].q < matchList[m].q and
          matchList[n].t < matchList[m].t and
          matchList[n].l + matchList[n].q >= 
          matchList[m].l + matchList[m].q and
          matchList[n].l + matchList[n].t >= 
          matchList[m].l + matchList[m].t) {
        matchList[m].l = 0;
        mergeFound = true;
      }
      n--;
    }
  }
  int numRemoved = RemoveZeroLengthAnchors(matchList);
  return numRemoved;
}

template<typename T_MatchList>
int SumAnchors(T_MatchList &pos, int start, int end) {
    int sum = 0;
    int i;
    for (i = start; i < end; i++) {
        sum += pos[i].l;
    }
    return sum;
}

template<typename T_MatchList,
         typename T_SequenceBoundaryDB>
void StoreLargestIntervals(
    T_MatchList &pos, 
    // End search for intervals at boundary positions
    // stored in seqBoundaries
	T_SequenceBoundaryDB & ContigStartPos,
	// parameters
	// How many values to search through for a max set.
	DNALength intervalLength,  
	// How many sets to keep track of
	int minSize,
	vector<DNALength> &start,
	vector<DNALength> &end) {

    if (pos.size() == 0) {
        return;
    }
    //
    // Search for clusters of intervals within the pos array within
    // pos[cur...next).  The value of 'next' should be the first anchor
    // outside the possible range to cluster, or the end of the anchor list.
    VectorIndex cur = 0;
    VectorIndex nPos = pos.size();
    VectorIndex next = cur + 1;
    DNALength curBoundary = 0, nextBoundary = 0;
    DNALength contigLength = ContigStartPos.Length(pos[cur].t);
    DNALength endOfCurrentInterval = curBoundary + contigLength;
    (void)(curBoundary); (void)(nextBoundary); (void)(endOfCurrentInterval);

    curBoundary = ContigStartPos(pos[cur].t);
    nextBoundary = ContigStartPos(pos[next].t);  

    //
    // Advance next until the anchor is outside the interval that
    // statrts at 'cur', and is inside the same contig that the anchor
    // at cur is in.
    //

    DNALength curIntervalLength = NumRemainingBases(pos[cur].q, intervalLength);
    (void)(curIntervalLength);

    AdvanceIndexToPastInterval(pos, nPos, intervalLength, contigLength, ContigStartPos,
            cur, curBoundary, next, nextBoundary);

    DNALength maxStart = cur, maxEnd = next;
    int maxSize = SumAnchors(pos, cur, next);
    int curSize = maxSize;
    if (curSize > minSize) {
        start.push_back(cur);
        end.push_back(next);
    }
    while ( cur < nPos ) {
        // 
        // This interval overlaps with a possible max start
        //

        if (pos[cur].t >= pos[maxStart].t and maxEnd > 0 and pos[cur].t < pos[maxEnd-1].t) {
            if (curSize > maxSize) {
                maxSize = curSize;
                maxStart = cur;
                maxEnd   = next;
            }
        }
        else {
            if (maxSize > minSize) {
                start.push_back(maxStart);
                end.push_back(maxEnd);
            }
            maxStart = cur;
            maxEnd   = next;
            maxSize  = curSize;
        }

        //
        // Done scoring current interval.  At this point the range
        // pos[cur...next) has been searched for a max increasing
        // interval.  Find a new range that will possibly yield a new
        // maximum interval.  
        // There are a few cases to consider:
        //
        //
        // genome  |---+----+------------+------+-----------------------|
        //  anchors  cur  cur+1        next   next+1
        //
        // Case 1.  The range on the target pos[ cur+1 ... next].t is a
        // valid interval (it is roughly the length of the read).  In this
        // case increase cur and next by 1, and search this range.
        //
        // genome  |---+----+------------+------+-----------------------|
        //            cur  cur+1        next   next+1
        // read interval   ====================
        //
        // Case 2.  The range on the target pos[cur+1 ... next] is not a
        // valid interval, and it is much longer than the length of the
        // read.  This implies that it is impossible to increase the score
        // of the read by including both 
        //
        // genome  |---+----+--------------------------------+-----+---|
        //            cur  cur+1                             next next+1
        // read interval   ==================== 
        //
        // Advance the interval until it includes the next anchor
        //
        // genome  |---+----+------------------+-------------+-----+---|
        //            cur  cur+1             cur+n          next next+1
        // read interval                     ==================== 
        // 

        // First advance pointer in anchor list.  If this advances to the
        // end, done and no need for further logic checking (break now).


        // 
        // If the next position is not within the same contig as the current,
        // advance the current to the next since it is impossible to find
        // any more intervals in the current pos.
        //
        bool recountInterval = false;
        if (curBoundary != nextBoundary) {
            cur = next;
            curBoundary = nextBoundary;

            //
            // Start the search for the first interval in the next contig
            // just after the current position.
            //
            if (next < nPos) {
                next = cur + 1;
            }
        }
        else {

            //
            // The next interval is in the same contig as the current interval.
            // Make sure not to double count the current interval.
            //

            curSize -= pos[cur].l;

            cur++;
            if (cur >= nPos)
                break;

            //
            // Advance the next to outside this interval.
            //
            curSize += pos[next].l;
            if (pos[next].t - pos[cur].t > intervalLength) {
                if (maxSize > minSize) {
                    start.push_back(maxStart);
                    end.push_back(maxEnd);
                }

                cur = next;
                recountInterval = true;
                maxSize = 0;
            }
            next++;
        }

        if (next > nPos) {
            //
            // Searched last interval, done.
            //
            break;
        }

        //
        // Next has advanced.  Check what contig it is in.
        //
        if (next < nPos) {
            nextBoundary = ContigStartPos(pos[next].t);
            //
            // Advance next to the maximum position within this contig that is
            // just after where the interval starting at cur is, or the first
            // position in the next contig.
            //
            int prevNext = next;
            AdvanceIndexToPastInterval(pos, nPos, intervalLength, contigLength, ContigStartPos,
                                       cur, curBoundary, next, nextBoundary);
            if (prevNext != next or recountInterval) {
                curSize = SumAnchors(pos, cur, next);
            }
        }

        // if next >= nPos, the boundary stays the same.

        //
        //  When searching multiple contigs, it is important to know the
        //  boundary of the contig that this anchor is in so that clusters
        //  do not span multiple contigs.  Find the (right hand side)
        //  boundary of the current contig.
        //

        curBoundary = ContigStartPos(pos[cur].t);
        contigLength  = ContigStartPos.Length(pos[cur].t);

        //
        // Previously tried to advance half.  This is being removed since
        // proper heuristics are making it not necessary to use.
        //
    }
    if (curSize > minSize) {
        start.push_back(maxStart);
        end.push_back(maxEnd);
    }
}

template<typename T_MatchList,
         typename T_PValueFunction, 
         typename T_WeightFunction,
         typename T_SequenceBoundaryDB,
         typename T_ReferenceSequence,
         typename T_Sequence>
int FindMaxIncreasingInterval(
    // Input
    // readDir is used to indicate if the interval that is being stored is 
    // in the forward or reverse strand.  This is important later when 
    // refining alignments so that the correct sequence is aligned back 
    // to the reference.
    int readDir, 
    T_MatchList &pos, 
    // How many values to search through for a max set.
    DNALength intervalLength,  
    // How many sets to keep track of
    VectorIndex nBest, 
    // End search for intervals at boundary positions
    // stored in seqBoundaries
    T_SequenceBoundaryDB & ContigStartPos,
    // First rand intervals by their p-value
    T_PValueFunction &MatchPValueFunction,  
    // When ranking intervals, sum over weights determined by 
    // MatchWeightFunction
    T_WeightFunction &MatchWeightFunction,  
    // Output.
    // The increasing interval coordinates, 
    // in order by queue weight.
    WeightedIntervalSet &intervalQueue, 
    T_ReferenceSequence &reference, 
    T_Sequence &query,
    IntervalSearchParameters &params,
    vector<BasicEndpoint<ChainedMatchPos> > *chainEndpointBuffer,
    ClusterList &clusterList,
    VarianceAccumulator<float> &accumPValue, 
    VarianceAccumulator<float> &accumWeight,
    VarianceAccumulator<float> &accumNumAnchorBases) {
    (void)(accumNumAnchorBases);

    int maxLISSize = 0;
    if (params.fastMaxInterval) {
        maxLISSize = FastFindMaxIncreasingInterval(
                readDir, pos, intervalLength,
                nBest, ContigStartPos, 
                MatchPValueFunction, MatchWeightFunction,
                intervalQueue, reference, query,
                params, chainEndpointBuffer, clusterList,
                accumPValue, accumWeight);
    } else {
        maxLISSize = ExhaustiveFindMaxIncreasingInterval(
                readDir, pos, intervalLength,
                nBest, ContigStartPos, 
                MatchPValueFunction, MatchWeightFunction,
                intervalQueue, reference, query,
                params, chainEndpointBuffer, clusterList,
                accumPValue, accumWeight);
    }

    if (params.aggressiveIntervalCut and intervalQueue.size() >= 3) {
        // aggressiveIntervalCut mode: 
        // only pick up the most promising intervals if we can classify
        // intervals into 'promising' and 'non-promising' clusters.
        WeightedIntervalSet::iterator it = intervalQueue.begin();
        int sz = intervalQueue.size();
        vector<float> pValues, ddPValues;
        pValues.resize(sz);
        ddPValues.resize(sz);
        float sumPValue = 0;
        int i = 0;
        for(; it != intervalQueue.end(); i++,it++) {
            sumPValue += (*it).pValue;
            pValues[i] = (*it).pValue;
        }
        //We will attemp to divide intervals into two clusters, promising 
        //and non-promising.
        float prevSumPValue = pValues[0];
        // ddPValue[i], where i in [1...n-2] is the difference of 
        //    (1) (mean pvalue of [0...i-1] minus pvalue[i])
        // and 
        //    (2) (pvalue[i] minus mean pvalue of [i+1...n-1])
        // pValues are all negative, the lower the better.
        // if ddPValue is negative, interval i is closer to cluster [i+1...n-1],
        // otherwise, interval i is closer to cluster [0..i-1].
        it = intervalQueue.begin();
        it ++; //both it and i should point to the second interval in intervalQueue.
        for(i = 1; i < sz - 1; i++,it++) {
            ddPValues[i] = (prevSumPValue / i) +
                           (sumPValue - prevSumPValue - pValues[i]) / (sz - i - 1) -
                           2 * pValues[i];
            if (ddPValues[i] <= params.ddPValueThreshold) {
                // PValue of interval i is much closer to cluster [i+1...n-1] than
                // to cluster [0...i-1]. Mean pValue of cluster [0..i-1] 
                // minus mean pvalue of cluster [i+1...n-1] < 2 * -500
                break;
            }
            prevSumPValue += pValues[i];
        }
        if (it != intervalQueue.end()) {
            // Erase intervals in the non-promising cluster.
            intervalQueue.erase(it, intervalQueue.end());
            // Recompute accumPValue, accmWeight, clusterList;
            accumPValue.Reset();
            accumWeight.Reset();
            clusterList.Clear();
            for(it = intervalQueue.begin(); it != intervalQueue.end(); it++) {
                accumPValue.Append((*it).pValue);
                accumWeight.Append((*it).size);
                clusterList.Store((*it).totalAnchorSize, (*it).start, (*it).end, (*it).nAnchors);
            }
        }
    }
    return maxLISSize;
}

template<typename T_MatchList,
         typename T_PValueFunction, 
         typename T_WeightFunction,
         typename T_SequenceBoundaryDB,
         typename T_ReferenceSequence,
         typename T_Sequence>
int FastFindMaxIncreasingInterval(
        // Input
        // readDir is used to indicate if the interval that is being stored is in the forward
        // or reverse strand.  This is important later when refining alignments so that the
        // correct sequene is aligned back to the reference.
        int readDir, 
        T_MatchList &pos, 
        // How many values to search through for a max set.
        DNALength intervalLength,  
        // How many sets to keep track of
        VectorIndex nBest, 
        // End search for intervals at boundary positions
        // stored in seqBoundaries
        T_SequenceBoundaryDB & ContigStartPos,
        // First rand intervals by their p-value
        T_PValueFunction &MatchPValueFunction,  
        // When ranking intervals, sum over weights determined by MatchWeightFunction
        T_WeightFunction &MatchWeightFunction,  
        // Output.
        // The increasing interval coordinates, 
        // in order by queue weight.
        WeightedIntervalSet &intervalQueue, 
        T_ReferenceSequence &reference, 
        T_Sequence &query,
        IntervalSearchParameters &params,
        vector<BasicEndpoint<ChainedMatchPos> > *chainEndpointBuffer,
        ClusterList &clusterList,
        VarianceAccumulator<float> &accumPValue, 
        VarianceAccumulator<float> &accumWeight) {
    (void)(nBest); (void)(query); (void)(reference);
 
    WeightedIntervalSet sdpiq;
    VectorIndex cur = 0;
    vector<VectorIndex> lisIndices;
    //
    // Initialize the first interval.
    //
    if (pos.size() == 0) {
        return 0;
    }

    int lisSize;
    float lisPValue;
    T_MatchList lis;
    int noOvpLisSize = 0;
    int noOvpLisNBases = 0;

    //
    // Search for clusters of intervals within the pos array within
    // pos[cur...next).  The value of 'next' should be the first anchor
    // outside the possible range to cluster, or the end of the anchor list.

    VectorIndex next = cur + 1;
    DNALength curBoundary = 0, nextBoundary = 0;
    DNALength contigLength = ContigStartPos.Length(pos[cur].t);
    DNALength endOfCurrentInterval = curBoundary + contigLength;
    (void)(curBoundary); (void)(nextBoundary); (void)(endOfCurrentInterval);

    vector<UInt> scores, prevOpt;
    vector<DNALength> start, end;

    StoreLargestIntervals(pos, ContigStartPos, intervalLength, 30, start, end);

    VectorIndex i;
    VectorIndex posi;
    int maxLISSize = 0;
    
    for (posi = 0; posi < start.size(); posi++) {
        lis.clear();
        lisIndices.clear();
        cur = start[posi];
        next = end[posi];
        if (next - cur == 1) {
            //
            // Just one match in this interval, don't invoke call to global chain since it is given.
            //
            lisSize = 0;
            lisIndices.push_back(0);
        }
        else {
            //
            // Find the largest set of increasing intervals that do not overlap.
            //
            if (params.globalChainType == 0) {
                lisSize = GlobalChain<ChainedMatchPos, BasicEndpoint<ChainedMatchPos> >(pos, cur, next, 
                        lisIndices, chainEndpointBuffer);
            }
            else {
                //
                //  A different call that allows for indel penalties.
                //
                lisSize = RestrictedGlobalChain(&pos[cur],next - cur, 0.1, lisIndices, scores, prevOpt);
            }
        }

        // Maybe this should become a function?
        for (i = 0; i < lisIndices.size(); i++) {
            lis.push_back(pos[lisIndices[i]+cur]);
        }

        // 
        // Compute pvalue of this match.
        //
        if (lis.size() > 0) {
            lisPValue = MatchPValueFunction.ComputePValue(lis, noOvpLisNBases, noOvpLisSize);
        }
        else {
            lisPValue = 0;
        }

        if (lisSize > maxLISSize) {
            maxLISSize  = lisSize;
        }

        //
        // Insert the interval into the interval queue maintaining only the 
        // top 'nBest' intervals. 
        //
        WeightedIntervalSet::iterator lastIt = intervalQueue.begin();
        MatchWeight lisWeight = MatchWeightFunction(lis);
        VectorIndex lisEnd = lis.size() - 1;

        accumPValue.Append(lisPValue);
        accumWeight.Append(lisWeight);

        if (lisPValue < params.maxPValue and lisSize > 0) {
            WeightedInterval weightedInterval(lisWeight, noOvpLisSize, noOvpLisNBases, 
                    lis[0].t, lis[lisEnd].t + lis[lisEnd].GetLength(), 
                    readDir, lisPValue, 
                    lis[0].q, lis[lisEnd].q + lis[lisEnd].GetLength(), 
                    lis);
            intervalQueue.insert(weightedInterval);
            if (weightedInterval.isOverlapping == false) {
                clusterList.Store((float)noOvpLisNBases, lis[0].t, lis[lis.size()-1].t, noOvpLisSize);
            }
            if (params.verbosity > 1) {
                cout << "Weighted Interval to insert:"<< endl << weightedInterval << endl;
                cout << "Interval Queue:"<< endl << intervalQueue << endl;
            }
        }
    }
    return maxLISSize;
}

template<typename T_MatchList,
         typename T_PValueFunction, 
         typename T_WeightFunction,
         typename T_SequenceBoundaryDB,
         typename T_ReferenceSequence,
         typename T_Sequence>
int ExhaustiveFindMaxIncreasingInterval(
        // Input
        // readDir is used to indicate if the interval that is being stored is in the forward
        // or reverse strand.  This is important later when refining alignments so that the
        // correct sequene is aligned back to the reference.
        int readDir, 
        T_MatchList &pos, 
        // How many values to search through for a max set.
        DNALength intervalLength,  
        // How many sets to keep track of
        VectorIndex nBest, 
        // End search for intervals at boundary positions
        // stored in seqBoundaries
        T_SequenceBoundaryDB & ContigStartPos,
        // First rand intervals by their p-value
        T_PValueFunction &MatchPValueFunction,  
        // When ranking intervals, sum over weights determined by MatchWeightFunction
        T_WeightFunction &MatchWeightFunction,  
        // Output.
        // The increasing interval coordinates, 
        // in order by queue weight.
        WeightedIntervalSet &intervalQueue, 
        T_ReferenceSequence &reference, 
        T_Sequence &query,
        IntervalSearchParameters &params,
        vector<BasicEndpoint<ChainedMatchPos> > *chainEndpointBuffer,
        ClusterList &clusterList,
        VarianceAccumulator<float> &accumPValue, 
        VarianceAccumulator<float> &accumWeight) {
    (void)(nBest); (void)(query); (void)(reference);

    WeightedIntervalSet sdpiq;
    VectorIndex cur = 0;
    VectorIndex nPos = pos.size();
    //
    // Initialize the first interval.
    //
    if (pos.size() == 0) {
        return 0;
    }

    int lisSize;
    float lisPValue;
    T_MatchList lis;
    int noOvpLisSize = 0;
    int noOvpLisNBases = 0;

    //
    // Search for clusters of intervals within the pos array within
    // pos[cur...next).  The value of 'next' should be the first anchor
    // outside the possible range to cluster, or the end of the anchor list.

    VectorIndex next = cur + 1;
    DNALength curBoundary = 0, nextBoundary = 0;
    DNALength contigLength = ContigStartPos.Length(pos[cur].t);
    DNALength endOfCurrentInterval = curBoundary + contigLength;
    (void)(curBoundary); (void)(nextBoundary); (void)(endOfCurrentInterval);

    curBoundary = ContigStartPos(pos[cur].t);
    nextBoundary = ContigStartPos(pos[next].t);  
    vector<UInt> scores, prevOpt;

    //
    // Advance next until the anchor is outside the interval that
    // statrts at 'cur', and is inside the same contig that the anchor
    // at cur is in.
    //

    DNALength curIntervalLength = NumRemainingBases(pos[cur].q, intervalLength);
    (void)(curIntervalLength);

    AdvanceIndexToPastInterval(pos, nPos, intervalLength, contigLength, ContigStartPos,
    cur, curBoundary, next, nextBoundary);


    vector<VectorIndex> lisIndices;
    VectorIndex i;

    //
    // Do some preprocessing.  If the number of anchors considered for this hit is 1, 
    // the global chain is this sole ancor.  Don't bother calling GlobalChain
    // since it allocates and deallocates extra memory.
    //

    int maxLISSize = 0;

    //
    // Search intervals until cur reaches the end of the list of
    // anchors.
    //
    while ( cur < nPos ) {
        //
        // Search the local interval for a LIS larger than a previous LIS.
        //
        lis.clear();
        lisIndices.clear();

        if (next - cur == 1) {
            //
            // Just one match in this interval, don't invoke call to global chain since it is given.
            //
            lisSize = 1;
            lisIndices.push_back(0);
        }
        else {
            //
            // Find the largest set of increasing intervals that do not overlap.
            //
            if (params.globalChainType == 0) {
                lisSize = GlobalChain<ChainedMatchPos, BasicEndpoint<ChainedMatchPos> >(pos, cur, next, 
                        lisIndices, chainEndpointBuffer);
            }
            else {
                //
                //  A different call that allows for indel penalties.
                //
                lisSize = RestrictedGlobalChain(&pos[cur],next - cur, 0.1, lisIndices, scores, prevOpt);
            }
        }

        // Maybe this should become a function?
        for (i = 0; i < lisIndices.size(); i++) {    lis.push_back(pos[lisIndices[i]+cur]); }


        // 
        // Compute pvalue of this match.
        //
        lisPValue = MatchPValueFunction.ComputePValue(lis, noOvpLisNBases, noOvpLisSize);

        if (lisSize > maxLISSize) {
            maxLISSize  = lisSize;
        }

        //
        // Insert the interval into the interval queue maintaining only the 
        // top 'nBest' intervals. 
        //

        WeightedIntervalSet::iterator lastIt = intervalQueue.begin();
        MatchWeight lisWeight = MatchWeightFunction(lis);
        VectorIndex lisEnd = lis.size() - 1;

        accumPValue.Append(lisPValue);
        accumWeight.Append(lisWeight);

        if (lisPValue < params.maxPValue and lisSize > 0) {
            WeightedInterval weightedInterval(lisWeight, noOvpLisSize, noOvpLisNBases, 
                    lis[0].t, lis[lisEnd].t + lis[lisEnd].GetLength(), 
                    readDir, lisPValue, 
                    lis[0].q, lis[lisEnd].q + lis[lisEnd].GetLength(), 
                    lis);
            intervalQueue.insert(weightedInterval);
            if (weightedInterval.isOverlapping == false) {
                clusterList.Store((float)noOvpLisNBases, lis[0].t, lis[lis.size()-1].t, noOvpLisSize);
            }
            if (params.verbosity > 1) {
                cout << "Weighted Interval to insert:"<< endl << weightedInterval << endl;
                cout << "Interval Queue:"<< endl << intervalQueue << endl;
            }
        }

        //
        // Done scoring current interval.  At this point the range
        // pos[cur...next) has been searched for a max increasing
        // interval.  Find a new range that will possibly yield a new
        // maximum interval.  
        // There are a few cases to consider:
        //
        //
        //genome  |---+----+------------+------+-----------------------|
        //  anchors  cur  cur+1        next   next+1
        //
        // Case 1.  The range on the target pos[ cur+1 ... next].t is a
        // valid interval (it is roughly the length of the read).  In this
        // case increase cur and next by 1, and search this range.
        //
        // genome  |---+----+------------+------+-----------------------|
        //            cur  cur+1        next   next+1
        // read interval   ====================
        //
        // Case 2.  The range on the target pos[cur+1 ... next] is not a
        // valid interval, and it is much longer than the length of the
        // read.  This implies that it is impossible to increase the score
        // of the read by including both 
        //
        // genome  |---+----+--------------------------------+-----+---|
        //            cur  cur+1                             next next+1
        // read interval   ==================== 
        //
        // Advance the interval until it includes the next anchor
        //
        // genome  |---+----+------------------+-------------+-----+---|
        //            cur  cur+1             cur+n          next next+1
        // read interval                     ==================== 
        // 

        // First advance pointer in anchor list.  If this advances to the
        // end, done and no need for further logic checking (break now).


        // 
        // If the next position is not within the same contig as the current,
        // advance the current to the next since it is impossible to find
        // any more intervals in the current pos.
        //
        if (curBoundary != nextBoundary) {
            cur = next;
            curBoundary = nextBoundary;

            //
            // Start the search for the first interval in the next contig
            // just after the current position.
            //
            if (next < nPos) {
                next = cur + 1;
            }
        }
        else {
            cur++;
            if (cur >= nPos)
                break;

            //
            //  Look for need to advance the interval within the same
            //  contig.  If the start of the next match is well past the
            //  length of this read, keep moving forward matches until it is
            //  possible to include the next interval in the matches for
            //  this read.
            // 
            if (params.warp) {
                while (cur < next and 
                        next < nPos and 
                        pos[next].t - pos[cur].t >= intervalLength) {
                    //
                    // It is impossible to increase the max interval weight any more
                    // using pos[cur] when pos[next] is too far away from
                    // pos[cur].  This is because the same set of anchors are used
                    // when clustering pos[cur+1] ... pos[next] as
                    // pos[cur] ... pos[next].  
                    // Advance cur until pos[cur] is close enough to matter again.
                    cur++;
                }
            }
            //
            // Advance the next to outside this interval.
            next++;
        }

        if (next > nPos) {
            //
            // Searched last interval, done.
            //
            break;
        }

        //
        // Next has advanced.  Check what contig it is in.
        //
        if (next < nPos) {
            nextBoundary = ContigStartPos(pos[next].t);
            //
            // Advance next to the maximum position within this contig that is
            // just after where the interval starting at cur is, or the first
            // position in the next contig.
            //
            AdvanceIndexToPastInterval(pos, nPos, intervalLength, contigLength, ContigStartPos,
                    cur, curBoundary, next, nextBoundary);
        }
        // if next >= nPos, the boundary stays the same.
        //
        //  When searching multiple contigs, it is important to know the
        //  boundary of the contig that this anchor is in so that clusters
        //  do not span multiple contigs.  Find the (right hand side)
        //  boundary of the current contig.
        //

        curBoundary = ContigStartPos(pos[cur].t);
        contigLength  = ContigStartPos.Length(pos[cur].t);

        //
        // Previously tried to advance half.  This is being removed since
        // proper heuristics are making it not necessary to use.
        //
    }

    return  maxLISSize;
}

#endif
