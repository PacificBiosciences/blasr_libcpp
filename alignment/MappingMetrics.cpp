#include <iomanip>
#include <unistd.h>
#include "MappingMetrics.hpp"
//In order to use clock_gettime in LINUX, add -lrt 

#ifdef __APPLE__
#pragma weak clock_gettime
int clock_gettime(clockid_t clk_id, struct timespec *tp) {
    kern_return_t   ret;
    clock_serv_t    clk;
    clock_id_t clk_serv_id;
    mach_timespec_t tm;

    uint64_t start, end, delta, nano;

    task_basic_info_data_t tinfo;
    task_thread_times_info_data_t ttinfo;
    mach_msg_type_number_t tflag;

    int retval = -1;
    switch (clk_id) {
        case CLOCK_REALTIME:
        case CLOCK_MONOTONIC:
            clk_serv_id = clk_id == CLOCK_REALTIME ? CALENDAR_CLOCK : SYSTEM_CLOCK;
            if (KERN_SUCCESS == (ret = host_get_clock_service(mach_host_self(), clk_serv_id, &clk))) {
                if (KERN_SUCCESS == (ret = clock_get_time(clk, &tm))) {
                    tp->tv_sec  = tm.tv_sec;
                    tp->tv_nsec = tm.tv_nsec;
                    retval = 0;
                }
            }
            if (KERN_SUCCESS != ret) {
                errno = EINVAL;
                retval = -1;
            }
            break;
        case CLOCK_PROCESS_CPUTIME_ID:
        case CLOCK_THREAD_CPUTIME_ID:
            start = mach_absolute_time();
            if (clk_id == CLOCK_PROCESS_CPUTIME_ID) {
                getpid();
            } else {
                sched_yield();
            }
            end = mach_absolute_time();
            delta = end - start;
            if (0 == __clock_gettime_inf.denom) {
                mach_timebase_info(&__clock_gettime_inf);
            }
            nano = delta * __clock_gettime_inf.numer / __clock_gettime_inf.denom;
            tp->tv_sec = nano * 1e-9;  
            tp->tv_nsec = nano - (tp->tv_sec * 1e9);
            retval = 0;
            break;
        default:
            errno = EINVAL;
            retval = -1;
    }
    return retval;
}
#endif // __APPLE__

Timer::Timer(std::string _header) {
    keepHistogram = false;
    keepList      = false;
    totalElapsedClock = 0;
    header        = _header;
    elapsedClockMsec = 0;
    elapsedTime   = 0.0;
}

int Timer::ListSize() {
    return msecList.size();
}

void Timer::PrintHeader(std::ostream &out) {
    if (msecList.size() > 0) {
        out << header << " ";
    }
}

void Timer::PrintListValue(std::ostream &out, int index) {
    if (msecList.size() > 0) {
        out << msecList[index] << " ";
    }
}
void Timer::Tick() {
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &cpuclock[0]);
}

void Timer::SetStoreElapsedTime(bool value) {
    keepList = value;
}

void Timer::SetStoreHistgram(bool value) {
    keepHistogram = value;
}

void Timer::Tock() {
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &cpuclock[1]);
    elapsedClockMsec   = (cpuclock[1].tv_nsec - cpuclock[0].tv_nsec)/1000;
    totalElapsedClock += elapsedClockMsec;
    elapsedTime  = ((1.0)*elapsedClockMsec);
    if (keepHistogram) {
        // keep a histogram in number of milliseconds per operation
        if (histogram.find(elapsedClockMsec) == histogram.end()) {
            histogram[elapsedClockMsec] = 1;
        }
        else {
            histogram[elapsedClockMsec]++;
        }
    }
    if (keepList) {
        msecList.push_back(elapsedClockMsec);
    }
}

void Timer::Add(const Timer &rhs) {
    elapsedClockMsec += rhs.elapsedClockMsec;
    elapsedTime  += rhs.elapsedTime;
    totalElapsedClock += rhs.totalElapsedClock;
    msecList.insert(msecList.end(), rhs.msecList.begin(), rhs.msecList.end());
}

void Timer::SetHeader(std::string _header) {
    header = _header;
}

void MappingClocks::AddCells(int nCells) {
    nCellsPerSample.push_back(nCells);
}

void MappingClocks::AddBases(int nBases) {
    nBasesPerSample.push_back(nBases);
}

int MappingClocks::GetSize() {
    return total.ListSize();
}

MappingClocks::MappingClocks() {
    total.SetHeader("Total");
    findAnchors.SetHeader("FindAnchors");
    mapToGenome.SetHeader("MapToGenome");
    sortMatchPosList.SetHeader("SortMatchPosList");
    findMaxIncreasingInterval.SetHeader("FindMaxIncreasingInterval");
    alignIntervals.SetHeader("AlignIntervals");
}

void MappingClocks::PrintHeader(std::ostream &out) {
    total.PrintHeader(out);
    findAnchors.PrintHeader(out);
    mapToGenome.PrintHeader(out);
    sortMatchPosList.PrintHeader(out);
    findMaxIncreasingInterval.PrintHeader(out);
    alignIntervals.PrintHeader(out); 
}

void MappingClocks::PrintList(std::ostream &out, int index) {

    total.PrintListValue(out,index);
    findAnchors.PrintListValue(out,index);
    mapToGenome.PrintListValue(out,index);
    sortMatchPosList.PrintListValue(out,index);
    findMaxIncreasingInterval.PrintListValue(out,index);
    alignIntervals.PrintListValue(out,index);
    if (nCellsPerSample.size() > 0) {
        out << nCellsPerSample[index] << " ";
    }
    if (nBasesPerSample.size() > 0) {
        out << nBasesPerSample[index] << " ";
    }
    out << std::endl;
}

void MappingClocks::SetStoreList(bool value) {
    total.SetStoreElapsedTime(value);
    findAnchors.SetStoreElapsedTime(value);
    mapToGenome.SetStoreElapsedTime(value);
    sortMatchPosList.SetStoreElapsedTime(value);
    findMaxIncreasingInterval.SetStoreElapsedTime(value);
    alignIntervals.SetStoreElapsedTime(value);
}

void MappingClocks::AddClockTime(const MappingClocks &rhs) {
    total.Add(rhs.total);
    findAnchors.Add(rhs.findAnchors);
    mapToGenome.Add(rhs.mapToGenome);
    sortMatchPosList.Add(rhs.sortMatchPosList);
    findMaxIncreasingInterval.Add(rhs.findMaxIncreasingInterval);
    alignIntervals.Add(rhs.alignIntervals);
}

MappingMetrics::MappingMetrics() {
    numReads = 0;
    numMappedReads = 0;
    numMappedBases = 0;
    anchorsPerRead = 0;
    totalAnchorsForMappedReads = 0;
    totalAnchors = 0;
}
void MappingMetrics::StoreSDPPoint(int nBases, int nSDPAnchors, int nClock) {
    sdpBases.push_back(nBases);
    sdpAnchors.push_back(nSDPAnchors);
    sdpClock.push_back(nClock);
}

void MappingMetrics::SetStoreList(bool value) {
    clocks.SetStoreList(value);
}

void MappingMetrics::PrintSeconds(std::ostream& out, long sec ){
    out << sec << " Msec";
}

void MappingMetrics::PrintFraction(std::ostream &out, float frac) {
    out << std::setprecision(2) << frac;
}

void MappingMetrics::RecordNumAlignedBases(int nBases) {
    mappedBases.push_back(nBases);
}

void MappingMetrics::RecordNumCells(int nCells) {
    cellsPerAlignment.push_back(nCells);
}

void MappingMetrics::Collect(MappingMetrics &rhs) {
    clocks.AddClockTime(rhs.clocks);
    totalAnchors += rhs.totalAnchors;
    numReads += rhs.numReads;
    numMappedReads += rhs.numMappedReads;
    totalAnchorsForMappedReads += rhs.totalAnchorsForMappedReads;
    mappedBases.insert(mappedBases.end(), rhs.mappedBases.begin(), rhs.mappedBases.end());
    cellsPerAlignment.insert(cellsPerAlignment.end(), rhs.cellsPerAlignment.begin(), rhs.cellsPerAlignment.end());
}

void MappingMetrics::CollectSDPMetrics(MappingMetrics &rhs) {
    sdpAnchors.insert(sdpAnchors.end(), rhs.sdpAnchors.begin(), rhs.sdpAnchors.end());
    sdpBases.insert(sdpBases.end(), rhs.sdpBases.begin(), rhs.sdpBases.end());
    sdpClock.insert(sdpClock.end(), rhs.sdpClock.begin(), rhs.sdpClock.end());
}

void MappingMetrics::PrintSDPMetrics(std::ostream &out) {
    out << "nbases ncells time" << std::endl;
    for (size_t i = 0; i < sdpAnchors.size(); i++) {
        out << sdpBases[i] << " " << sdpAnchors[i] << " " << sdpClock[i] << std::endl;
    }
}

void MappingMetrics::PrintFullList(std::ostream &out) {
    // 
    // Print the full header
    //
    clocks.PrintHeader(out);
    out << " MappedBases Cells " << std::endl;
    //
    // Print all values: clocks + bases and cells.
    //
    int i;
    for (i = 0; i < clocks.GetSize(); i++) {
        clocks.PrintList(out,i);
        //      out << mappedBases[i] << " " << cellsPerAlignment[i] << endl;
    }
}

void MappingMetrics::PrintSummary(std::ostream &out) {
    out << "Examined " << numReads << std::endl;
    out << "Mapped   " << numMappedReads << std::endl;
    out << "Total mapping time\t";
    PrintSeconds(out, clocks.total.elapsedClockMsec);
    out << " \t";
    PrintSeconds(out, (1.0*clocks.total.elapsedClockMsec)/numReads);
    out << " /read" << std::endl;
    out << "      find anchors\t";
    PrintSeconds(out, clocks.mapToGenome.elapsedClockMsec);
    out << " \t";
    PrintSeconds(out, (1.0*clocks.mapToGenome.elapsedClockMsec)/numReads);
    out << std::endl;
    out << "      sort anchors\t";
    PrintSeconds(out, clocks.sortMatchPosList.elapsedClockMsec);
    out << " \t";
    PrintSeconds(out, (1.0*clocks.sortMatchPosList.elapsedClockMsec)/numReads);
    out << std::endl;
    out << " find max interval\t";
    PrintSeconds(out, clocks.findMaxIncreasingInterval.elapsedClockMsec);
    out << " \t";
    PrintSeconds(out, (1.0*clocks.findMaxIncreasingInterval.elapsedClockMsec)/numReads);
    out << std::endl;
    out << "Total anchors: " << totalAnchors << std::endl;
    out << "   Anchors per read: " << (1.0*totalAnchors) / numReads << std::endl;
    out << "Total mapped: " << totalAnchorsForMappedReads << std::endl;
    out << "   Anchors per mapped read: " << (1.0*totalAnchorsForMappedReads) / numMappedReads << std::endl;
}

void MappingMetrics::AddClock(MappingClocks &clocks) {
    clocks.AddClockTime(clocks);
}
