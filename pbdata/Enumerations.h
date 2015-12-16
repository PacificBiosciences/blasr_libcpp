// Author: Mark Chaisson

#ifndef _BLASR_ENUMERATIONS_HPP_
#define _BLASR_ENUMERATIONS_HPP_

enum class FileType {Fasta,
    Fastq,
    HDFPulse,
    Fourbit,
    HDFBase,
    HDFCCS,
    HDFCCSONLY,
    PBBAM,
    PBDATASET,
    None
};

typedef enum T_Strand {Forward, Reverse} Strand;

typedef enum T_PlatformType {
    Astro=1,
    Springfield=2,
    Sequel=4,
    NoPlatform
} PlatformId;

typedef enum T_RegionType {
    Adapter,
    Insert, 
    HQRegion,
    BarCode,
    UnknownRegionType
} RegionType;

typedef enum T_PulseMetricType {
    QualityValueMetric,ClassifierQVMetric,StartTimeMetric,PulseWidthMetric,WidthInFramesMetric,
    pkmidMetric,IPDMetric,LightMetric, PreBaseFramesMetric
} PulseMetricType;


typedef enum T_AlignMode { NoAlignMode, Fullread, Subread, CCSDeNovo, CCSFullPass, CCSAllPass, ZmwSubreads } AlignMode;

#endif
