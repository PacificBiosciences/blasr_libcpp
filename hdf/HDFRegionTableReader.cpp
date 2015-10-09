#include <cassert>
#include "HDFRegionTableReader.hpp"

using namespace std;

int HDFRegionTableReader::Initialize(string &regionTableFileName, 
        const H5::FileAccPropList & fileAccPropList) {
    /*
     * Initialize access to the HDF file.
     */
    try {
        regionTableFile.Open(regionTableFileName.c_str(), H5F_ACC_RDONLY, fileAccPropList);
    }
    catch (H5::Exception &e) {
        cout << e.getDetailMsg() << endl;
        return 0;
    }
    if (pulseDataGroup.Initialize(regionTableFile.rootGroup, "PulseData") == 0) {
        return 0;
    }
    if (pulseDataGroup.ContainsObject("Regions") == 0) {
        fileContainsRegionTable = false;
        return 0;
    }
    else {
        fileContainsRegionTable = true;
    }

    if (regions.Initialize(pulseDataGroup, "Regions") == 0) {
        return 0;
    }

    if (columnNames.Initialize(regions, "ColumnNames") == 0) {
        return 0;
    }
    if (regionTypes.Initialize(regions, "RegionTypes") == 0) {
        return 0;
    }
    if (regionDescriptions.Initialize(regions, "RegionDescriptions") == 0) {
        return 0;
    }
    if (regionSources.Initialize(regions,  "RegionSources") == 0) {
        return 0;
    }

    nRows = regions.GetNRows();
    isInitialized_ = true;
    curRow = 0;
    return 1;
}

bool HDFRegionTableReader::IsInitialized(void) const {
    return isInitialized_;
}

bool HDFRegionTableReader::HasRegionTable(void) const {
    assert(IsInitialized() or false == "HDFRegionTable is not initialize!");
    return fileContainsRegionTable;
}

int HDFRegionTableReader::GetNext(RegionAnnotation &annotation) {
    assert(IsInitialized() or false == "HDFRegionTable is not initialize!");
    //
    // Bail with no-op if this is the last row.
    //
    if (fileContainsRegionTable == false) {
        return 0;
    }

    if (curRow == nRows) {
        return 0;
    }

    regions.Read(curRow, curRow+1, annotation.row);
    ++curRow;
    return 1;
}	


void HDFRegionTableReader::Close() {
    isInitialized_ = false;
    fileContainsRegionTable = false;
    columnNames.Close();
    regionTypes.Close();
    regionDescriptions.Close();
    regionSources.Close();
    pulseDataGroup.Close();
    regions.Close();
    regionTableFile.Close();
}

// Note that (1) there is NO GUARANTEE that region annotations in hdf5 
// `Regions` dataset be sorted in any order, so we cannot iterate over
// `Regions` in order to traverse zmws in order.
// (2) region table of a million zmws is approximately 5M.
void HDFRegionTableReader::ReadTable(RegionTable & table) {
    assert(IsInitialized() or false == "HDFRegionTable is not initialize!");
    table.Reset();

    if (fileContainsRegionTable) {
        // Read attributes.
        std::vector<std::string> names, types, descs, sources;
        if (columnNames.IsInitialized()) columnNames.Read(names);
        if (regionTypes.IsInitialized()) regionTypes.Read(types);
        else {
            cout << "ERROR MUST HAVE REGIONTYPES" << endl;
            exit(1);
        }
        if (regionDescriptions.IsInitialized()) regionDescriptions.Read(descs);
        if (regionSources.IsInitialized()) regionSources.Read(sources);

        // Read region annotations
        std::vector<RegionAnnotation> ras;
        ras.resize(nRows);
        assert(curRow == 0);
        for (; curRow < nRows; curRow++) {
            regions.Read(curRow, curRow+1, ras[curRow].row);
        }

        // Reconstruct table
        table.ConstructTable(ras, types);
        table.ColumnNames(names);
        table.RegionDescriptions(descs);
        table.RegionSources(sources);
    }
}

void HDFRegionTableReader::GetMinMaxHoleNumber(UInt &minHole,
                                               UInt &maxHole) {
    assert(IsInitialized() or false == "HDFRegionTable is not initialize!");
    // Hole numbers may not be sorted ascendingly, so do not
    // return the first and last hole numbers as the min and max.
    UInt saveCurRow = curRow;
    curRow = 0;
    bool init = false;
    RegionAnnotation annotation;
    while (GetNext(annotation) == 1) {
        UInt curHole = annotation.GetHoleNumber();
        if (not init) {
            minHole = maxHole = curHole;
            init = true;
        } else {
            minHole = (minHole > curHole)?(curHole):(minHole);
            maxHole = (maxHole < curHole)?(curHole):(maxHole);
        }
    }
    curRow = saveCurRow;
}
