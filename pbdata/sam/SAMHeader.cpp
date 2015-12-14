#include "SAMHeader.hpp"
#include <algorithm>

void SAMHeader::StoreValues(std::vector<SAMKeywordValuePair> &kvPairs,
                            int lineNumber) {
  for (size_t i = 0; i < kvPairs.size(); i++) {
    if (kvPairs[i].key == "VN") {
      formatVersion = kvPairs[i].value;
    }
    else if (kvPairs[i].key == "SO") {
      string value = kvPairs[i].value;
      std::transform(value.begin(), value.end(), value.begin(), ::tolower);
      if (value == "unknown" || value == "unsorted") {
        sortingOrder = unknown;
      }
      else if (value == "sorted") {
        sortingOrder = sorted;
      }
      else if (value == "queryname") {
        sortingOrder =queryname;
      }
      else if (value == "coordinate") {
        sortingOrder = coordinate;
      }
      else {
        std::cout << "Invalid sorting order " << kvPairs[i].value << " at line " << lineNumber;
      }
    }
  }
}
