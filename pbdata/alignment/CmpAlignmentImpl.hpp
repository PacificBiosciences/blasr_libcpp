
template<typename T_Field>
void CmpAlignment::StoreField(std::string fieldName, T_Field* fieldValues, DSLength length) {
    fields[fieldName].resize(length);
    memcpy(&fields[fieldName][0], fieldValues, length * sizeof(T_Field));
}    
