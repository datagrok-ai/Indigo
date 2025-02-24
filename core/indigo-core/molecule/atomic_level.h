#include "molecule.h"
#include "ket_document.h"

using namespace indigo;

const float PI = 3.141592f;
const float STD_BOND = 1.54f;

const int V3K_COUNTS_SHIFT = 14;
const int V3K_IDX_SHIFT= 7;
const std::string V3K_HEADER_FIRST_LINE = "\n Indigo \n\n";
const std::string V3K_HEADER_SECOND_LINE = "  0  0  0  0  0  0            999 V3000\n";
const std::string V3K_BEGIN_CTAB_BLOCK = "M  V30 BEGIN CTAB\n";
const std::string V3K_END_CTAB_BLOCK = "M  V30 END CTAB\n";
const std::string V3K_BEGIN_COUNTS_LINE = "M  V30 COUNTS ";
const std::string V3K_COUNTS_LINE_ENDING = " 0 0 0\n";
const std::string V3K_BEGIN_ATOM_BLOCK = "M  V30 BEGIN ATOM\n";
const std::string V3K_END_ATOM_BLOCK = "M  V30 END ATOM\n";
const std::string V3K_BEGIN_BOND_BLOCK = "M  V30 BEGIN BOND\n";
const std::string V3K_END_BOND_BLOCK = "M  V30 END BOND\n";
const std::string V3K_BOND_CONFIG = " CFG=";
const std::string V3K_BEGIN_DATA_LINE = "M  V30 ";
const std::string V3K_END = "M  END";

std::string  toAtomicLevel(const KetDocument& ket, MonomerTemplateLibrary& library);