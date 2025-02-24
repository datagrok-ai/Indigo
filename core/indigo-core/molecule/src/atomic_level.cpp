#include "molecule/atomic_level.h"
using monomers_map = std::map<std::string, std::unique_ptr<KetBaseMonomer>>;
using attach_map = std::map<std::string, KetAttachmentPoint>;

std::string strReplace(std::string str, const std::string& substr, const std::string& tostr, const bool first = false) {
    if (substr.size() > str.size())
      return str;
  
    size_t detectionDomain = str.size() - substr.size() + 1;
    detectionDomain = (substr == "\\t" || substr == "\\r" || substr == "\n") ? detectionDomain + 1 : detectionDomain;
  
    bool stop = false;
  
    if (substr.size() == 0 || detectionDomain <= 0)
      return str;
  
    std::string res;
    for (size_t i = 0; i < str.size(); ++i) {
      if (stop) {
        res += str[i];
        continue;
      }
  
      if (str[i] != substr[0] || i >= detectionDomain)
        res += str[i];
      else {
        bool equal = true;
        for (size_t j = 1; j < substr.size(); ++j) {
          if (str[i + j] != substr[j]) {
            equal = false;
            break;
          }
        }
  
        if (equal) {
          res += tostr;
          i += substr.size() - 1;
          if (first)
            stop = true;
        }
        else
          res += str[i];
      }
    }
  
    return res;
}
  
float findAngle(float dx, float dy) {
    float angle;
    if (dx == 0)
      angle = dy > 0 ? 0 : PI;
    else if (dy == 0)
      angle = dx > 0 ? -PI / 2 : PI / 2;
    else {
      float tan = dy / dx;
      float at = atan(tan);
      angle = (dx < 0) ? PI / 2 + at : -PI / 2 + at;
    }
    return angle;
}

std::string toAtomicLevel(const KetDocument& ket, MonomerTemplateLibrary& library) {
    //create all essential monomer based objects objects
    auto& monomers = ket.monomers();
    auto& connections = ket.connections();

    size_t nMonomers = monomers.size();
    size_t nConnections = connections.size();

    std::vector<std::string> templatesNumToName(0);
    std::map<std::string, size_t> templatesNameToNum;

    std::vector<size_t> monomerNumToName(nMonomers);
    std::vector<size_t> monomerNumToTemplate(nMonomers);
    std::vector<float> monomerNumToX(nMonomers);
    std::vector<float> monomerNumToY(nMonomers);
    std::vector<size_t> monomerROne(nMonomers);
    std::vector<size_t> monomerRTwo(nMonomers);

    std::vector<size_t> monomerFirstAtom(nMonomers);
    std::vector<size_t> monomerLastAtom(nMonomers);

    monomers_map::const_iterator mIt;
    size_t mI = 0;
    for (mIt = monomers.begin(); mIt != monomers.end(); ++mIt) {
        monomerNumToName[mI] = stoi(mIt->first);
        std::string temp = mIt->second->templateId();
        size_t num;
        if (templatesNameToNum.find(temp) != templatesNameToNum.end()) {
            num = templatesNameToNum[temp];
        } else {
            num = templatesNameToNum.size();
            templatesNameToNum[temp] = num;
            templatesNumToName.push_back(temp);
        }

        monomerNumToTemplate[mI] = num;
        monomerNumToX[mI] = mIt->second->position().value().x;
        monomerNumToY[mI] = mIt->second->position().value().y;

        auto& attachmentPoints = mIt->second->attachmentPoints();

        monomerROne[mI] = attachmentPoints.at("R1").attachment_atom();
        monomerRTwo[mI] = attachmentPoints.at("R2").attachment_atom();

       ++mI;
    }

    //create all essential atom based objects objects
    size_t nAtoms = 0;
    size_t nBonds = 0;

    for (size_t i = 0; i < nMonomers; ++i) {
        const std::string code = templatesNumToName[monomerNumToTemplate[i]];
        nAtoms += library.getMonomerTemplateById(code).atoms().size();
        nBonds += library.getMonomerTemplateById(code).bonds().size();
    }

    std::vector<std::string> symbolNumToSymbol(0);
    std::map<std::string, size_t> symbolToSymbolNum;

    std::vector<size_t> atomNumToSymbol(nAtoms);
    std::vector<size_t> atomNumToMonomerNum(nAtoms);
    std::vector<size_t> atomNumToNumInMonomer(nAtoms);
    std::vector<float> atomX(nAtoms);
    std::vector<float> atomY(nAtoms);

    std::vector<size_t> bondNumToMonomerNum(nBonds);
    std::vector<size_t> bondNumToNumInMonomer(nBonds);
    std::vector<size_t> bondNumToFirst(nBonds);
    std::vector<size_t> bondNumToSecond(nBonds);
    std::vector<indigo::KetBond::bond_types> bondNumToOrder(nBonds);

    size_t atomCounter = 0;
    size_t previousAddition = 0;
    size_t bondCounter = 0;
    for (size_t i = 0; i < nMonomers; ++i) {
        const std::string code = templatesNumToName[monomerNumToTemplate[i]];

        auto& atoms = library.getMonomerTemplateById(code).atoms();
        auto& bonds = library.getMonomerTemplateById(code).bonds();
        monomerFirstAtom[i] = atomCounter; 

        for (size_t j = 0; j < atoms.size(); ++j) {
            auto atom = (dynamic_cast<KetAtom*>(atoms[j].get()));

            std::string label = atom->label();

            size_t num;
            if (symbolToSymbolNum.find(label) != symbolToSymbolNum.end()) {
                num = symbolToSymbolNum[label];
            } else {
                num = symbolToSymbolNum.size();
                symbolToSymbolNum[label] = num;
                symbolNumToSymbol.push_back(label);
            }
    
            atomNumToSymbol[atomCounter] = num;
            atomNumToMonomerNum[atomCounter] = i;
            atomNumToNumInMonomer[atomCounter] = j;
            atomX[atomCounter] = atom->location()->x;
            atomY[atomCounter] = atom->location()->y;

            ++atomCounter;
        }
    
        monomerLastAtom[i] = atomCounter;

        for (size_t j = 0; j < bonds.size(); ++j) { 
            auto bond = bonds[j];

            bondNumToMonomerNum[bondCounter] = i;
            bondNumToNumInMonomer[bondCounter] = j;
            bondNumToFirst[bondCounter] = bond.atoms().first + previousAddition;
            bondNumToSecond[bondCounter] = bond.atoms().second + previousAddition;
            bondNumToOrder[bondCounter] = bond.getType();

            ++bondCounter;
        }

        previousAddition = atomCounter;
    }

    // std::vector<size_t> monomerNumToName(nMonomers);
    // std::vector<size_t> monomerNumToTemplate(nMonomers);
    // std::vector<float> monomerNumToX(nMonomers);
    // std::vector<float> monomerNumToY(nMonomers);
    // std::vector<size_t> monomerROne(nMonomers);
    // std::vector<size_t> monomerRTwo(nMonomers);
    // std::vector<size_t> monomerFirstAtom(nMonomers);
    // std::vector<size_t> monomerLastAtom(nMonomers);

    // std::vector<float> atomX(nAtoms);
    // std::vector<float> atomY(nAtoms);
    //arrangement
    float lastX = 0;
    for (size_t i = 0; i < nMonomers; ++i) {
        size_t first = monomerFirstAtom[i];
        size_t last = monomerLastAtom[i];
        size_t left = monomerROne[i] + first;
        size_t right = monomerRTwo[i] + first;
        float dy = atomY[right] - atomY[left];
        float dx = atomX[right] - atomX[left]; 
        float xCentre = (atomX[right] + atomX[left])/2;
        float yCentre = (atomX[right] + atomX[left])/2;
        float interAttachSize = sqrt(dy*dy + dx*dx);

        //transition
        float shift = lastX + interAttachSize/2 + STD_BOND;
        for (size_t at = first; at < last; ++at)
            atomX[at] += shift;

        xCentre += shift; 
        
        lastX =  atomX[right];
    }

    //cut and link
    std::vector<size_t> odds;
    for (size_t i = 0; i < nConnections; ++i) {
        std::string m1 = connections[i].ep1().getStringProp("monomerId");//"attachmentPointId"
        size_t m1_num = stoi(strReplace(m1, "monomer", ""));
        std::string m2 = connections[i].ep2().getStringProp("monomerId");//"attachmentPointId"
        size_t m2_num = stoi(strReplace(m2, "monomer", ""));

        const std::string t1 = templatesNumToName[monomerNumToTemplate[m1_num]];
        const std::string t2 = templatesNumToName[monomerNumToTemplate[m2_num]];

        size_t leftLinkAtom = monomerFirstAtom[m1_num] + monomerRTwo[m1_num];
        size_t rightLinkAtom = monomerFirstAtom[m2_num] + monomerROne[m2_num];

        //bondNumToMonomerNum;
        //bondNumToNumInMonomer;
        ++nBonds;
        bondNumToFirst.push_back(leftLinkAtom);
        bondNumToSecond.push_back(rightLinkAtom);;
        bondNumToOrder.push_back(indigo::KetBond::bond_types::single);

        auto oddsLeft = library.getMonomerTemplateById(t1).attachmentPoints().at("R2").leavingGroup().value();
        auto oddsRight = library.getMonomerTemplateById(t2).attachmentPoints().at("R1").leavingGroup().value();

        for (size_t n : oddsLeft)
            odds.push_back(monomerFirstAtom[m1_num] + n);
        for (size_t n : oddsRight)
            odds.push_back(monomerFirstAtom[m2_num] + n);
    }

    sort(odds.begin(), odds.end());

    for (size_t i = 0; i < odds.size(); ++i) {
        size_t odd = odds[i] - i;
        atomNumToSymbol.erase(atomNumToSymbol.begin() + odd);
        atomNumToMonomerNum.erase(atomNumToMonomerNum.begin() + odd);
        atomNumToNumInMonomer.erase(atomNumToNumInMonomer.begin() + odd);
        atomX.erase(atomX.begin() + odd);
        atomY.erase(atomY.begin() + odd);

        for (size_t j = bondNumToFirst.size(); j > 0; --j) {
            if (bondNumToFirst[j - 1] == odd || bondNumToSecond[j - 1] == odd) {
                bondNumToFirst.erase(bondNumToFirst.begin() + j - 1);
                bondNumToSecond.erase(bondNumToSecond.begin() + j - 1);
                bondNumToOrder.erase(bondNumToOrder.begin() + j - 1);
            } else{
                bondNumToFirst[j - 1] = bondNumToFirst[j - 1] > odd ? bondNumToFirst[j - 1] - 1:  bondNumToFirst[j - 1];
                bondNumToSecond[j - 1] = bondNumToSecond[j - 1] > odd ? bondNumToSecond[j - 1] - 1:  bondNumToSecond[j - 1];
            }
        }
    }

    //recover mol block
    std::string molfileCountsLine = V3K_BEGIN_COUNTS_LINE + std::to_string(atomNumToSymbol.size()) + ' ' + std::to_string(bondNumToOrder.size()) + V3K_COUNTS_LINE_ENDING;

    std::string result = "";
    result += V3K_HEADER_FIRST_LINE;
    result += V3K_HEADER_SECOND_LINE;
    result += V3K_BEGIN_CTAB_BLOCK;
    result += molfileCountsLine;
    result += V3K_BEGIN_ATOM_BLOCK;

    for (size_t i = 0; i < atomNumToSymbol.size(); ++i) {
        result += V3K_BEGIN_DATA_LINE + std::to_string(i + 1) +
            " " + symbolNumToSymbol[atomNumToSymbol[i]] +
            " " + std::to_string(atomX[i]) +
            " " + std::to_string(atomY[i]) +
            " 0.0 0\n";
    }

    result += V3K_END_ATOM_BLOCK;
    result += V3K_BEGIN_BOND_BLOCK;

    for (size_t i = 0; i < bondNumToOrder.size(); ++i) {
        result += V3K_BEGIN_DATA_LINE + std::to_string(i + 1) +
            " " + std::to_string(static_cast<int>(bondNumToOrder[i])) +
            " " + std::to_string(bondNumToFirst[i] + 1) +
            " " + std::to_string(bondNumToSecond[i] + 1) +
            "\n";
    }

    result += V3K_END_BOND_BLOCK;
    //steabs collection
    result += V3K_END_CTAB_BLOCK;
    result += V3K_END;

    return result;
}