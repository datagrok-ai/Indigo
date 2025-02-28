/****************************************************************************
 * Copyright (C) from 2009 to Present EPAM Systems.
 *
 * This file is part of Indigo toolkit.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ***************************************************************************/

#ifdef _MSC_VER
#pragma warning(push)
#endif

#include <gtest/gtest.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

#include "molecule/ket_document.h"
#include "molecule/ket_document_json_loader.h"
#include "molecule/ket_document_json_saver.h"
#include "molecule/monomers_template_library.h"
#include "molecule/sequence_loader.h"
#include "molecule/sequence_saver.h"
#include <base_cpp/output.h>
#include <base_cpp/scanner.h>
#include <molecule/cmf_loader.h>
#include <molecule/cmf_saver.h>
#include <molecule/cml_saver.h>
#include <molecule/molecule_cdxml_saver.h>
#include <molecule/molecule_json_loader.h>
#include <molecule/molecule_json_saver.h>
#include <molecule/molecule_mass.h>
#include <molecule/molecule_substructure_matcher.h>
#include <molecule/molfile_loader.h>
#include <molecule/molfile_saver.h>
#include <molecule/query_molecule.h>
#include <molecule/sdf_loader.h>
#include <molecule/smiles_loader.h>
#include <molecule/smiles_saver.h>

#include <molecule/atomic_level.h>

#include "common.h"

#include <algorithm>

using namespace indigo;

class IndigoCoreFormatsTest : public IndigoCoreTest
{
};

TEST_F(IndigoCoreFormatsTest, load_targets_cmf)
{
    FileScanner sc(dataPath("molecules/resonance/resonance.sdf").c_str());

    SdfLoader sdf(sc);
    QueryMolecule qmol;

    Array<char> qbuf;
    qbuf.readString("N(#C)=C(C)C", false);
    BufferScanner sm_scanner(qbuf);
    SmilesLoader smiles_loader(sm_scanner);
    smiles_loader.loadQueryMolecule(qmol);

    sdf.readAt(138);
    try
    {
        BufferScanner bsc(sdf.data);
        MolfileLoader loader(bsc);
        Molecule mol;
        loader.loadMolecule(mol);
        Array<char> buf;
        ArrayOutput buf_out(buf);
        CmfSaver cmf_saver(buf_out);

        cmf_saver.saveMolecule(mol);

        Molecule mol2;
        BufferScanner buf_in(buf);
        CmfLoader cmf_loader(buf_in);
        cmf_loader.loadMolecule(mol2);

        MoleculeSubstructureMatcher matcher(mol2);
        matcher.use_pi_systems_matcher = true;
        matcher.setQuery(qmol);
        matcher.find();
    }
    catch (Exception& e)
    {
        ASSERT_STREQ("", e.message());
    }
}

TEST_F(IndigoCoreFormatsTest, save_cdxml)
{
    Molecule t_mol;

    loadMolecule("c1ccccc1N", t_mol);

    Array<char> out;
    ArrayOutput std_out(out);
    MoleculeCdxmlSaver saver(std_out);
    saver.saveMolecule(t_mol);
    loadMolecule("c1ccccc1", t_mol);
    saver.saveMolecule(t_mol);

    ASSERT_TRUE(out.size() > 2000);
}

TEST_F(IndigoCoreFormatsTest, save_cml)
{
    Molecule t_mol;

    loadMolecule("c1ccccc1N", t_mol);

    Array<char> out;
    ArrayOutput std_out(out);
    CmlSaver saver(std_out);
    saver.saveMolecule(t_mol);
    loadMolecule("c1ccccc1", t_mol);
    saver.saveMolecule(t_mol);

    ASSERT_TRUE(out.size() > 1000);
}

TEST_F(IndigoCoreFormatsTest, smiles_data_sgroups)
{
    Molecule t_mol;

    loadMolecule("c1ccccc1N |SgD:3,2,1,0:name:data:like:unit:t:(-1)|", t_mol);
    ASSERT_EQ(t_mol.sgroups.getSGroupCount(), 1);
    SGroup& sg = t_mol.sgroups.getSGroup(0);
    ASSERT_EQ(sg.sgroup_type, SGroup::SG_TYPE_DAT);
    DataSGroup& dsg = static_cast<DataSGroup&>(sg);
    ASSERT_STREQ(dsg.name.ptr(), "name");
    ASSERT_STREQ(dsg.data.ptr(), "data");
    ASSERT_STREQ(dsg.queryoper.ptr(), "like");
    ASSERT_STREQ(dsg.description.ptr(), "unit");
    ASSERT_EQ(dsg.tag, 't');
    ASSERT_EQ(dsg.display_pos.x, 0.0f);
    ASSERT_EQ(dsg.display_pos.y, 0.0f);
    ASSERT_EQ(dsg.atoms.size(), 4);
    ASSERT_EQ(dsg.atoms.at(0), 3);
    ASSERT_EQ(dsg.atoms.at(1), 2);
    ASSERT_EQ(dsg.atoms.at(2), 1);
    ASSERT_EQ(dsg.atoms.at(3), 0);
    Array<char> out;
    ArrayOutput std_out(out);
    SmilesSaver saver(std_out);
    saver.saveMolecule(t_mol);
    ASSERT_EQ(out.size(), 48);
    std::string str{out.ptr(), static_cast<std::size_t>(out.size())};
    ASSERT_STREQ(str.c_str(), "c1c(N)cccc1 |SgD:3,2,1,0:name:data:like:unit:t:|");
}

TEST_F(IndigoCoreFormatsTest, smiles_data_sgroups_coords)
{
    Molecule t_mol;

    loadMolecule("c1ccccc1 |SgD:1,2,0:::::s:(-1.5,7.8)|", t_mol);
    ASSERT_EQ(t_mol.sgroups.getSGroupCount(), 1);
    SGroup& sg = t_mol.sgroups.getSGroup(0);
    ASSERT_EQ(sg.sgroup_type, SGroup::SG_TYPE_DAT);
    DataSGroup& dsg = static_cast<DataSGroup&>(sg);
    ASSERT_STREQ(dsg.name.ptr(), "");
    ASSERT_STREQ(dsg.data.ptr(), "");
    ASSERT_STREQ(dsg.queryoper.ptr(), "");
    ASSERT_STREQ(dsg.description.ptr(), "");
    ASSERT_EQ(dsg.tag, 's');
    ASSERT_EQ(dsg.display_pos.x, -1.5f);
    ASSERT_EQ(dsg.display_pos.y, 7.8f);
    ASSERT_EQ(dsg.atoms.size(), 3);
    ASSERT_EQ(dsg.atoms.at(0), 1);
    ASSERT_EQ(dsg.atoms.at(1), 2);
    ASSERT_EQ(dsg.atoms.at(2), 0);
    Array<char> out;
    ArrayOutput std_out(out);
    SmilesSaver saver(std_out);
    saver.saveMolecule(t_mol);
    ASSERT_EQ(out.size(), 27);
    std::string str{out.ptr(), static_cast<std::size_t>(out.size())};
    ASSERT_STREQ(str.c_str(), "c1ccccc1 |SgD:1,2,0:::::s:|");
}

TEST_F(IndigoCoreFormatsTest, smiles_data_sgroups_short)
{
    Molecule t_mol;

    loadMolecule("c1ccccc1 |SgD:1,2,0:name|", t_mol);
    ASSERT_EQ(t_mol.sgroups.getSGroupCount(), 1);
    SGroup& sg = t_mol.sgroups.getSGroup(0);
    ASSERT_EQ(sg.sgroup_type, SGroup::SG_TYPE_DAT);
    DataSGroup& dsg = static_cast<DataSGroup&>(sg);
    ASSERT_STREQ(dsg.name.ptr(), "name");
    ASSERT_EQ(dsg.data.size(), 0);
    ASSERT_EQ(dsg.queryoper.size(), 0);
    ASSERT_EQ(dsg.description.size(), 0);
    ASSERT_EQ(dsg.tag, ' ');
    ASSERT_EQ(dsg.display_pos.x, 0.0f);
    ASSERT_EQ(dsg.display_pos.y, 0.0f);
    ASSERT_EQ(dsg.atoms.size(), 3);
    ASSERT_EQ(dsg.atoms.at(0), 1);
    ASSERT_EQ(dsg.atoms.at(1), 2);
    ASSERT_EQ(dsg.atoms.at(2), 0);
    Array<char> out;
    ArrayOutput std_out(out);
    SmilesSaver saver(std_out);
    saver.saveMolecule(t_mol);
    ASSERT_EQ(out.size(), 31);
    std::string str{out.ptr(), static_cast<std::size_t>(out.size())};
    ASSERT_STREQ(str.c_str(), "c1ccccc1 |SgD:1,2,0:name:::: :|");
}

TEST_F(IndigoCoreFormatsTest, smiles_pol_sgroups_conn_and_flip)
{
    Molecule t_mol;

    loadMolecule("*CC*C*N* |$star;;;star;;star;;star$,Sg:n:6,1,2,4::hh&#44;f:6,0,:4,2,|", t_mol);
    ASSERT_EQ(t_mol.sgroups.getSGroupCount(), 1);
    SGroup& sg = t_mol.sgroups.getSGroup(0);
    ASSERT_EQ(sg.sgroup_type, SGroup::SG_TYPE_SRU);
    RepeatingUnit& ru = static_cast<RepeatingUnit&>(sg);
    ASSERT_EQ(ru.atoms.size(), 4);
    ASSERT_EQ(ru.atoms.at(0), 6);
    ASSERT_EQ(ru.atoms.at(1), 1);
    ASSERT_EQ(ru.atoms.at(2), 2);
    ASSERT_EQ(ru.atoms.at(3), 4);
    ASSERT_EQ(ru.connectivity, RepeatingUnit::HEAD_TO_HEAD);
    Array<char> out;
    ArrayOutput std_out(out);
    SmilesSaver saver(std_out);
    saver.saveMolecule(t_mol);
    ASSERT_EQ(out.size(), 53);
    std::string str{out.ptr(), static_cast<std::size_t>(out.size())};
    ASSERT_STREQ(str.c_str(), "*CC*C*N* |$star;;;star;;star;;star$,Sg:n:6,1,2,4::hh|");
}

TEST_F(IndigoCoreFormatsTest, smiles_pol_sgroups_bracket)
{
    Molecule t_mol;

    loadMolecule("C1CCCCC1 |Sg:n:0,5,4,3,2,1:::::(d,s,-7.03,2.12,-2.21,2.12,-2.21,-3.11,-7.03,-3.11,)|", t_mol);
    ASSERT_EQ(t_mol.sgroups.getSGroupCount(), 1);
    SGroup& sg = t_mol.sgroups.getSGroup(0);
    ASSERT_EQ(sg.sgroup_type, SGroup::SG_TYPE_SRU);
    RepeatingUnit& ru = static_cast<RepeatingUnit&>(sg);
    ASSERT_EQ(ru.atoms.size(), 6);
    ASSERT_EQ(ru.atoms.at(0), 0);
    ASSERT_EQ(ru.atoms.at(1), 5);
    ASSERT_EQ(ru.atoms.at(2), 4);
    ASSERT_EQ(ru.atoms.at(3), 3);
    ASSERT_EQ(ru.atoms.at(4), 2);
    ASSERT_EQ(ru.atoms.at(5), 1);
    Array<char> out;
    ArrayOutput std_out(out);
    SmilesSaver saver(std_out);
    saver.saveMolecule(t_mol);
    ASSERT_EQ(out.size(), 31);
    std::string str{out.ptr(), static_cast<std::size_t>(out.size())};
    ASSERT_STREQ(str.c_str(), "C1CCCCC1 |Sg:n:0,5,4,3,2,1::eu|");
}

TEST_F(IndigoCoreFormatsTest, smiles_pol_sgroups_gen)
{
    Molecule t_mol;

    loadMolecule("CCCC |Sg:gen:0,1,2:|", t_mol);
    ASSERT_EQ(t_mol.sgroups.getSGroupCount(), 1);
    SGroup& sg = t_mol.sgroups.getSGroup(0);
    ASSERT_EQ(sg.sgroup_type, SGroup::SG_TYPE_GEN);
    Array<char> out;
    ArrayOutput std_out(out);
    SmilesSaver saver(std_out);
    saver.saveMolecule(t_mol);
    ASSERT_EQ(out.size(), 20);
    std::string str{out.ptr(), static_cast<std::size_t>(out.size())};
    ASSERT_STREQ(str.c_str(), "CCCC |Sg:gen:0,1,2:|");
}

TEST_F(IndigoCoreFormatsTest, mol_saver_issue_1200)
{
    Molecule t_mol;

    const char* mol = R"(
  -INDIGO-07262316452D

  6  6  0  0  0  0  0  0  0  0999 V2000
   -1.4617   -0.6508    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3856   -0.8000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.5747    0.2706    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.9239    1.7323    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3327    1.5650    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  4  0  0  0  0
  3  4  4  0  0  0  0
  4  5  4  0  0  0  0
  5  6  4  0  0  0  0
  6  2  4  0  0  0  0
M  STY  4   1 DAT   2 DAT   3 DAT   4 DAT
M  SLB  4   1   1   2   2   3   3   4   4
M  SAL   1  1   3
M  SDT   1 MRV_IMPLICIT_H                                                       
M  SDD   1     0.0000    0.0000    DA    ALL  1       1  
M  SED   1 IMPL_H1
M  SAL   2  1   4
M  SDT   2 MRV_IMPLICIT_H                                                       
M  SDD   2     0.0000    0.0000    DA    ALL  1       1  
M  SED   2 IMPL_H1
M  SAL   3  1   3
M  SDT   3 MRV_IMPLICIT_H                                                       
M  SDD   3     0.0000    0.0000    DA    ALL  1       1  
M  SED   3 IMPL_H1
M  SAL   4  1   4
M  SDT   4 MRV_IMPLICIT_H                                                       
M  SDD   4     0.0000    0.0000    DA    ALL  1       1  
M  SED   4 IMPL_H1
M  END
)";
    loadMolecule(mol, t_mol);
    ASSERT_EQ(t_mol.sgroups.getSGroupCount(), 0);
    Array<char> out;
    ArrayOutput std_out(out);
    MolfileSaver saver(std_out);
    saver.saveMolecule(t_mol);
    ASSERT_EQ(t_mol.sgroups.getSGroupCount(), 0);
    saver.mode = MolfileSaver::MODE_2000;
    saver.saveMolecule(t_mol);
    ASSERT_EQ(t_mol.sgroups.getSGroupCount(), 0);
    saver.mode = MolfileSaver::MODE_3000;
    saver.saveMolecule(t_mol);
    ASSERT_EQ(t_mol.sgroups.getSGroupCount(), 0);
}

TEST_F(IndigoCoreFormatsTest, smarts_load_save)
{
    QueryMolecule q_mol;

    std::string smarts_in{"([#8].[#6]).([#6].[#8])"};
    BufferScanner scanner(smarts_in.c_str());
    SmilesLoader loader(scanner);
    loader.smarts_mode = true;
    loader.loadQueryMolecule(q_mol);
    Array<char> out;
    ArrayOutput std_out(out);
    SmilesSaver saver(std_out);
    saver.smarts_mode = true;
    saver.saveQueryMolecule(q_mol);
    std::string smarts_out{out.ptr(), static_cast<std::size_t>(out.size())};
    ASSERT_EQ(smarts_in, smarts_out);
}

TEST_F(IndigoCoreFormatsTest, json_load_save)
{
    QueryMolecule q_mol;

    FileScanner sc(dataPath("molecules/basic/ket_with_query_properties.ket").c_str());
    std::string json;
    sc.readAll(json);
    rapidjson::Document data;
    if (!data.Parse(json.c_str()).HasParseError())
    {
        if (data.HasMember("root"))
        {
            MoleculeJsonLoader loader(data);
            loader.loadMolecule(q_mol);
        }
    }

    Array<char> out;
    ArrayOutput std_out(out);
    MoleculeJsonSaver saver(std_out);
    saver.pretty_json = true;
    saver.saveMolecule(q_mol);
    std::string json_out{out.ptr(), static_cast<std::size_t>(out.size())};
    // ASSERT_EQ(json, json_out);
}

TEST_F(IndigoCoreFormatsTest, idt_load)
{
    const char* idt = "ARAS";
    BufferScanner scanner(idt);
    MonomerTemplateLibrary library;
    FileScanner sc(dataPath("molecules/basic/monomer_library.ket").c_str());
    std::string json;
    sc.readAll(json);
    rapidjson::Document data;
    if (!data.Parse(json.c_str()).HasParseError())
    {
        if (data.HasMember("root"))
        {
            MoleculeJsonLoader loader(data);
            loader.loadMonomerLibrary(library);
        }
    }

    SequenceLoader loader(scanner, library);
    KetDocument document;
    loader.loadIdt(document);

    Array<char> out;
    ArrayOutput std_out(out);
    KetDocumentJsonSaver saver(std_out);
    saver.pretty_json = true;
    saver.saveKetDocument(document);
    std::string json_out{out.ptr(), static_cast<std::size_t>(out.size())};
    // printf("%s", json_out.c_str());
    Array<char> buf;
    ArrayOutput output(buf);
    SequenceSaver idt_saver(output, library);
    FileScanner ket(dataPath("molecules/basic/idt_mixed_std.ket").c_str());
    std::string json2;
    ket.readAll(json2);
    json2.erase(std::remove(json2.begin(), json2.end(), '\r'), json2.end());
    ASSERT_EQ(json2, json_out);
}

TEST_F(IndigoCoreFormatsTest, idt_save)
{
    MonomerTemplateLibrary library;
    FileScanner sc(dataPath("molecules/basic/monomer_library.ket").c_str());
    std::string json;
    sc.readAll(json);
    rapidjson::Document data;
    if (!data.Parse(json.c_str()).HasParseError())
    {
        if (data.HasMember("root"))
        {
            MoleculeJsonLoader loader(data);
            loader.loadMonomerLibrary(library);
        }
    }

    Array<char> buf;
    ArrayOutput output(buf);
    SequenceSaver idt_saver(output, library);
    FileScanner ket(dataPath("molecules/basic/idt_mixed_std.ket").c_str());
    std::string json2;
    ket.readAll(json2);
    KetDocument ket_document;
    KetDocumentJsonLoader kdloader;
    kdloader.parseJson(json2, ket_document);

    idt_saver.saveKetDocument(ket_document, SequenceSaver::SeqFormat::IDT);

    std::string json_out{buf.ptr(), static_cast<std::size_t>(buf.size())};
    // printf("%s", json_out.c_str());
    ASSERT_EQ("ARAS", json_out);
}

TEST_F(IndigoCoreFormatsTest, mol_to_document)
{
    MonomerTemplateLibrary library;
    FileScanner sc(dataPath("molecules/basic/chem_peptide.mol").c_str());
    MolfileLoader loader(sc);
    Molecule mol;
    loader.loadMolecule(mol);

    KetDocument& document = mol.getKetDocument();
    std::vector<std::deque<std::string>> sequences;
    // parse connections between monomers and return backbone sequences
    // non-backbone connection stored in document.nonSequenceConnections()
    document.parseSimplePolymers(sequences);
    int idx = 1;
    for (auto& sequence : sequences)
    {
        printf("Backbone %d\n", idx++);
        for (auto monomer_id : sequence)
        {

            const std::unique_ptr<KetBaseMonomer>& monomer = document.monomers().at(monomer_id);
            MonomerClass monomer_class = document.getMonomerClass(*monomer);
            KetBaseMonomer::MonomerType monomer_type = monomer->monomerType();
            const std::optional<Vec2f>& position = monomer->position();
            Vec2f pos = position.has_value() ? position.value() : Vec2f{0, 0};
            printf("monomer %s\tclass=%s\tposition: %f,%f\n", monomer->alias().c_str(), MonomerTemplate::MonomerClassToStr(monomer_class).c_str(), pos.x,
                   pos.y);
            if (monomer_type == KetBaseMonomer::MonomerType::AmbiguousMonomer)
            {
                const KetAmbiguousMonomerTemplate& templ = document.ambiguousTemplates().at(monomer->templateId());
                printf("ambigous monomer, templates: ");
                for (const KetAmbiguousMonomerOption& option : templ.options())
                {
                    printf("%s ", option.templateId().c_str());
                }
                printf("\n");
            }
            else if (monomer_type == KetBaseMonomer::MonomerType::Monomer)
            {
                const MonomerTemplate& templ = document.templates().at(monomer->templateId());
                printf("monomer, template: %s\n", templ.getStringProp("alias").c_str());
                int atom_idx = 0;
                printf("Atoms:\n");
                for (const std::shared_ptr<KetBaseAtomType>& batom : templ.atoms())
                {
                    printf("%d: ", atom_idx++);
                    if (batom->getType() == KetBaseAtom::atype::atom)
                    {
                        const KetAtom& atom = static_cast<KetAtom&>(*batom);
                        printf("atom %s", atom.label().c_str());
                    }
                    else if (batom->getType() == KetBaseAtom::atype::atom_list)
                    {
                        const KetAtomList& atom_list = static_cast<KetAtomList&>(*batom);
                        printf("atom list: ");
                        for (std::string atom : atom_list.atomList())
                        {
                            printf("%s ", atom.c_str());
                        }
                    }
                    const std::optional<Vec3f>& location = batom->location();
                    if (location.has_value())
                        printf(" (%f, %f, %f)", location->x, location->y, location->z);
                    printf("\n");
                }
                printf("Bonds:\n");
                for (const KetBond& bond : templ.bonds())
                {
                    const std::pair<int, int>& atoms = bond.atoms();
                    printf("bond: type=%d, atoms=(%d, %d)\n", bond.getType(), atoms.first, atoms.second);
                }
            }
        }
    }
    if (document.nonSequenceConnections().size() > 0)
        printf("Non-standard connections:\n");
    for (const KetConnection& connection : document.nonSequenceConnections())
    {
        auto fill_conn_info = [](const KetConnectionEndPoint& ep, std::string& info) {
            if (ep.hasStringProp("monomerId"))
            {
                info += "monomer " + ep.getStringProp("monomerId");
                if (ep.hasStringProp("attachmentPointId"))
                {
                    info += " attachment point " + ep.getStringProp("attachmentPointId");
                }
            }
            else if (ep.hasStringProp("moleculeId"))
            {
                info += "molecule " + ep.getStringProp("moleculeId");
                if (ep.hasStringProp("atomId"))
                {
                    info += " atom " + ep.getStringProp("atomId");
                }
            }
        };
        std::string ep1_info, ep2_info;
        fill_conn_info(connection.ep1(), ep1_info);
        fill_conn_info(connection.ep2(), ep2_info);
        printf("%s connection from %s to %s\n", connection.connectionType().c_str(), ep1_info.c_str(), ep2_info.c_str());
    }
}

TEST_F(IndigoCoreFormatsTest, atomic_level)
{
   MonomerTemplateLibrary library;
   std::string json = "{\"root\":{\"nodes\":[{\"$ref\":\"monomer0\"},{\"$ref\":\"monomer1\"},{\"$ref\":\"monomer2\"},{\"$ref\":\"monomer3\"}],\"connections\":[{\"connectionType\":\"single\",\"endpoint1\":{\"monomerId\":\"monomer0\",\"attachmentPointId\":\"R2\"},\"endpoint2\":{\"monomerId\":\"monomer1\",\"attachmentPointId\":\"R1\"}},{\"connectionType\":\"single\",\"endpoint1\":{\"monomerId\":\"monomer1\",\"attachmentPointId\":\"R2\"},\"endpoint2\":{\"monomerId\":\"monomer2\",\"attachmentPointId\":\"R1\"}},{\"connectionType\":\"single\",\"endpoint1\":{\"monomerId\":\"monomer2\",\"attachmentPointId\":\"R2\"},\"endpoint2\":{\"monomerId\":\"monomer3\",\"attachmentPointId\":\"R1\"}}],\"templates\":[{\"$ref\":\"monomerTemplate-A___Alanine\"}]},\"monomer0\":{\"type\":\"monomer\",\"id\":\"0\",\"position\":{\"x\":1.25,\"y\":-1.25},\"alias\":\"A\",\"templateId\":\"A___Alanine\"},\"monomerTemplate-A___Alanine\":{\"type\":\"monomerTemplate\",\"atoms\":[{\"label\":\"N\",\"location\":[-1.2549,-0.392,0]},{\"label\":\"C\",\"location\":[-0.272,0.2633,0],\"stereoLabel\":\"abs\"},{\"label\":\"C\",\"location\":[-0.3103,1.7393,0]},{\"label\":\"C\",\"location\":[1.0523,-0.392,0]},{\"label\":\"O\",\"location\":[1.0829,-1.5722,0]},{\"label\":\"O\",\"location\":[2.0353,0.2633,0]},{\"label\":\"H\",\"location\":[-2.3334,0.0905,0]}],\"bonds\":[{\"type\":1,\"atoms\":[1,0]},{\"type\":1,\"atoms\":[1,2],\"stereo\":1},{\"type\":1,\"atoms\":[1,3]},{\"type\":2,\"atoms\":[3,4]},{\"type\":1,\"atoms\":[3,5]},{\"type\":1,\"atoms\":[0,6]}],\"class\":\"AminoAcid\",\"classHELM\":\"PEPTIDE\",\"id\":\"A___Alanine\",\"fullName\":\"Alanine\",\"alias\":\"A\",\"attachmentPoints\":[{\"attachmentAtom\":0,\"leavingGroup\":{\"atoms\":[6]},\"type\":\"left\"},{\"attachmentAtom\":3,\"leavingGroup\":{\"atoms\":[5]},\"type\":\"right\"}],\"naturalAnalogShort\":\"A\"},\"monomer1\":{\"type\":\"monomer\",\"id\":\"1\",\"position\":{\"x\":2.75,\"y\":-1.25},\"alias\":\"A\",\"templateId\":\"A___Alanine\"},\"monomer2\":{\"type\":\"monomer\",\"id\":\"2\",\"position\":{\"x\":4.25,\"y\":-1.25},\"alias\":\"A\",\"templateId\":\"A___Alanine\"},\"monomer3\":{\"type\":\"monomer\",\"id\":\"3\",\"position\":{\"x\":5.75,\"y\":-1.25},\"alias\":\"A\",\"templateId\":\"A___Alanine\"}}";
   rapidjson::Document data;
   if (!data.Parse(json.c_str()).HasParseError())
   {
       if (data.HasMember("root"))
       {
           MoleculeJsonLoader loader(data);
           loader.loadMonomerLibrary(library);
       }
   }

   KetDocument ket_document;
   KetDocumentJsonLoader kdloader;
   kdloader.parseJson(json, ket_document);

   std:: string res = toAtomicLevel(ket_document, library);
   //std::cout << res << std::endl;

   std::string check = "\n Indigo \n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 22 21 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 N 1.438700 -0.392000 0.0 0\nM  V30 2 C 2.421600 0.263300 0.0 0\nM  V30 3 C 2.383300 1.739300 0.0 0\nM  V30 4 C 3.745900 -0.392000 0.0 0\nM  V30 5 O 3.776500 -1.572200 0.0 0\nM  V30 6 H 0.360200 0.090500 0.0 0\nM  V30 7 N 5.184600 -0.392000 0.0 0\nM  V30 8 C 6.167500 0.263300 0.0 0\nM  V30 9 C 6.129200 1.739300 0.0 0\nM  V30 10 C 7.491800 -0.392000 0.0 0\nM  V30 11 O 7.522400 -1.572200 0.0 0\nM  V30 12 N 8.930500 -0.392000 0.0 0\nM  V30 13 C 9.913400 0.263300 0.0 0\nM  V30 14 C 9.875100 1.739300 0.0 0\nM  V30 15 C 11.237700 -0.392000 0.0 0\nM  V30 16 O 11.268300 -1.572200 0.0 0\nM  V30 17 N 12.676399 -0.392000 0.0 0\nM  V30 18 C 13.659299 0.263300 0.0 0\nM  V30 19 C 13.620999 1.739300 0.0 0\nM  V30 20 C 14.983599 -0.392000 0.0 0\nM  V30 21 O 15.014199 -1.572200 0.0 0\nM  V30 22 O 15.966599 0.263300 0.0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 2 1\nM  V30 2 1 2 3\nM  V30 3 1 2 4\nM  V30 4 2 4 5\nM  V30 5 1 1 6\nM  V30 6 1 8 7\nM  V30 7 1 8 9\nM  V30 8 1 8 10\nM  V30 9 2 10 11\nM  V30 10 1 13 12\nM  V30 11 1 13 14\nM  V30 12 1 13 15\nM  V30 13 2 15 16\nM  V30 14 1 18 17\nM  V30 15 1 18 19\nM  V30 16 1 18 20\nM  V30 17 2 20 21\nM  V30 18 1 20 22\nM  V30 19 1 4 7\nM  V30 20 1 10 12\nM  V30 21 1 15 17\nM  V30 END BOND\nM  V30 END CTAB\nM  END";

   ASSERT_EQ(res, check);
}

#ifdef _MSC_VER
#pragma warning(pop)
#endif
