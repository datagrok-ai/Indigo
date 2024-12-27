#ifdef _DEBUG
#include <crtdbg.h>
#define _CRTDBG_MAP_ALLOC
#endif
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

#include "common.h"

#include <algorithm>

using namespace indigo;

int main()
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

    return 0;
}
