#include "TopologyReader.h"

void TopologyReader::ReadFile(std::string& name, std::vector<std::string>& contents)
{
    contents.clear();

    // First read everything in the file
    std::ifstream ifs;

    // open the file and make sure it is opened
    ifs.open(name);
    ASSERT((ifs.is_open()), "The name " << name << " is not opened.");

    // start parsing the file 
    std::string sentence;

    // start parsing the sentences 
    while (std::getline(ifs,sentence)){
        // skip the empty lines
        if (! StringTools::CheckIfOnlyWhiteSpace(sentence)){
            std::vector<std::string> words = StringTools::split(sentence);

            // if we have a include file, then we read that file
            if (words[0] == "#include"){
                // remove double quotes from the name
                words[1].erase(std::remove(words[1].begin(), words[1].end(), '\"'), words[1].end());

                // read the file --> always search from the provided path
                std::vector<std::string> content_include;
                std::string path = FileSystem::joinPath(absolute_path_, words[1]);
                ReadFile(path, content_include);

                contents.insert(contents.end(), content_include.begin(), content_include.end());
            }
            else{
                // no # or ; symbol in the first letter of the sentence 
                bool NotComment = words[0].find("#") == std::string::npos && words[0].find(";") == std::string::npos;
                if (NotComment){
                    contents.push_back(sentence);
                }
            }
        }
    }
}

void TopologyReader::Parse(std::string& name)
{
    // extract the absolute path
    std::size_t botdirPos = name.find_last_of("/");
    absolute_path_ = name.substr(0, botdirPos);

    // read the file
    std::vector<std::string> contents;
    ReadFile(name, contents);

    // contents is the content of the file --> all the lines inside the file 
    int moleculeIndices;

    // read all the contents starting with [ ] --> usually indicates a starting section 
    for (int i=0;i<contents.size();i++){
        // see if we can access [ atoms ]
        // find the first one not space 
        std::string stripped = StringTools::strip(contents[i]);
        if (stripped[0] == '['){
            // slit the sentence
            std::vector<std::string> words = StringTools::split(stripped);
            ASSERT((words.size() != 0), "Read an empty line in topology for some reason");

            // if the word is atoms, then we start parsing
            if (words[1] == "atoms"){
                AtomIndices_.push_back(i);
            }

            if (words[1] == "molecules"){
                moleculeIndices = i;
            }

            if (words[1] == "atomtypes"){
                AtomtypeIndices_.push_back(i);
            }

            if (words[1] == "moleculetype"){
                MoleculeTypeIndices_.push_back(i);
            }
        }
    }

    ////////////////////////////////////////   [ moleculetype ]  ///////////////////////////////////////////////////
    for (int i=0;i<MoleculeTypeIndices_.size();i++){
        int start = MoleculeTypeIndices_[i];
        for (int j=start+1;j<contents.size();j++){
            int index = contents[j].find_first_not_of(" ");
            if (contents[j][index] == '['){
                break;
            }
            else{
                std::vector<std::string> words = StringTools::split(contents[j], comment_str_);
                ASSERT(words.size() == 2, "Moleculetype should only have ; name nrexl.");
                unique_resnames_.push_back(words[0]);
            }
        }
    }

    ////////////////////////////////////////    [ molecules ] //////////////////////////////////////////////
    for (int i = moleculeIndices+1;i<contents.size();i++){
        std::string stripped = StringTools::strip(contents[i]);

        if (stripped[0] == '['){
            break;
        }
        else{
            std::vector<std::string> words = StringTools::split(stripped, comment_str_);
            ASSERT((words.size() == 2), "For molecule directives, we must only have 2 entries, [molecuname, number], while we have " << words.size());

            // map resname to the number of residues
            auto it = MapResnameToNumberResidues_.find(words[0]);
            int num = StringTools::StringToType<int>(words[1]);
            MapResnameToNumberResidues_.insert(std::make_pair(words[0], num));

            // record the resnames 
            std::vector<std::string> Residues(num, words[0]);
            resnames_.insert(resnames_.end(), Residues.begin(), Residues.end());
        }
    }


    ASSERT((MoleculeTypeIndices_.size() == AtomIndices_.size()), "There must be the same number of moleculetypes and atoms directives.");

    ////////////////////////////////////////    [ atomtypes ] //////////////////////////////////////////////
    for (int i=0;i<AtomtypeIndices_.size();i++){
        int start = AtomtypeIndices_[i];

        for (int j=start+1;j<contents.size();j++){
            int index = contents[j].find_first_not_of(" ");
            if (contents[j][index] == '['){
                break;
            }
            else{
                std::vector<std::string> words =  StringTools::split(contents[j], comment_str_);

                Molecule::AtomType a;
                if (words.size() == 8){
                    a.type_     = words[0];
                    a.charge_   = StringTools::StringToType<Real>(words[4]);
                    a.mass_     = StringTools::StringToType<Real>(words[3]);
                    a.sigma_    = StringTools::StringToType<Real>(words[6]);
                    a.epsilon_  = StringTools::StringToType<Real>(words[7]);
                }
                else if (words.size() == 7){
                    a.type_     = words[0];
                    a.mass_     = StringTools::StringToType<Real>(words[2]);
                    a.charge_   = StringTools::StringToType<Real>(words[3]);
                    a.sigma_    = StringTools::StringToType<Real>(words[5]);
                    a.epsilon_  = StringTools::StringToType<Real>(words[6]);
                }
                else if (words.size() == 6){
                    a.type_     = words[0];
                    a.mass_     = StringTools::StringToType<Real>(words[1]);
                    a.charge_   = StringTools::StringToType<Real>(words[2]);
                    a.sigma_    = StringTools::StringToType<Real>(words[4]);
                    a.epsilon_  = StringTools::StringToType<Real>(words[5]);
                }
                else{ASSERT((true == false), "Something went wrong in the atomtypes reading.");}
                bool mapped = Algorithm::InsertInMap(MapTypenameToAtomType_, a.type_, a);
                std::cout << "mapped = " << mapped << "\n";
                ASSERT(mapped, "The typename " << a.type_ << " is repeated more than once.");
            }

        }
    }

    // read the [ atoms ] directives --> see if there's any differences in charge and mass
    for (int i =0;i<AtomIndices_.size();i++){
        int start = AtomIndices_[i];
        for (int j=start+1;j<contents.size();j++){
            int index = contents[j].find_first_not_of(" ");
            if (contents[j][index] == '['){break;}
            else{
                std::vector<std::string> words = StringTools::split(contents[j], comment_str_);

                // in [ atoms ] directive, sometimes we also have information about the atom types
                // we check if the current mass or charge in atom type is zero
                std::string rname = unique_resnames_[i];
                std::string tname= words[TopologyItemIndex::atomtype];

                // if words size >= 7, that means we have charge info --> compare to the atomtype
                if (words.size() >= 7){
                    Real charge = StringTools::StringToType<Real>(words[TopologyItemIndex::charge]);

                    if (MapTypenameToAtomType_.find(tname)->second.charge_ == 0){
                        MapTypenameToAtomType_.find(tname)->second.charge_ = charge;
                    }
                }

                // if words size >= 8, then we have mass info
                if (words.size() >= 8){
                    Real mass = StringTools::StringToType<Real>(words[TopologyItemIndex::mass]);

                    if (MapTypenameToAtomType_.find(tname)->second.mass_ == 0){
                        MapTypenameToAtomType_.find(tname)->second.mass_ = mass;
                    }
                }

                // make a map of residue name to type names 
                auto it = MapResnameToTypename_.find(rname);
                if (it != MapResnameToTypename_.end()){
                    it -> second.push_back(tname);
                }
                else{
                    std::vector<std::string> tname = {words[TopologyItemIndex::atomtype]};
                    MapResnameToTypename_.insert(std::make_pair(rname, tname));
                }

                auto itatom = MapResnameToAtomname_.find(rname);
                if (itatom != MapResnameToAtomname_.end()){
                    itatom -> second.push_back(words[TopologyItemIndex::atomName]);
                }
                else{
                    std::vector<std::string> aname = {words[TopologyItemIndex::atomName]};
                    MapResnameToAtomname_.insert(std::make_pair(rname, aname));
                }
            }
        }        
    }

    MapIndicesToAtom();
}

void TopologyReader::MapIndicesToAtom(){
    int numatoms = 1;
    for (int i=0;i<resnames_.size();i++){
        std::string name = resnames_[i];

        std::vector<std::string> atomtypename;
        bool mapped = Algorithm::FindInMap(MapResnameToTypename_, name, atomtypename);
        ASSERT(mapped, "The resname " << name << " is not found.");

        auto atomname = MapResnameToAtomname_.find(name)->second;

        std::vector<Molecule::atom> AtomVec;
        Molecule::residue res;
        for (int j=0;j<atomname.size();j++){
            auto it = MapTypenameToAtomType_.find(atomtypename[j]);
            ASSERT((it != MapTypenameToAtomType_.end()), "The type name " << atomtypename[j] << " in resname " << name << " is not found.");
            auto atype = it->second;

            Molecule::atom a;
            a.atomName_ = atomname[j];
            a.atomNumber_ = numatoms;  
            numatoms++;
            a.charge_   = atype.charge_;
            a.mass_     = atype.mass_;
            a.resnum_ = i+1;
            a.resname_ = name;

            AtomVec.push_back(a);
        }

        res.atoms_ = AtomVec;
        atoms_.insert(atoms_.end(), AtomVec.begin(), AtomVec.end());
        residues_.push_back(res);
    }
}
