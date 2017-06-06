/* 
 * File:   GenomeFactory.h
 * Author: veraalva
 *
 * Created on May 5, 2016, 1:41 PM
 */

#ifndef GENOMEFACTORY_H
#define GENOMEFACTORY_H

#include <unordered_map>
#include <set>

#include "TextParser.h"
#include "bstring.h"


namespace genome {

    template<typename T>
    class Feature {
    public:

        Feature(std::string type) {
            this->type = type;
            this->start = -1;
            this->end = -1;
            this->processed = false;
        }

        Feature(std::string type, unsigned int start, unsigned int end) {
            this->type = type;
            this->start = start;
            this->end = end;
            this->hash = (unsigned long long) start << 32 | end;
            this->processed = false;
        }

        virtual ~Feature() {
        }

        T& getData() {
            return data;
        }

        void setData(T data) {
            this->data = data;
        }

        unsigned int getEnd() const {
            return end;
        }

        void setEnd(unsigned int end) {
            this->end = end;
            this->hash = (unsigned long long) this->start << 32 | this->end;
        }

        unsigned int getStart() const {
            return start;
        }

        void setStart(unsigned int start) {
            this->start = start;
            this->hash = (unsigned long long) this->start << 32 | this->end;
        }

        std::string getType() const {
            return type;
        }

        void setType(std::string type) {
            this->type = type;
        }

        unsigned int getLength() const {
            return this->end - this->start + 1;
        }

        unsigned long long getHash() const {
            return hash;
        }

        bool isProcessed() const {
            return processed;
        }

        void setProcessed(bool processed) {
            this->processed = processed;
        }

        char getStrand() const {
            return strand;
        }

        void setStrand(char strand) {
            this->strand = strand;
        }

        bool isInside(unsigned int start, unsigned int end, unsigned int overlap) {
            if ((start >= this->start && start <= (this->end - overlap)) ||
                    (end >= (this->start + overlap) && end <= this->end) ||
                    (start <= this->start && end >= this->end))return true;
            return false;
        }

        bool operator==(const Feature& right) const {
            return this->type.compare(right.getType()) == 0 &&
                    this->hash == right.getHash();
        }

        bool operator!=(const Feature& right) const {
            return !(*this == right);
        }

        bool operator>(const Feature& right) const {
            return this->start > right.getStart();
        }

        bool operator<(const Feature& right) const {
            return right > * this;
        }

        friend std::ostream& operator<<(std::ostream& os, Feature *obj) {
            os << obj->getType() << " " << obj->getStrand();
            os << " coord: " << obj->getStart() << "-" << obj->getEnd();
            return os;
        }

    private:
        bool processed;
        std::string type;
        char strand;
        unsigned int start;
        unsigned int end;
        T data;
        unsigned long long hash;
    };

    template<typename T>
    using SPtrFeature = std::shared_ptr<Feature<T>>;

    template<typename T>
    class Isoform : public Feature<T> {
    public:

        Isoform(std::string id, unsigned int start, unsigned int end) : Feature<T>("isoform", start, end), id(id) {
            this->length = 0;
        }

        virtual ~Isoform() {
        }

        struct FeatureComp {

            bool operator()(const SPtrFeature<T> l, const SPtrFeature<T> r) {
                return *l < *r;
            }
        };

        std::set<SPtrFeature<T>, FeatureComp>& getFeatures() {
            return features;
        }

        std::unordered_map<std::string, std::string>& getFields() {
            return fields;
        }

        std::string getId() const {
            return id;
        }

        unsigned int getLength() const {
            return this->length;
        }

        std::unordered_map<std::string, std::string> getFields() const {
            return fields;
        }

        void setFields(std::unordered_map<std::string, std::string> fields) {
            this->fields = fields;
        }

        void insertFeature(SPtrFeature<T> f) {
            std::pair <typename std::set<SPtrFeature<T>, FeatureComp>::iterator, bool> res = features.insert(f);
            if (res.second) {
                if (this->getStart() > f.get()->getStart()) {
                    this->setStart(f.get()->getStart());
                }
                if (this->getEnd() < f.get()->getEnd()) {
                    this->setEnd(f.get()->getEnd());
                }
                if (f->getType().compare("exon") == 0) {
                    length += f->getLength();
                }
            }
        }

        void insertFeaturesInGaps(std::string featToUse, std::string featToInsert) {
            SPtrFeature<T> f;
            for (auto it = features.begin(); it != features.end(); ++it) {
                if ((*it).get()->getType().compare(featToUse) == 0) {
                    for (auto it1 = it; it1 != features.end(); ++it1) {
                        if (*it != *it1) {
                            if ((*it1).get()->getType().compare(featToUse) != 0) {
                                if (!(*it).get()->isInside((*it1).get()->getStart(), (*it1).get()->getEnd(), 0)) break;
                            } else {
                                f = std::make_shared<Feature < T >> (featToInsert, (*it).get()->getEnd() + 1, (*it1).get()->getStart() - 1);
                                f.get()->setStrand((*it).get()->getStrand());
                                insertFeature(f);
                                break;
                            }
                        }
                    }
                }
            }
        }

        typename std::set<SPtrFeature<T>, FeatureComp>::iterator findFeatureUpperBound(unsigned int start, unsigned int end) {
            typename std::set < SPtrFeature<T>, typename Isoform<T>::FeatureComp>::iterator it;
            SPtrFeature<T> f = std::make_shared<Feature < T >> ("", start, end);
            it = features.lower_bound(f);
            if (it == features.end()) --it;
            return it;
        }

        /**
         * Process an array of words coming from a single line of a GTF file format
         * @param words vector of strings
         */
        void processGTFLine(std::vector<std::string> &words) {
            SPtrFeature<T> f = std::make_shared<Feature < T >> (words[2], atoi(words[3].c_str()) - 1, atoi(words[4].c_str()) - 1);
            f->setStrand(words[6][0]);
            insertFeature(f);
        }
    private:
        std::string id;
        std::unordered_map<std::string, std::string> fields;
        std::set<SPtrFeature<T>, FeatureComp> features;
        unsigned int length;
    };

    template<typename T>
    using FeatureSetItr = typename std::set<SPtrFeature<T>, typename Isoform<T>::FeatureComp>::iterator;

    template<typename T>
    using SPtrIsoform = std::shared_ptr<Isoform<T>>;

    template<typename T>
    using IsoformUnMap = std::unordered_map<std::string, SPtrIsoform<T>>;

    template<typename T>
    using IsoformUnMapItr = typename IsoformUnMap<T>::iterator;

    template<typename T>
    using IsoformMultiSet = std::multiset<SPtrIsoform<T>, typename Isoform<T>::FeatureComp>;

    template<typename T>
    using IsoformMultiSetItr = typename IsoformMultiSet<T>::iterator;

    template<typename T>
    class Gene : public Feature<T> {
    public:

        Gene(std::string id, unsigned int start, unsigned int end) : Feature<T>("gene", start, end), id(id) {
            this->currentIsoform = nullptr;
        }

        virtual ~Gene() {
        }

        std::string getId() const {
            return id;
        }

        SPtrIsoform<T> getCurrentIsoform() const {
            return currentIsoform;
        }

        void setCurrentIsoform(std::string isoformName) {
            if (!currentIsoform || currentIsoform->getId().compare(isoformName) != 0) {
                this->currentIsoform = findIsoform(isoformName);
            }
        }

        IsoformUnMap<T>& getIsoformsNameIndex() {
            return isoformsNameIndex;
        }

        IsoformMultiSet<T>& getIsoforms() {
            return isoforms;
        }

        std::set<SPtrFeature<T>, typename Isoform<T>::FeatureComp>& getUniquefeatures() {
            return uniquefeatures;
        }

        SPtrIsoform<T> findIsoform(std::string isoformName) {
            IsoformUnMapItr<T> it;
            it = isoformsNameIndex.find(isoformName);
            if (it == isoformsNameIndex.end()) {
                throw exceptions::NotFoundException("Isoform with name: " + isoformName + " does not exist");
            }
            return it->second;
        }

        IsoformMultiSetItr<T> findIsoformUpperBound(unsigned int start, unsigned int end) {
            IsoformMultiSetItr<T> it;
            std::shared_ptr<Isoform < T>> i = std::make_shared<Isoform < T >> ("isoform", start, end);
            it = isoforms.lower_bound(i);
            if (it == isoforms.end()) --it;
            return it;
        }

        /**
         * Process an array of words coming from a single line of a GTF file format
         * @param words vector of strings
         * @param isoformName name of the isoform
         * @param fieldsMap map of attributes created from a single line of a GTF file (words[9] splited by '"; ')
         */
        void processGTFLine(std::vector<std::string> &words, std::string isoformName, std::unordered_map<std::string, std::string>& fieldsMap) {
            IsoformUnMapItr<T> it;

            try {
                setCurrentIsoform(isoformName);
            } catch (exceptions::NotFoundException ex) {
                std::pair < IsoformUnMapItr<T>, bool> res = isoformsNameIndex.insert(std::make_pair(isoformName, std::make_shared<Isoform < T >> (isoformName, atoi(words[3].c_str()), atoi(words[4].c_str()))));
                if (!res.second) {
                    std::cerr << "Error inserting isoform" << std::endl;
                    exit(-1);
                }
                currentIsoform = res.first->second;
                currentIsoform->setStrand(words[6][0]);
            }
            currentIsoform->processGTFLine(words);
            currentIsoform->setFields(fieldsMap);
        }

        void createUniqueIsoformFeatures(int intronCutOff) {
            //std::cout << "createUniqueIsoformFeatures - START" << std::endl;
            for (IsoformMultiSetItr<T> iIt = isoforms.begin(); iIt != isoforms.end(); ++iIt) {
                SPtrIsoform<T> isoform = *iIt;
                for (FeatureSetItr<T> fIt = isoform->getFeatures().begin(); fIt != isoform->getFeatures().end(); ++fIt) {
                    SPtrFeature<T> feature = *fIt;
                    std::set<std::pair<unsigned int, unsigned int>> segments;
                    std::pair<unsigned int, unsigned int> seg = std::make_pair(feature->getStart(), feature->getEnd());
                    for (IsoformMultiSetItr<T> iIt2 = isoforms.begin(); iIt2 != isoforms.end(); ++iIt2) {
                        if (iIt != iIt2) {
                            for (FeatureSetItr<T> fIt2 = (*iIt2)->getFeatures().begin(); fIt2 != (*iIt2)->getFeatures().end(); ++fIt2) {
                                if (feature->getType() == "exon" && (*fIt2)->getType() == "exon") {
                                    if ((*fIt2)->getStart() <= seg.first && (*fIt2)->getEnd() > seg.first && (*fIt2)->getEnd() <= seg.second) {
                                        seg.first = (*fIt2)->getStart();
                                    }
                                    if ((*fIt2)->getStart() > seg.first && (*fIt2)->getStart() < seg.second && (*fIt2)->getEnd() >= seg.second) {
                                        seg.second = (*fIt2)->getEnd();
                                    }
                                }
                                if (feature->getType() == "intron" && (*fIt2)->getType() == "exon") {
                                    if ((*fIt2)->getStart() <= seg.first && (*fIt2)->getEnd() > seg.first && (*fIt2)->getEnd() < seg.second) {
                                        seg.first = (*fIt2)->getEnd();
                                    } else if ((*fIt2)->getStart() > seg.first && (*fIt2)->getStart() < seg.second && (*fIt2)->getEnd() >= seg.second) {
                                        seg.second = (*fIt2)->getStart();
                                    } else if ((*fIt2)->getStart() > seg.first && (*fIt2)->getStart() < seg.second && (*fIt2)->getEnd() < seg.second) {
                                        seg.second = (*fIt2)->getStart();
                                        segments.insert(seg);
                                        seg.first = (*fIt2)->getEnd();
                                        seg.second = feature->getEnd();
                                    } else if ((*fIt2)->getStart() <= seg.first && (*fIt2)->getEnd() >= seg.second) {
                                        seg.first = 0;
                                        seg.second = 0;
                                    }
                                }
                                if ((*fIt2)->getStart() > seg.second || (seg.first == 0 && seg.second == 0)) break;
                            }
                        }
                    }
                    if (seg.first != feature->getStart() || seg.second != feature->getEnd() || !segments.empty()) {
                        if (!segments.empty()) {
                            for (auto it = segments.begin(); it != segments.end(); ++it) {
                                if ((it->second - it->first + 1) > intronCutOff) {
                                    SPtrFeature<T> f = std::make_shared<Feature < T >> ("intron", it->first, it->second);
                                    f->setStrand(this->getStrand());
                                    uniquefeatures.insert(f);
                                }
                            }
                            if (seg.first != 0 && seg.second != 0) {
                                bool add = true;
                                if (feature->getType() == "intron" && (seg.second - seg.first + 1) < intronCutOff) add = false;
                                if (add) {
                                    SPtrFeature<T> f = std::make_shared<Feature < T >> (feature->getType(), seg.first, seg.second);
                                    f->setStrand(this->getStrand());
                                    uniquefeatures.insert(f);
                                }
                            }
                        } else if (seg.first != 0 && seg.second != 0) {
                            bool add = true;
                            if (feature->getType() == "intron" && (seg.second - seg.first + 1) < intronCutOff) add = false;
                            if (add) {
                                SPtrFeature<T> f = std::make_shared<Feature < T >> (feature->getType(), seg.first, seg.second);
                                f->setStrand(this->getStrand());
                                uniquefeatures.insert(f);
                            }
                        }
                    } else {
                        bool add = true;
                        if (feature->getType() == "intron" && feature->getLength() < intronCutOff) add = false;
                        if (add) {
                            SPtrFeature<T> f = std::make_shared<Feature < T >> (feature->getType(), feature->getStart(), feature->getEnd());
                            f->setStrand(this->getStrand());
                            uniquefeatures.insert(f);
                        }
                    }
                }
            }
            //std::cout << "createUniqueIsoformFeatures - END" << std::endl;
        }

        void createUniqueIntronFeatures(std::shared_ptr<Gene<T>> gene, int intronCutOff) {
            //std::cout << "createUniqueIntronFeatures - START" << std::endl;
            for (FeatureSetItr<T> fIt = uniquefeatures.begin(); fIt != uniquefeatures.end();) {
                SPtrFeature<T> feature = *fIt;
                if (feature->getType() == "intron") {
                    std::set<std::pair<unsigned int, unsigned int>> segments;
                    std::pair<unsigned int, unsigned int> seg = std::make_pair(feature->getStart(), feature->getEnd());
                    for (FeatureSetItr<T> fIt2 = gene->getUniquefeatures().begin(); fIt2 != gene->getUniquefeatures().end(); ++fIt2) {
                        if ((*fIt2)->getType() == "exon") {
                            if ((*fIt2)->getStart() <= seg.first && (*fIt2)->getEnd() > seg.first && (*fIt2)->getEnd() < seg.second) {
                                seg.first = (*fIt2)->getEnd();
                            } else if ((*fIt2)->getStart() > seg.first && (*fIt2)->getStart() < seg.second && (*fIt2)->getEnd() >= seg.second) {
                                seg.second = (*fIt2)->getStart();
                            } else if ((*fIt2)->getStart() > seg.first && (*fIt2)->getStart() < seg.second && (*fIt2)->getEnd() < seg.second) {
                                seg.second = (*fIt2)->getStart();
                                segments.insert(seg);
                                seg.first = (*fIt2)->getEnd();
                                seg.second = feature->getEnd();
                            } else if ((*fIt2)->getStart() <= seg.first && (*fIt2)->getEnd() >= seg.second) {
                                seg.first = 0;
                                seg.second = 0;
                            }
                            if ((*fIt2)->getStart() > seg.second || (seg.first == 0 && seg.second == 0)) break;
                        }
                    }
                    if (seg.first != feature->getStart() || seg.second != feature->getEnd() || !segments.empty()) {
                        if (!segments.empty()) {
                            fIt = uniquefeatures.erase(fIt);
                            for (auto it = segments.begin(); it != segments.end(); ++it) {
                                if ((it->second - it->first + 1) > intronCutOff) {
                                    SPtrFeature<T> f = std::make_shared<Feature < T >> ("intron", it->first, it->second);
                                    f->setStrand(this->getStrand());
                                    uniquefeatures.insert(f);
                                }
                            }
                            if (seg.first != 0 && seg.second != 0) {
                                if ((seg.second - seg.first + 1) > intronCutOff) {
                                    SPtrFeature<T> f = std::make_shared<Feature < T >> ("intron", seg.first, seg.second);
                                    f->setStrand(this->getStrand());
                                    uniquefeatures.insert(f);
                                }
                            }
                        } else if (seg.first != 0 && seg.second != 0) {
                            fIt = uniquefeatures.erase(fIt);
                            if ((seg.second - seg.first + 1) > intronCutOff) {
                                SPtrFeature<T> f = std::make_shared<Feature < T >> ("intron", seg.first, seg.second);
                                f->setStrand(this->getStrand());
                                uniquefeatures.insert(f);
                            }
                        } else {
                            fIt = uniquefeatures.erase(fIt);
                        }
                    } else {
                        if (feature->getLength() < intronCutOff) {
                            fIt = uniquefeatures.erase(fIt);
                        } else {
                            ++fIt;
                        }
                    }
                } else {
                    ++fIt;
                }
            }
            //std::cout << "createUniqueIntronFeatures - END" << std::endl;
        }

        void printGeneUniqueIsoformFeaturesGTF(std::ofstream& gtfFile, std::string chr, bool include_introns) {
            SPtrFeature<T> f;
            for (auto itr = uniquefeatures.begin(); itr != uniquefeatures.end(); ++itr) {
                f = *itr;
                bool print = true;
                if (!include_introns && f->getType() == "intron") print = false;
                if (print) {
                    gtfFile << chr << "\tunknown\t" << f->getType() << "\t"
                            << (f->getStart() + 1) << "\t"
                            << (f->getEnd() + 1) << "\t.\t"
                            << f->getStrand() << "\t.\tgene_id \""
                            << this->getId() << "\"; gene_name \""
                            << this->getId() << "\"; transcript_id \""
                            << this->getId() << "\""
                            << std::endl;
                }
            }
        }
    private:
        std::string id;
        IsoformMultiSet<T> isoforms;
        IsoformUnMap<T> isoformsNameIndex;
        SPtrIsoform<T> currentIsoform;
        std::set<SPtrFeature<T>, typename Isoform<T>::FeatureComp> uniquefeatures;
    };

    template<typename T>
    using SPtrGene = std::shared_ptr<Gene<T>>;

    template<typename T>
    using GeneUnMap = std::unordered_map<std::string, SPtrGene<T>>;

    template<typename T>
    using GeneUnMapItr = typename GeneUnMap<T>::iterator;

    template<typename T>
    using GeneMultiSet = std::multiset<SPtrGene<T>, typename Isoform<T>::FeatureComp>;

    template<typename T>
    using GeneMultiSetItr = typename GeneMultiSet<T>::iterator;

    template<typename T>
    class Chromosome {
    public:

        Chromosome(std::string id) : id(id) {
            this->currentGene = nullptr;
        }

        virtual~Chromosome() {
        }

        GeneMultiSet<T>& getGenes() {
            return genes;
        }

        GeneUnMap<T>& getGenesNameIndex() {
            return genesNameIndex;
        }

        std::string getId() const {
            return id;
        }

        SPtrGene<T> getCurrentGene() const {
            return currentGene;
        }

        void setCurrentGene(std::string geneName) {
            if (!currentGene || currentGene->getId().compare(geneName) != 0) {
                this->currentGene = findGene(geneName);
            }
        }

        SPtrGene<T> findGene(std::string geneName) {
            GeneUnMapItr<T> it;
            it = genesNameIndex.find(geneName);
            if (it == genesNameIndex.end()) {
                throw exceptions::NotFoundException("Gene with name: " + geneName + " does not exist");
            }
            return it->second;
        }

        GeneMultiSetItr<T> findGeneUpperBound(unsigned int start, unsigned int end) {
            GeneMultiSetItr<T> it;
            std::shared_ptr<Gene < T>> g = std::make_shared<Gene < T >> ("gene", start, end);
            it = genes.lower_bound(g);
            if (it == genes.end()) --it;
            return it;
        }

        void createUniqueIntronsPerGene(int intronCutOff) {
            for (GeneMultiSetItr<T> gIt = genes.begin(); gIt != genes.end(); ++gIt) {
                SPtrGene<T> g = *gIt;
                GeneMultiSetItr<T> overlapDown = genes.end();
                GeneMultiSetItr<T> overlapUp = genes.end();
                for (auto it = gIt; it != genes.begin(); --it) {
                    if (it != gIt) {
                        if ((*it)->getEnd() < g->getStart()) break;
                        overlapDown = it;
                    }
                }
                for (auto it = gIt; it != genes.end(); ++it) {
                    if (it != gIt) {
                        if (g->getEnd() < (*it)->getStart()) break;
                        overlapUp = it;
                    }
                }
                if (overlapDown != genes.end() && overlapUp != genes.end()) {
                    ++overlapUp;
                    for (auto it = overlapDown; it != overlapUp; ++it) {
                        if (it != gIt) {
                            g->createUniqueIntronFeatures(*it, intronCutOff);
                        }
                    }
                } else if (overlapDown != genes.end()) {
                    for (auto it = overlapDown; it != gIt; ++it) {
                        g->createUniqueIntronFeatures(*it, intronCutOff);
                    }
                } else if (overlapUp != genes.end()) {
                    ++overlapUp;
                    for (auto it = gIt; it != overlapUp; ++it) {
                        if (it != gIt) {
                            g->createUniqueIntronFeatures(*it, intronCutOff);
                        }
                    }
                }
            }
        }

        /**
         * Process an array of words coming from a single line of a GTF file format
         * @param words vector of strings
         * @param geneIdKey string to be used to identify the genes id
         * @param isoformIdKey string to be used to identify the isoforms id
         */
        void processGTFLine(std::vector<std::string> &words, std::string geneIdKey, std::string isoformIdKey) {
            std::vector<std::string> fields;
            std::string geneName;
            std::string isoformName;
            unsigned int wStart, wEnd;
            std::unordered_map<std::string, std::string> fieldsMap;
            GeneUnMapItr<T> it;

            wStart = atoi(words[3].c_str()) - 1;
            wEnd = atoi(words[4].c_str()) - 1;
            BString::split(words[8], ";", fields);
            geneName = "";
            isoformName = "";
            for (size_t i = 0; i < fields.size(); i++) {
                std::vector<std::string> fieldsIn;
                BString::split(fields[i], "\"", fieldsIn);
                if (fieldsIn.size() == 2) {
                    std::string key = BString::trim(fieldsIn[0]);
                    std::string value = BString::trim(fieldsIn[1]);
                    if (geneIdKey.compare(key) == 0) {
                        geneName = value;
                    } else if (isoformIdKey.compare(key) == 0) {
                        isoformName = value;
                    } else {
                        fieldsMap.insert(std::make_pair(key, value));
                    }
                }

            }
            if (geneName.empty()) {
                throw exceptions::NotFoundException("Key " + geneIdKey + " for gene name was not found on GTF line.");
            }
            if (isoformName.empty()) {
                throw exceptions::NotFoundException("Key " + isoformName + " for isoform name was not found on GTF line.");
            }
            try {
                setCurrentGene(geneName);
            } catch (exceptions::NotFoundException ex) {
                std::pair < GeneUnMapItr<T>, bool> res = genesNameIndex.insert(make_pair(geneName, std::make_shared<Gene < T >> (geneName, wStart, wEnd)));
                if (!res.second) {
                    std::cerr << "Error inserting gene" << std::endl;
                    exit(-1);
                }
                currentGene = res.first->second;
                currentGene->setStrand(words[6][0]);
            }
            if (currentGene->getStart() > wStart) {
                currentGene->setStart(wStart);
            }
            if (currentGene->getEnd() < wEnd) {

                currentGene->setEnd(wEnd);
            }
            currentGene->processGTFLine(words, isoformName, fieldsMap);
        }
    private:
        std::string id;
        GeneMultiSet<T> genes;
        GeneUnMap<T> genesNameIndex;
        SPtrGene<T> currentGene;
    };

    template<typename T>
    using SPtrChromosome = std::shared_ptr<Chromosome<T>>;

    template<typename T>
    using ChromosomeUnMap = std::unordered_map<std::string, SPtrChromosome<T>>;

    template<typename T>
    using ChromosomeUnMapItr = typename ChromosomeUnMap<T>::iterator;

    template<typename T>
    using GeneIsoformUnMap = std::unordered_map<std::string, std::pair<SPtrGene<T>, SPtrIsoform<T>>>;

    template<typename T>
    using GeneIsoformUnMapItr = typename GeneIsoformUnMap<T>::iterator;

    template<typename T>
    class GenomeFactory {
    public:

        GenomeFactory() {
            intronCutOff = 16;
            currentChr = nullptr;
            currentIso = nullptr;
            currentGene = nullptr;
        }

        GenomeFactory(int intronCutOff) {
            this->intronCutOff = intronCutOff;
            currentChr = nullptr;
            currentIso = nullptr;
            currentGene = nullptr;
        }

        virtual ~GenomeFactory() {
        }

        ChromosomeUnMap<T>& getChromosomes() {

            return chromosomes;
        }

        unsigned int size() {

            return chromosomes.size();
        }

        int getIntronCutOff() {
            return intronCutOff;
        }

        void setIntronCutOff(int intronCutOff) {
            this->intronCutOff = intronCutOff;
        }

        Chromosome<T>* getCurrentChr() {

            return currentChr;
        }

        void setCurrentChr(std::string chrName) {
            if (!currentChr || currentChr->getId().compare(chrName) != 0) {

                this->currentChr = findChromosome(chrName);
            }
        }

        Chromosome<T> *findChromosome(std::string chrName) {
            ChromosomeUnMapItr<T> it;
            it = chromosomes.find(chrName);
            if (it == chromosomes.end()) {

                throw exceptions::NotFoundException("Chromosome with name: " + chrName + " does not exist");
            }
            return it->second.get();
        }

        SPtrGene<T> findGene(std::string chrName, std::string geneName) {
            if (!currentGene || currentGene->getId().compare(geneName) != 0) {

                setCurrentChr(chrName);
                currentGene = currentChr->findGene(geneName);
            }
            return currentGene;
        }

        SPtrGene<T> getCurrentGene() const {

            return currentGene;
        }

        SPtrIsoform<T> getCurrentIso() const {

            return currentIso;
        }

        GeneMultiSetItr<T> findGeneUpperBound(std::string chrName, unsigned int start, unsigned int end) {
            GeneMultiSetItr<T> it;
            setCurrentChr(chrName);
            SPtrGene<T> g = std::make_shared<Gene < T >> ("gene", start, end);
            it = currentChr->getGenes().lower_bound(g);
            if (it == currentChr->getGenes().end()) --it;

            return it;
        }

        SPtrIsoform<T> findIsoform(std::string isoName) {
            if (!currentIso || currentIso->getId().compare(isoName) != 0) {
                GeneIsoformUnMapItr<T> transIt = transcript2Chr.find(isoName);
                if (transIt == transcript2Chr.end()) {

                    throw exceptions::NotFoundException("Isoform with name: " + isoName + " does not exist");
                }
                currentGene = transIt->second.first;
                currentIso = transIt->second.second;
            }
            return currentIso;
        }

        /**
         * Process a GTF file and creates a genome data structure in RAM
         * @param gtfFileName GTF file name
         * @param geneIdKey string to be used to identify the genes id
         * @param isoformIdKey string to be used to identify the isoforms id
         * @param features a set of features to include in the genome
         * @param featuresToCreate unordered_map of features to create. Key is the features parsed 
         * from GTF and the value is the feature to create between them 
         */
        void processGTFFile(std::string gtfFileName, std::string geneIdKey, std::string isoformIdKey, std::set<std::string>& features, std::unordered_map<std::string, std::string>& featuresToCreate) {
            SPtrChromosome<T> c;
            SPtrGene<T> g;
            SPtrIsoform<T> i;

            try {
                parsers::TextParser fParser;
                fParser.setFileToParse(gtfFileName);
                while (fParser.iterate("#", "\t")) {
                    if (fParser.getWords().size() != 9) {
                        std::cerr << "GTF file with wrong number of fields. It should be 9 tab separated fields" << std::endl;
                        exit(-1);
                    }
                    if (features.find(fParser.getWords()[2]) != features.end()) {
                        try {
                            setCurrentChr(fParser.getWords()[0]);
                        } catch (exceptions::NotFoundException ex) {
                            std::pair < ChromosomeUnMapItr<T>, bool> res = chromosomes.insert(std::make_pair(fParser.getWords()[0], std::make_shared<Chromosome < T >> (fParser.getWords()[0])));
                            if (!res.second) {
                                std::cerr << "Error inserting new Chromosome" << std::endl;
                                exit(-1);
                            }
                            currentChr = res.first->second.get();
                        }
                        try {
                            currentChr->processGTFLine(fParser.getWords(), geneIdKey, isoformIdKey);
                        } catch (exceptions::NotFoundException ex) {
                            std::cerr << fParser.getLine() << std::endl;
                            exit(-1);
                        }
                    }
                }
            } catch (exceptions::FileHandledException ex) {
                std::cerr << "Error parsing file: " << gtfFileName << std::endl;
                exit(-1);
            }

            for (auto it = chromosomes.begin(); it != chromosomes.end(); ++it) {
                c = it->second;
                for (auto it1 = c->getGenesNameIndex().begin(); it1 != c->getGenesNameIndex().end(); ++it1) {
                    g = it1->second;
                    for (auto it2 = g->getIsoformsNameIndex().begin(); it2 != g->getIsoformsNameIndex().end(); ++it2) {
                        i = it2->second;
                        for (auto it3 = featuresToCreate.begin(); it3 != featuresToCreate.end(); ++it3) {
                            i->insertFeaturesInGaps(it3->first, it3->second);
                        }
                        transcript2Chr.insert(make_pair(i->getId(), std::make_pair(g, i)));
                        g->getIsoforms().insert(i);
                    }
                    g->createUniqueIsoformFeatures(intronCutOff);
                    c->getGenes().insert(g);
                }
                c->createUniqueIntronsPerGene(intronCutOff);
            }
        }
    private:
        ChromosomeUnMap<T> chromosomes;
        GeneIsoformUnMap<T> transcript2Chr;
        SPtrIsoform<T> currentIso;
        SPtrGene<T> currentGene;
        Chromosome<T> *currentChr;
        int intronCutOff;
    };

}


#endif /* GENOMEFACTORY_H */

