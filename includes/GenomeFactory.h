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
#include <vector>
#include <string>

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
            if (this->start != right.getStart())
                return this->start > right.getStart();
            else if (this->end != right.getEnd())
                return this->end > right.getEnd();
            return false;
        }

        bool operator<(const Feature& right) const {
            return right > * this;
        }

        friend std::ostream& operator<<(std::ostream& os, Feature *obj) {
            os << obj->getType() << " " << obj->getStrand();
            os << " coord: " << obj->getStart() << "-" << obj->getEnd() << "\tlength: " << obj->getLength();
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

        friend std::ostream& operator<<(std::ostream& os, std::shared_ptr<Isoform<T>> obj) {
            os << obj->getId() << " Length: " << obj->getLength() << " " << obj->getStart() << "-" << obj->getEnd() << " Features: " << obj->getFeatures().size() << std::endl;
            for (auto f : obj->getFeatures()) {
                os << "\t\t\t" << f;
                os << std::endl;
            }
            return os;
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
            this->length = 0;
        }

        virtual ~Gene() {
        }

        std::string getId() const {
            return id;
        }

        unsigned int getLength() const {
            return this->length;
        }

        void setLength(unsigned int length) {
            this->length = length;
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

        void setIsoforms(IsoformMultiSet<T> isoforms) {
            this->isoforms = isoforms;
        }

        std::set<SPtrFeature<T>, typename Isoform<T>::FeatureComp>& getFeatures() {
            return features;
        }

        std::set<SPtrFeature<T>, typename Isoform<T>::FeatureComp>& getUniqueFeatures() {
            return uniqueFeatures;
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

        void createGeneFeatures(unsigned int intronCutOff) {
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
                                    if (((*fIt2)->getEnd() > seg.first && (*fIt2)->getEnd() <= seg.second) ||
                                            ((*fIt2)->getStart() >= seg.first && (*fIt2)->getStart() < seg.second) ||
                                            ((*fIt2)->getStart() < seg.first && (*fIt2)->getEnd() > seg.second)) {
                                        seg.first = std::min(seg.first, (*fIt2)->getStart());
                                        seg.second = std::max(seg.second, (*fIt2)->getEnd());
                                    }
                                }
                                if (feature->getType() == "intron" && (*fIt2)->getType() == "exon") {
                                    if ((*fIt2)->getStart() <= seg.first && (*fIt2)->getEnd() > seg.first && (*fIt2)->getEnd() < seg.second) {
                                        seg.first = (*fIt2)->getEnd() + 1;
                                    } else if ((*fIt2)->getStart() > seg.first && (*fIt2)->getStart() < seg.second && (*fIt2)->getEnd() >= seg.second) {
                                        seg.second = (*fIt2)->getStart() - 1;
                                    } else if ((*fIt2)->getStart() > seg.first && (*fIt2)->getStart() < seg.second && (*fIt2)->getEnd() < seg.second) {
                                        seg.second = (*fIt2)->getStart() - 1;
                                        segments.insert(seg);
                                        seg.first = (*fIt2)->getEnd() + 1;
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
                                std::pair<unsigned int, unsigned int> seg_in = std::make_pair(it->first, it->second);
                                for (IsoformMultiSetItr<T> iIt2 = isoforms.begin(); iIt2 != isoforms.end(); ++iIt2) {
                                    for (FeatureSetItr<T> fIt2 = (*iIt2)->getFeatures().begin(); fIt2 != (*iIt2)->getFeatures().end(); ++fIt2) {
                                        if ((*fIt2)->getType() == "exon") {
                                            if ((*fIt2)->getStart() <= seg_in.first && (*fIt2)->getEnd() > seg_in.first && (*fIt2)->getEnd() < seg_in.second) {
                                                seg_in.first = (*fIt2)->getEnd() + 1;
                                            } else if ((*fIt2)->getStart() > seg_in.first && (*fIt2)->getStart() < seg_in.second && (*fIt2)->getEnd() >= seg_in.second) {
                                                seg_in.second = (*fIt2)->getStart() - 1;
                                            } else if ((*fIt2)->getStart() <= seg_in.first && (*fIt2)->getEnd() >= seg_in.second) {
                                                seg_in.first = 0;
                                                seg_in.second = 0;
                                            }
                                        }
                                        if ((*fIt2)->getStart() > seg_in.second || (seg_in.first == 0 && seg_in.second == 0)) break;
                                    }
                                }
                                if (seg_in.first != 0 && seg_in.second != 0 && (seg_in.second - seg_in.first + 1) > intronCutOff) {
                                    SPtrFeature<T> f = std::make_shared<Feature < T >> ("intron", seg_in.first, seg_in.second);
                                    f->setStrand(this->getStrand());
                                    features.insert(f);
                                }
                            }
                            if (seg.first != 0 && seg.second != 0) {
                                bool add = true;
                                if (feature->getType() == "intron" && (seg.second - seg.first + 1) < intronCutOff) add = false;
                                if (add) {
                                    SPtrFeature<T> f = std::make_shared<Feature < T >> (feature->getType(), seg.first, seg.second);
                                    f->setStrand(this->getStrand());
                                    features.insert(f);
                                }
                            }
                        } else if (seg.first != 0 && seg.second != 0) {
                            bool add = true;
                            if (feature->getType() == "intron" && (seg.second - seg.first + 1) < intronCutOff) add = false;
                            if (add) {
                                SPtrFeature<T> f = std::make_shared<Feature < T >> (feature->getType(), seg.first, seg.second);
                                f->setStrand(this->getStrand());
                                features.insert(f);
                            }
                        }
                    } else {
                        bool add = true;
                        if (feature->getType() == "intron" && feature->getLength() < intronCutOff) add = false;
                        if (add) {
                            SPtrFeature<T> f = std::make_shared<Feature < T >> (feature->getType(), feature->getStart(), feature->getEnd());
                            f->setStrand(this->getStrand());
                            features.insert(f);
                        }
                    }
                }
            }

            for (FeatureSetItr<T> fIt = this->features.begin(); fIt != this->features.end();) {
                SPtrFeature<T> feature = *fIt;
                std::pair<unsigned int, unsigned int> seg = std::make_pair(feature->getStart(), feature->getEnd());
                for (FeatureSetItr<T> fIt2 = this->features.begin(); fIt2 != this->features.end(); ++fIt2) {
                    if (fIt2 != fIt) {
                        if (((*fIt2)->getEnd() > seg.first && (*fIt2)->getEnd() <= seg.second) ||
                                ((*fIt2)->getStart() >= seg.first && (*fIt2)->getStart() < seg.second) ||
                                ((*fIt2)->getStart() < seg.first && (*fIt2)->getEnd() > seg.second)) {
                            seg.first = std::min(seg.first, (*fIt2)->getStart());
                            seg.second = std::max(seg.second, (*fIt2)->getEnd());
                        }
                    }
                }
                if (seg.first != feature->getStart() || seg.second != feature->getEnd()) {
                    fIt = this->features.erase(fIt);
                    SPtrFeature<T> f = std::make_shared<Feature < T >> (feature->getType(), seg.first, seg.second);
                    f->setStrand(this->getStrand());
                    this->features.insert(f);
                } else {
                    ++fIt;
                }
            }
            this->length = 0;
            for (FeatureSetItr<T> fIt = this->features.begin(); fIt != this->features.end(); ++fIt) {
                if ((*fIt)->getType() == "exon") this->length += (*fIt)->getLength();
            }
        }

        void createUniqueIntronFeatures(std::shared_ptr<Gene<T>> gene, unsigned int intronCutOff) {
            std::cout << "createUniqueIntronFeatures - START: " << this->id << " compared with " << gene->getId() << std::endl;
            for (FeatureSetItr<T> fIt = features.begin(); fIt != features.end();) {
                SPtrFeature<T> feature = *fIt;
                if (feature->getType() == "intron") {
                    std::set<std::pair<unsigned int, unsigned int>> segments;
                    std::pair<unsigned int, unsigned int> seg = std::make_pair(feature->getStart(), feature->getEnd());
                    for (FeatureSetItr<T> fIt2 = gene->getFeatures().begin(); fIt2 != gene->getFeatures().end(); ++fIt2) {
                        if ((*fIt2)->getType() == "exon") {
                            if ((*fIt2)->getStart() <= seg.first && (*fIt2)->getEnd() > seg.first && (*fIt2)->getEnd() < seg.second) {
                                seg.first = (*fIt2)->getEnd() + 1;
                            } else if ((*fIt2)->getStart() >= seg.first && (*fIt2)->getStart() < seg.second && (*fIt2)->getEnd() >= seg.second) {
                                seg.second = (*fIt2)->getStart() - 1;
                            } else if ((*fIt2)->getStart() > seg.first && (*fIt2)->getStart() < seg.second && (*fIt2)->getEnd() < seg.second) {
                                seg.second = (*fIt2)->getStart() - 1;
                                segments.insert(seg);
                                seg.first = (*fIt2)->getEnd() + 1;
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
                            fIt = features.erase(fIt);
                            for (auto it = segments.begin(); it != segments.end(); ++it) {
                                std::pair<unsigned int, unsigned int> seg_in = std::make_pair(it->first, it->second);
                                for (FeatureSetItr<T> fIt2 = gene->getFeatures().begin(); fIt2 != gene->getFeatures().end(); ++fIt2) {
                                    if ((*fIt2)->getType() == "exon") {
                                        if ((*fIt2)->getStart() <= seg_in.first && (*fIt2)->getEnd() > seg_in.first && (*fIt2)->getEnd() < seg_in.second) {
                                            seg_in.first = (*fIt2)->getEnd() + 1;
                                        } else if ((*fIt2)->getStart() >= seg_in.first && (*fIt2)->getStart() < seg_in.second && (*fIt2)->getEnd() >= seg_in.second) {
                                            seg_in.second = (*fIt2)->getStart() - 1;
                                        } else if ((*fIt2)->getStart() <= seg_in.first && (*fIt2)->getEnd() >= seg_in.second) {
                                            seg_in.first = 0;
                                            seg_in.second = 0;
                                        }
                                        if ((*fIt2)->getStart() > seg_in.second || (seg_in.first == 0 && seg_in.second == 0)) break;
                                    }
                                }
                                if (seg_in.first != 0 && seg_in.second != 0 && (seg_in.second - seg_in.first + 1) > intronCutOff) {
                                    SPtrFeature<T> f = std::make_shared<Feature < T >> ("intron", seg_in.first, seg_in.second);
                                    f->setStrand(this->getStrand());
                                    features.insert(f);
                                }
                            }
                            if (seg.first != 0 && seg.second != 0) {
                                if ((seg.second - seg.first + 1) > intronCutOff) {
                                    SPtrFeature<T> f = std::make_shared<Feature < T >> ("intron", seg.first, seg.second);
                                    f->setStrand(this->getStrand());
                                    features.insert(f);
                                }
                            }
                        } else if (seg.first != 0 && seg.second != 0) {
                            fIt = features.erase(fIt);
                            if ((seg.second - seg.first + 1) > intronCutOff) {
                                SPtrFeature<T> f = std::make_shared<Feature < T >> ("intron", seg.first, seg.second);
                                f->setStrand(this->getStrand());
                                features.insert(f);
                            }
                        } else {
                            fIt = features.erase(fIt);
                        }
                    } else {
                        if (feature->getLength() < intronCutOff) {
                            fIt = features.erase(fIt);
                        } else {
                            ++fIt;
                        }
                    }
                } else {
                    ++fIt;
                }
            }
        }

        void printGeneFeaturesGTF(std::ofstream& gtfFile, std::string chr, bool include_introns) {
            SPtrFeature<T> f;
            for (auto itr = features.begin(); itr != features.end(); ++itr) {
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

        void printGeneUniqueFeaturesGTF(std::ofstream& gtfFile, std::string chr, bool include_introns) {
            SPtrFeature<T> f;
            for (auto itr = uniqueFeatures.begin(); itr != uniqueFeatures.end(); ++itr) {
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

        friend std::ostream& operator<<(std::ostream& os, std::shared_ptr<Gene<T>> obj) {
            os << obj->getId() << " Length: " << obj->getLength() << " " << obj->getStart() << "-" << obj->getEnd();
            os << " Isoforms: " << obj->getIsoforms().size();
            os << " Unique: " << obj->getFeatures().size();
            os << std::endl;
            os << "\t\tIsoforms: " << std::endl;
            for (auto i : obj->getIsoforms()) {
                os << "\t\t" << i;
                os << std::endl;
            }
            os << "\t\tUnique: " << std::endl;
            for (auto i : obj->getFeatures()) {
                os << "\t\t" << i;
                os << std::endl;
            }
            return os;
        }
    private:
        std::string id;
        unsigned int length;
        IsoformMultiSet<T> isoforms;
        IsoformUnMap<T> isoformsNameIndex;
        SPtrIsoform<T> currentIsoform;
        std::set<SPtrFeature<T>, typename Isoform<T>::FeatureComp> features;
        std::set<SPtrFeature<T>, typename Isoform<T>::FeatureComp> uniqueFeatures;
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
    using GeneMultiSetRevItr = typename std::multiset<SPtrGene<T>, typename Isoform<T>::FeatureComp>::reverse_iterator;

    template<typename T>
    using GeneIsoformReadUnMap = typename std::unordered_map<std::string, std::shared_ptr<IsoformUnMap<T>>>;

    template<typename T>
    using GeneIsoformUnMap = std::unordered_map<std::string, std::pair<SPtrGene<T>, SPtrIsoform<T>>>;

    template<typename T>
    using GeneIsoformUnMapItr = typename GeneIsoformUnMap<T>::iterator;

    template<typename T>
    SPtrFeature<T> feature_overlap(SPtrFeature<T> f1, SPtrFeature<T> f2) {
        if (static_cast<int> (f1->getEnd()) - static_cast<int> (f2->getStart()) >= 0
                && static_cast<int> (f2->getEnd()) - static_cast<int> (f1->getStart()) >= 0) {
            unsigned int start = std::max(f1->getStart(), f2->getStart());
            unsigned int end = std::min(f1->getEnd(), f2->getEnd());
            SPtrFeature<T> f = std::make_shared<Feature < T >> (f1->getType(), start, end);
            f->setStrand(f1->getStrand());
            return f;
        }
        return nullptr;
    }

    template<typename T>
    SPtrFeature<T> feature_overlap(SPtrFeature<T> f1, SPtrGene<T> f2) {
        if (static_cast<int> (f1->getEnd()) - static_cast<int> (f2->getStart()) >= 0
                && static_cast<int> (f2->getEnd()) - static_cast<int> (f1->getStart()) >= 0) {
            unsigned int start = std::max(f1->getStart(), f2->getStart());
            unsigned int end = std::min(f1->getEnd(), f2->getEnd());
            SPtrFeature<T> f = std::make_shared<Feature < T >> (f1->getType(), start, end);
            f->setStrand(f1->getStrand());
            return f;
        }
        return nullptr;
    }

    template<typename T>
    SPtrFeature<T> feature_union(SPtrFeature<T> f1, SPtrFeature<T> f2) {
        if (static_cast<int> (f1->getEnd()) - static_cast<int> (f2->getStart()) >= 0
                && static_cast<int> (f2->getEnd()) - static_cast<int> (f1->getEnd()) >= 0) {
            unsigned int start = std::min(f1->getStart(), f2->getStart());
            unsigned int end = std::max(f1->getEnd(), f2->getEnd());
            SPtrFeature<T> f = std::make_shared<Feature < T >> (f1->getType(), start, end);
            f->setStrand(f1->getStrand());
            return f;
        }
        return nullptr;
    }

    template<typename T>
    std::set<SPtrFeature<T>, typename Isoform<T>::FeatureComp> feature_disunion(SPtrFeature<T> f1, SPtrGene<T> f2) {
        std::set < SPtrFeature<T>, typename Isoform<T>::FeatureComp> features;
        SPtrFeature<T> fi;
        SPtrFeature<T> f = feature_overlap(f1, f2);
        if (f && !(f->getStart() == f1->getStart() && f->getEnd() == f1->getEnd())) {
            if (f1->getEnd() == f->getEnd()) {
                fi = std::make_shared<Feature < T >> (f1->getType(), f1->getStart(), f->getStart() - 1);
                fi->setStrand(f1->getStrand());
                features.insert(fi);
            } else if (f1->getStart() == f->getStart()) {
                fi = std::make_shared<Feature < T >> (f1->getType(), f->getEnd() + 1, f1->getEnd());
                fi->setStrand(f1->getStrand());
                features.insert(fi);
            } else {
                fi = std::make_shared<Feature < T >> (f1->getType(), f1->getStart(), f->getStart() - 1);
                fi->setStrand(f1->getStrand());
                features.insert(fi);

                fi = std::make_shared<Feature < T >> (f1->getType(), f->getEnd() + 1, f1->getEnd());
                fi->setStrand(f1->getStrand());
                features.insert(fi);
            }
        }
        return features;
    }

    template<typename T>
    std::set<SPtrFeature<T>, typename Isoform<T>::FeatureComp> features_union(std::set < SPtrFeature<T>, typename Isoform<T>::FeatureComp> segments,
            std::set < SPtrFeature<T>, typename Isoform<T>::FeatureComp> feats) {
        std::set < SPtrFeature<T>, typename Isoform<T>::FeatureComp> features;

        if (feats.empty()) return segments;

        for (auto sIt = segments.begin(); sIt != segments.end(); ++sIt) {
            for (auto s2It = feats.begin(); s2It != feats.end(); ++s2It) {
                SPtrFeature<T> f = feature_overlap((*sIt), (*s2It));
                if (f) {
                    features.insert(f);
                }
            }
        }

        return features;
    }

    template<typename T>
    class Chromosome {
    public:

        Chromosome(std::string id) : id(id) {
            this->currentGene = nullptr;
            this->currentIsoform = nullptr;
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

        GeneIsoformReadUnMap<T> &getGeneIsoforms() {
            return geneIsoforms;
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

        void createGeneUniqueFeatures() {
            for (GeneMultiSetItr<T> gIt = genes.begin(); gIt != genes.end(); ++gIt) {
                SPtrGene<T> g = *gIt;
                for (FeatureSetItr<T> fIt = g->getFeatures().begin(); fIt != g->getFeatures().end(); ++fIt) {
                    SPtrFeature<T> f = *fIt;
                    std::set < SPtrFeature<T>, typename Isoform<T>::FeatureComp> segments;
                    bool overlapped = false;
                    for (GeneMultiSetItr<T> gIt2 = genes.begin(); gIt2 != genes.end(); ++gIt2) {
                        if (gIt != gIt2) {
                            SPtrGene<T> g2 = *gIt2;
                            if (g2->getStart() <= f->getStart() && g2->getEnd() >= f->getEnd()) {
                                segments.clear();
                                overlapped = true;
                                break;
                            }
                            std::set < SPtrFeature<T>, typename Isoform<T>::FeatureComp> feats = feature_disunion(f, g2);
                            
                            if (segments.empty()) {
                                segments.insert(feats.begin(), feats.end());
                            } else {
                                segments = features_union(segments, feats);
                                if (segments.empty()) {
                                    overlapped = true;
                                    break;
                                }

                            }
                        }
                    }
                    if (!overlapped) {
                        if (!segments.empty()) {
                            g->getUniqueFeatures().insert(segments.begin(), segments.end());                            
                        } else {
                            g->getUniqueFeatures().insert(f);
                        }
                    }
                }
            }
        }

        void createGeneMultiSetItr(int intronCutOff, std::unordered_map<std::string, std::string>& featuresToCreate, GeneIsoformUnMap<T> transcript2Chr) {
            genes.clear();
            for (auto gIt : geneIsoforms) {
                std::string geneName = gIt.first;
                std::shared_ptr<IsoformUnMap < T>> iMap = gIt.second;
                if (iMap->size() == 1) {
                    SPtrIsoform<T> i = iMap->begin()->second;
                    for (auto feat : featuresToCreate) {
                        i->insertFeaturesInGaps(feat.first, feat.second);
                    }
                    currentGene = std::make_shared<Gene < T >> (geneName, i->getStart(), i->getEnd());
                    currentGene->setStrand(i->getStrand());                    
                    currentGene->getIsoforms().insert(i);                    
                    unsigned int length = 0;
                    for (auto fIt : i->getFeatures()) {
                        if ((*fIt).getType().compare("exon") == 0) length += (*fIt).getLength();
                        currentGene->getFeatures().insert(fIt);
                    }
                    currentGene->setLength(length);

                    std::pair < GeneUnMapItr<T>, bool> res = genesNameIndex.insert(make_pair(geneName, currentGene));
                    if (!res.second) {
                        std::cerr << "Error inserting gene with single transcript" << std::endl;
                        exit(-1);
                    }
                    transcript2Chr.insert(make_pair(i->getId(), std::make_pair(currentGene, i)));
                    genes.insert(currentGene);
                } else {
                    std::map < std::pair<unsigned int, unsigned int>, IsoformMultiSet < T>> segments;
                    for (auto iIt : *(iMap)) {
                        std::string isoformName = iIt.first;
                        SPtrIsoform<T> i = iIt.second;
                        for (auto feat : featuresToCreate) {
                            i->insertFeaturesInGaps(feat.first, feat.second);
                        }
                        if (segments.empty()) {
                            IsoformMultiSet<T> iSet;
                            iSet.insert(i);
                            segments.insert(std::make_pair(std::make_pair(i->getStart(), i->getEnd()), iSet));
                        } else {
                            bool inserted = false;
                            for (auto sIt = segments.begin(); sIt != segments.end(); ++sIt) {
                                std::pair<unsigned int, unsigned int> seg = sIt->first;
                                IsoformMultiSet<T> iSet = sIt->second;
                                if ((i->getEnd() > seg.first && i->getEnd() <= seg.second) ||
                                        (i->getStart() >= seg.first && i->getStart() < seg.second) ||
                                        (i->getStart() < seg.first && i->getEnd() > seg.second)) {
                                    iSet.insert(i);
                                    seg.first = std::min(seg.first, i->getStart());
                                    seg.second = std::max(seg.second, i->getEnd());
                                    inserted = true;
                                    segments.erase(sIt);
                                    segments.insert(std::make_pair(seg, iSet));
                                    break;
                                }
                            }
                            if (!inserted) {
                                IsoformMultiSet<T> iSet;
                                iSet.insert(i);
                                segments.insert(std::make_pair(std::make_pair(i->getStart(), i->getEnd()), iSet));
                            }
                        }
                    }
                    int gene_copy = 1;
                    for (auto s : segments) {
                        std::pair<unsigned int, unsigned int> seg = s.first;
                        IsoformMultiSet<T> iSet = s.second;
                        std::string geneNamewithCopy = geneName;
                        if (segments.size() > 1) geneNamewithCopy += "#" + std::to_string(gene_copy);
                        gene_copy++;
                        currentGene = std::make_shared<Gene < T >> (geneNamewithCopy, seg.first, seg.second);
                        currentGene->setIsoforms(iSet);
                        SPtrIsoform<T> i = *(iSet.begin());
                        currentGene->setStrand(i->getStrand());
                        currentGene->createGeneFeatures(intronCutOff);
                        std::pair < GeneUnMapItr<T>, bool> res = genesNameIndex.insert(make_pair(geneNamewithCopy, currentGene));
                        if (!res.second) {
                            std::cerr << "Error inserting gene from segments" << std::endl;
                            exit(-1);
                        }
                        for (auto iIt : iSet) {
                            i = iIt;
                            transcript2Chr.insert(make_pair(i->getId(), std::make_pair(currentGene, i)));
                        }
                        genes.insert(currentGene);
                    }
                }
            }
            geneIsoforms.clear();
            this->createGeneUniqueFeatures();
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
            std::unordered_map<std::string, std::string> fieldsMap;
            GeneUnMapItr<T> it;
            int wStart, wEnd;

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
            if (geneName.find("-AS") == std::string::npos) {
                auto gItr = geneIsoforms.find(geneName);
                if (gItr == geneIsoforms.end()) {
                    std::shared_ptr<IsoformUnMap < T>> iMap = std::make_shared<IsoformUnMap < T >> ();
                    currentIsoform = std::make_shared<Isoform < T >> (isoformName, wStart, wEnd);
                    currentIsoform->setStrand(words[6][0]);
                    currentIsoform->processGTFLine(words);
                    currentIsoform->setFields(fieldsMap);
                    std::pair < IsoformUnMapItr<T>, bool> res = iMap->insert(std::make_pair(isoformName, currentIsoform));
                    if (!res.second) {
                        std::cerr << "Error inserting isoform" << std::endl;
                        exit(-1);
                    }
                    geneIsoforms.insert(std::make_pair(geneName, iMap));
                } else {
                    std::shared_ptr<IsoformUnMap < T>> iMap = gItr->second;
                    try {
                        if (!currentIsoform || currentIsoform->getId().compare(isoformName) != 0) {
                            IsoformUnMapItr<T> it;
                            it = iMap->find(isoformName);
                            if (it == iMap->end()) {
                                throw exceptions::NotFoundException("Isoform with name: " + isoformName + " does not exist");
                            }
                            currentIsoform = it->second;
                        }
                    } catch (exceptions::NotFoundException ex) {
                        currentIsoform = std::make_shared<Isoform < T >> (isoformName, wStart, wEnd);
                        std::pair < IsoformUnMapItr<T>, bool> res = iMap->insert(std::make_pair(isoformName, currentIsoform));
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
            }
        }

        friend std::ostream& operator<<(std::ostream& os, std::shared_ptr<Chromosome<T>> obj) {
            os << obj->getId() << " Genes: " << obj->getGenes().size() << std::endl;
            for (auto g : obj->getGenes()) {
                os << "\t" << g;
                os << std::endl;
            }
            return os;
        }

    private:
        std::string id;
        GeneMultiSet<T> genes;
        GeneUnMap<T> genesNameIndex;
        SPtrGene<T> currentGene;
        GeneIsoformReadUnMap<T> geneIsoforms;
        SPtrIsoform<T> currentIsoform;
    };

    template<typename T>
    using SPtrChromosome = std::shared_ptr<Chromosome<T>>;

    template<typename T>
    using ChromosomeUnMap = std::unordered_map<std::string, SPtrChromosome<T>>;

    template<typename T>
    using ChromosomeUnMapItr = typename ChromosomeUnMap<T>::iterator;

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
            SPtrFeature<T> f;

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
                c->createGeneMultiSetItr(intronCutOff, featuresToCreate, transcript2Chr);
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

