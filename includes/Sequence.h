/* 
 * File:   Sequence.h
 * Author: veraalva
 *
 * Created on May 5, 2016, 9:25 AM
 */

#ifndef SEQUENCE_H
#define SEQUENCE_H

namespace sequence {

    class Sequence {
    public:

        Sequence() {
            this->seq = "";
            this->id = "";
            this->description = "";
        }

        virtual ~Sequence() {
        }

        Sequence(const Sequence& other) :
        id(other.id), description(other.description), seq(other.seq) {
        }

        std::string getSegment(unsigned long int pos, unsigned long int length) {
            return seq.substr(pos, length);
        }

        std::string getId() {
            return id;
        }

        void setId(std::string id) {
            this->id = id;
        }

        std::string &getSeq() {
            return seq;
        }

        void setSeq(std::string seq) {
            this->seq = seq;
        }

        std::string getDescription() {
            return description;
        }

        void setDescription(std::string desc) {
            description = desc;
        }

        unsigned long int getLength() {
            return seq.size();
        }

        void shuffle() {
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::shuffle(seq.begin(), seq.end(), std::default_random_engine(seed));
        }

        void reverse() {
            std::reverse(seq.begin(), seq.end());
        }

        Sequence newSegment(unsigned long int pos, unsigned long int length) {
            Sequence s;
            s.setId(this->id);
            s.setDescription(this->description);
            s.setSeq(this->getSegment(pos, length));
            return s;
        }

    private:
        std::string id;
        std::string description;
        std::string seq;
    };

    class DNA : public Sequence {
    public:

        DNA() :
        Sequence() {
        }

        virtual ~DNA() {
        }
        
//        DNA shuffle() {
//            DNA s(*this);
//            s.shuffle();
//            return s;
//        }

        DNA complement() {
            DNA s(*this);
            for (auto it = s.getSeq().begin(); it < s.getSeq().end(); ++it) {
                switch (*it) {
                    case 'A': *it = 'T';
                        break;
                    case 'T': *it = 'A';
                        break;
                    case 'C': *it = 'G';
                        break;
                    case 'G': *it = 'C';
                        break;
                    case 'U': *it = 'A';
                        break;
                    case 'R': *it = 'Y';
                        break;
                    case 'Y': *it = 'R';
                        break;
                    case 'K': *it = 'M';
                        break;
                    case 'M': *it = 'K';
                        break;
                    case 'B': *it = 'V';
                        break;
                    case 'V': *it = 'B';
                        break;
                    case 'D': *it = 'H';
                        break;
                    case 'H': *it = 'D';
                        break;

                    case 'a': *it = 't';
                        break;
                    case 't': *it = 'a';
                        break;
                    case 'c': *it = 'g';
                        break;
                    case 'g': *it = 'c';
                        break;
                    case 'u': *it = 'a';
                        break;
                    case 'r': *it = 'y';
                        break;
                    case 'y': *it = 'r';
                        break;
                    case 'k': *it = 'm';
                        break;
                    case 'm': *it = 'k';
                        break;
                    case 'b': *it = 'v';
                        break;
                    case 'v': *it = 'b';
                        break;
                    case 'd': *it = 'h';
                        break;
                    case 'h': *it = 'd';
                        break;
                }
            }
            return s;
        }

        DNA reverseComplement() {
            DNA s = complement();
            s.reverse();
            return s;
        }

        DNA newSegment(unsigned long int pos, unsigned long int length) {
            DNA s;
            s.setId(this->getId());
            s.setDescription(this->getDescription());
            s.setSeq(this->getSegment(pos, length));
            return s;
        }
    };
    
    typedef std::shared_ptr<sequence::DNA> SPtrDNA;
    typedef std::unordered_map<std::string, SPtrDNA> TDNAMap;

    /**
     * Class to store and manipulate sequences
     */
    class DNAContainer {
    public:

        DNAContainer() {
        }

        virtual ~DNAContainer() {
        }

        void clearContainer() {
            sequences.clear();
        }

        TDNAMap &getContainer() {
            return sequences;
        }

        SPtrDNA getFirstElement() {
            TDNAMap::iterator it = sequences.begin();
            if (it == sequences.end()) {
                throw exceptions::NotFoundException("Not sequences on the container");
            }
            return it->second;
        }

        SPtrDNA getDNAFromID(std::string id) {
            TDNAMap::iterator it = sequences.find(id);
            if (it == sequences.end()) {
                throw exceptions::NotFoundException("Id " + id + " was not found in the sequence container");
            }
            return it->second;
        }

        std::pair<SPtrDNA, bool> addElement(std::string id) {
            std::pair < TDNAMap::iterator, bool> result;
            result = sequences.insert(std::make_pair(id, std::make_unique<DNA>()));
            return std::make_pair(result.first->second, result.second);
        }

        unsigned long int size() {
            return sequences.size();
        }
    private:
        TDNAMap sequences;
    };

}

#endif /* SEQUENCE_H */

