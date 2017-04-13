/* 
 * File:   FastaFactory.h
 * Author: veraalva
 *
 * Created on February 10, 2016, 3:41 PM
 */

#ifndef FASTAFACTORY_H
#define FASTAFACTORY_H

namespace formats {

    class FastaFactory {
    public:

        FastaFactory() {
        }

        virtual ~FastaFactory() {
        }

        static void parseDNAFastaInDirectory(sequence::DNAContainer &seqContainer, std::string dirName, std::string prefix, std::string sufix, bool binary);
        static long unsigned int parseDNAFastaFile(sequence::DNAContainer &seqContainer, std::string fName, bool binary);
        static void writeDNASequencesToFile(sequence::DNAContainer &seqContainer, std::string fileName, bool binary);
    private:

    };

}

#endif /* FASTAFACTORY_H */

