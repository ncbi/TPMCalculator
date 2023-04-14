/* 
 * File:   ReadFactory.h
 * Author: veraalva
 *
 * Created on May 6, 2016, 10:49 AM
 */

#ifndef READFACTORY_H
#define READFACTORY_H

#include "GenomeFactory.h"


namespace ngs {

    class SampleData {
    public:

        SampleData() {
            this->reads = 0;
            this->TPM = 0.0;

            this->exonReads = 0;
            this->exonTPM = 0.0;

            this->intronReads = 0;
            this->intronTPM = 0.0;

            this->uniqueReads = 0;
            this->uniqueTPM = 0.0;

            this->uniqueExonReads = 0;
            this->uniqueExonTPM = 0.0;

            this->uniqueIntronReads = 0;
            this->uniqueIntronTPM = 0.0;

            this->bridgesReads = 0;
            this->bridgesTPM = 0.0;
        }

        SampleData(int reads) {
            this->reads = reads;
            this->TPM = 0.0;

            this->exonReads = 0;
            this->exonTPM = 0.0;

            this->intronReads = 0;
            this->intronTPM = 0.0;

            this->uniqueReads = 0;
            this->uniqueTPM = 0.0;

            this->uniqueExonReads = 0;
            this->uniqueExonTPM = 0.0;

            this->uniqueIntronReads = 0;
            this->uniqueIntronTPM = 0.0;

            this->bridgesReads = 0;
            this->bridgesTPM = 0.0;
        }

        virtual ~SampleData() {
        }

        int getReads() const {
            return reads;
        }

        void increaseReads() {
            this->reads++;
        }

        void increaseReads(int reads) {
            this->reads += reads;
        }

        double getTPM() const {
            return TPM;
        }

        void setTPM(double TPM) {
            this->TPM = TPM;
        }

        int getExonReads() const {
            return exonReads;
        }

        void increaseExonReads() {
            this->exonReads++;
        }

        void increaseExonReads(int exonCount) {
            this->exonReads += exonCount;
        }

        double getExonTPM() const {
            return exonTPM;
        }

        void setExonTPM(double exonTPM) {
            this->exonTPM = exonTPM;
        }

        int getIntronReads() const {
            return intronReads;
        }

        void increaseIntronReads(int intronCount) {
            this->intronReads += intronCount;
        }

        void increaseIntronReads() {
            this->intronReads++;
        }

        double getIntronTPM() const {
            return intronTPM;
        }

        void setIntronTPM(double intronTPM) {
            this->intronTPM = intronTPM;
        }

        double getUniqueTPM() const {
            return uniqueTPM;
        }

        void setUniqueTPM(double uniqueTPM) {
            this->uniqueTPM = uniqueTPM;
        }

        int getUniqueReads() const {
            return uniqueReads;
        }

        void increaseUniqueReads() {
            this->uniqueReads++;
        }

        void increaseUniqueReads(int uniqueReads) {
            this->uniqueReads += uniqueReads;
        }

        int getUniqueExonReads() {
            return uniqueExonReads;
        }

        void increaseUniqueExonReads() {
            this->uniqueExonReads++;
        }

        void increaseUniqueExonReads(int uniqueReadsExon) {
            this->uniqueExonReads += uniqueReadsExon;
        }

        double getUniqueExonTPM() const {
            return uniqueExonTPM;
        }

        void setUniqueExonTPM(double uniqueExonTPM) {
            this->uniqueExonTPM = uniqueExonTPM;
        }

        int getUniqueIntronReads() const {
            return uniqueIntronReads;
        }

        void increaseUniqueIntronReads() {
            this->uniqueIntronReads++;
        }

        void increaseUniqueIntronReads(int uniqueReadsIntron) {
            this->uniqueIntronReads += uniqueReadsIntron;
        }

        double getUniqueIntronTPM() const {
            return uniqueIntronTPM;
        }

        void setUniqueIntronTPM(double uniqueIntronTPM) {
            this->uniqueIntronTPM = uniqueIntronTPM;
        }

        int getBridgesReads() const {
            return bridgesReads;
        }

        void increaseBridgesReads() {
            this->bridgesReads++;
        }

        double getBridgesTPM() const {
            return bridgesTPM;
        }

        void setBridgesTPM(double bridgesTPM) {
            this->bridgesTPM = bridgesTPM;
        }

        double getValueFromColumn(std::string column) {
            if (column.compare("Reads") == 0) {
                return this->getReads();
            } else if (column.compare("TPM") == 0) {
                return this->getTPM();
            } else if (column.compare("ExonReads") == 0) {
                return this->getExonReads();
            } else if (column.compare("ExonTPM") == 0) {
                return this->getExonTPM();
            } else if (column.compare("IntronReads") == 0) {
                return this->getIntronReads();
            } else if (column.compare("IntronTPM") == 0) {
                return this->getIntronTPM();
            } else if (column.compare("UniqueReads") == 0) {
                return this->getUniqueReads();
            } else if (column.compare("UniqueTPM") == 0) {
                return this->getUniqueTPM();
            } else if (column.compare("UniqueExonReads") == 0) {
                return this->getUniqueExonReads();
            } else if (column.compare("UniqueExonTPM") == 0) {
                return this->getUniqueExonTPM();
            } else if (column.compare("UniqueIntronReads") == 0) {
                return this->getUniqueIntronReads();
            } else if (column.compare("UniqueIntronTPM") == 0) {
                return this->getUniqueIntronTPM();
            } else if (column.compare("BridgesReads") == 0) {
                return this->getBridgesReads();
            } else if (column.compare("BridgesTPM") == 0) {
                return this->getBridgesTPM();
            }
            return 0.0;
        }

    private:
        int reads;
        double TPM;

        int exonReads;
        double exonTPM;

        int intronReads;
        double intronTPM;

        int uniqueReads;
        double uniqueTPM;

        int uniqueExonReads;
        double uniqueExonTPM;

        int uniqueIntronReads;
        double uniqueIntronTPM;

        int bridgesReads;
        float bridgesTPM;

    };

    typedef std::shared_ptr<SampleData> SPtrSampleData;
    typedef std::unordered_map<std::string, SPtrSampleData> SampleDataUnMap;

    class ReadData {
    public:

        ReadData() {
        }

        virtual ~ReadData() {
        }

        SampleDataUnMap getData() const {
            return data;
        }

        SPtrSampleData createSampleData(std::string sampleName) {
            SPtrSampleData sampleData;
            try {
                sampleData = getSampleData(sampleName);
            } catch (exceptions::NotFoundException) {
                sampleData = std::make_shared<SampleData>(SampleData(0));
                data.insert(std::pair<std::string, SPtrSampleData> (sampleName, sampleData));
            }
            return sampleData;
        }

        SPtrSampleData getSampleData(std::string sampleName) {
            SampleDataUnMap::iterator it = data.find(sampleName);
            if (it == data.end()) {
                throw exceptions::NotFoundException("Sample with name: " + sampleName + " does not exist");
            }
            return it->second;
        }

        void increaseReads(std::string sampleName) {
            SPtrSampleData sampleData = createSampleData(sampleName);
            sampleData->increaseReads();
        }

    private:
        SampleDataUnMap data;
    };

    typedef genome::GenomeFactory<ReadData> GenomeFactoryNGS;
    typedef genome::SPtrChromosome<ReadData> SPtrChromosomeNGS;
    typedef genome::SPtrGene<ReadData> SPtrGeneNGS;
    typedef genome::SPtrIsoform<ReadData> SPtrIsoformNGS;
    typedef genome::SPtrFeature<ReadData> SPtrFeatureNGS;
    typedef genome::GeneMultiSet<ReadData> GeneMultiSetNGS;

    class ReadFactory {
    public:
        ReadFactory();
        virtual ~ReadFactory();

        GenomeFactoryNGS& getGenomeFactory() {
            return genomeFactory;
        }

        std::vector<std::string>& getSamples() {
            return samples;
        }

        struct coordinateLessCMP {

            bool operator()(const std::pair<unsigned int, unsigned int> a, const std::pair<unsigned int, unsigned int> b) {
                return a.first < b.first;
            }
        };

        int processBAMSAMFromDir(std::string dirName, bool onlyProperlyPaired, uint16_t minMAPQ, uint16_t minOverlap);
        int processReadsFromBAM(std::string bamFileName, std::string sampleName, bool onlyProperlyPaired, uint16_t minMAPQ, uint16_t minOverlap);
        std::vector<BamTools::CigarOp> processCigar(std::string cigar);
        void printResults(bool singleFile, bool extendedOutput, bool all_feat, std::string output_path);
        void printResultsMatrix(std::string output_name, std::vector<std::string> tpmColumns);

        void processReadAtGenomeLevel(std::string chrName, std::string sampleName, std::set < std::pair<unsigned int, unsigned int>, coordinateLessCMP> read_coords, uint16_t minOverlap);
        //        void processReadAtGenomeLevelUnique(std::string chrName, std::string sampleName, unsigned int start, unsigned int end, uint16_t minOverlap);
        void processReadAtGeneLevel(SPtrGeneNGS gene, std::string sampleName, std::set < std::pair<unsigned int, unsigned int>, coordinateLessCMP> read_coords, uint16_t minOverlap);
        //        void processReadAtGeneLevelUnique(SPtrGeneNGS gene, std::string sampleName, unsigned int start, unsigned int end, uint16_t minOverlap);
        //        void processReadAtIsoformLevel(SPtrIsoformNGS isoform, std::string sampleName, unsigned int start, unsigned int end, uint16_t minOverlap);

        void loadTPMCalculatorGenesOutput(std::string dirName);

        void createSIMSingleReadsIR(std::string outFileName,
                sequence::DNAContainer seqContainer,
                unsigned int numberFeat, unsigned int intronNumber, unsigned int len);
    private:
        BamTools::BamReader reader;
        BamTools::SamHeader header;
        BamTools::RefVector references;
        GenomeFactoryNGS genomeFactory;
        std::vector<std::string> samples;

        void calculateTPMperSample(std::string sampleName);
        //void processReadAtGenomeLevel(std::string chrName, std::string sampleName, unsigned int start, unsigned int end);
        //        void processReadAtGeneLevel(std::shared_ptr<genome::Gene<ReadData>> gene, std::string sampleName, unsigned int start, unsigned int end);
        //        void processReadAtIsoformLevel(std::shared_ptr<genome::Isoform<ReadData>> isoform, std::string sampleName, unsigned int start, unsigned int end);
        //        void PopulateReads(std::string sampleName);
    };
}


#endif /* READFACTORY_H */

