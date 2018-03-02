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
            this->TPM = 0.0;
            this->uniqueTPM = 0.0;
            this->uniqueTPMIntron = 0.0;
            this->uniqueReads = 0;
            this->uniqueLength = 0;
            this->uniqueReadsIntron = 0;
            this->uniqueIntronLength = 0;
            this->TPMExon = 0.0;
            this->TPMIntron = 0.0;
            this->exonReads = 0;
            this->exonLength = 0;
            this->intronReads = 0;
            this->intronLength = 0;
            this->reads = 0;            
            this->TPMBridges = 0.0;
            this->bridgeReads = 0;
        }

        SampleData(int reads) {
            this->TPM = 0.0;
            this->uniqueTPM = 0.0;
            this->uniqueTPMIntron = 0.0;
            this->uniqueReads = 0;
            this->uniqueLength = 0;
            this->uniqueReadsIntron = 0;
            this->uniqueIntronLength = 0;
            this->TPMExon = 0.0;
            this->TPMIntron = 0.0;
            this->exonReads = 0;
            this->exonLength = 0;
            this->intronReads = 0;
            this->intronLength = 0;
            this->reads = reads;            
            this->TPMBridges = 0.0;
            this->bridgeReads = 0;
        }

        virtual ~SampleData() {
        }

        double getTPM() const {
            return TPM;
        }

        void setTPM(double TPM) {
            this->TPM = TPM;
        }

        double getUniqueTPM() const {
            return uniqueTPM;
        }

        void setUniqueTPM(double uniqueTPM) {
            this->uniqueTPM = uniqueTPM;
        }

        double getUniqueTPMIntron() const {
            return uniqueTPMIntron;
        }

        void setUniqueTPMIntron(double uniqueTPMIntron) {
            this->uniqueTPMIntron = uniqueTPMIntron;
        }

        double getTPMExon() const {
            return TPMExon;
        }

        void setTPMExon(double TPMExon) {
            this->TPMExon = TPMExon;
        }

        double getTPMIntron() const {
            return TPMIntron;
        }

        void setTPMIntron(double TPMIntron) {
            this->TPMIntron = TPMIntron;
        }

        int getExonReads() const {
            return exonReads;
        }

        void increaseExonReads(int exonCount) {
            this->exonReads += exonCount;
        }

        int getExonLength() const {
            return exonLength;
        }

        void increaseExonLength(int exonLength) {
            this->exonLength += exonLength;
        }

        int getIntronReads() const {
            return intronReads;
        }

        void increaseIntronReads(int intronCount) {
            this->intronReads += intronCount;
        }

        int getIntronLength() const {
            return intronLength;
        }

        void increaseIntronLength(int intronLength) {
            this->intronLength += intronLength;
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

        int getUniqueReads() const {
            return uniqueReads;
        }

        void increaseUniqueReads() {
            this->uniqueReads++;
        }

        void increaseUniqueReads(int uniqueReads) {
            this->uniqueReads += uniqueReads;
        }

        int getUniqueReadsIntron() const {
            return uniqueReadsIntron;
        }

        void increaseUniqueReadsIntron() {
            this->uniqueReadsIntron++;
        }

        void increaseUniqueReadsIntron(int uniqueReadsIntron) {
            this->uniqueReadsIntron += uniqueReadsIntron;
        }

        int getUniqueIntronLength() const {
            return uniqueIntronLength;
        }

        void setUniqueIntronLength(int uniqueIntronLength) {
            this->uniqueIntronLength = uniqueIntronLength;
        }

        void increaseUniqueIntronLength(int uniqueIntronLength) {
            this->uniqueIntronLength += uniqueIntronLength;
        }

        int getUniqueLength() const {
            return uniqueLength;
        }

        void setUniqueLength(int uniqueLength) {
            this->uniqueLength = uniqueLength;
        }

        void increaseUniqueLength(int uniqueLength) {
            this->uniqueLength += uniqueLength;
        }

        float getTPMBridges() const {
            return TPMBridges;
        }

        void setTPMBridges(float TPMBridges) {
            this->TPMBridges = TPMBridges;
        }

        int getBridgeReads() const {
            return bridgeReads;
        }

        void increaseBridgeReads() {
            this->bridgeReads++;
        }

    private:
        int reads;
        int uniqueReads;
        int uniqueLength;
        int uniqueReadsIntron;
        int uniqueIntronLength;
        double TPM;
        double uniqueTPM;
        double uniqueTPMIntron;
        int exonLength;
        int exonReads;
        int intronReads;
        int intronLength;
        double TPMExon;
        double TPMIntron;
        float TPMBridges;
        int bridgeReads;
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

        int processBAMSAMFromDir(std::string dirName, bool onlyProperlyPaired);
        int processReadsFromBAM(std::string bamFileName, std::string sampleName, bool onlyProperlyPaired);
        int processReadsFromIsoformBAM(std::string bamFileName, std::string sample);
        std::vector<BamTools::CigarOp> processCigar(std::string cigar);
        void printResults(bool singleFile);

        void processReadAtGenomeLevel(std::string chrName, std::string sampleName, unsigned int start, unsigned int end);
        void processReadAtGenomeLevelUnique(std::string chrName, std::string sampleName, unsigned int start, unsigned int end);
        void processReadAtGeneLevel(SPtrGeneNGS gene, std::string sampleName, unsigned int start, unsigned int end);
        void processReadAtGeneLevelUnique(SPtrGeneNGS gene, std::string sampleName, unsigned int start, unsigned int end);
        void processReadAtIsoformLevel(SPtrIsoformNGS isoform, std::string sampleName, unsigned int start, unsigned int end);

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
        void PopulateReads(std::string sampleName);
    };
}


#endif /* READFACTORY_H */

