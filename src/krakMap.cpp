#include <fstream>
#include <iostream>
#include <string>
#include <algorithm> // std::sort
#include "CLI/CLI.hpp"
#include "krakMap.h"

#define LEFT true
#define RIGHT true

struct krakMapOpts {
    std::string taxonomyTree_filename;
    std::string refId2TaxId_filename;
    std::string mapperOutput_filename;
    std::string output_filename;
    std::string level = "species";
    double filterThreshold = 0;
};

void TaxaNode::addInterval(uint64_t begin, uint64_t len, bool isLeft) {
    if (isLeft)
        lintervals.emplace_back(begin, begin+len);
    else
        rintervals.emplace_back(begin, begin+len);
}

void TaxaNode::updateScore() {
    score = 0;
    for (auto& it : lintervals) {
        score += it.end - it.begin;
    }
    for (auto& it : rintervals) {
        score += it.end - it.begin;
    }
    //std::cout << " score: " << score << "\n";
}
void TaxaNode::updateIntervals(TaxaNode* child, bool isLeft) {

    std::vector<Interval>* intervals;
    
    // merge two sorted interval lists into parent
    // update parent score
    auto& childIntervals = child->getIntervals(isLeft);
    if (isLeft)
        intervals = &lintervals;
    else
        intervals = &rintervals;
    std::vector<Interval> parentIntervals(intervals->size());
    std::copy(intervals->begin(), intervals->end(), parentIntervals.begin());
    intervals->clear();
    intervals->reserve(parentIntervals.size()+childIntervals.size());
    
    std::vector<Interval>::iterator pit = parentIntervals.begin();
    std::vector<Interval>::iterator cit = childIntervals.begin();
    std::vector<Interval>::iterator fit = intervals->begin();

    // add the smallest interval as the first interval
    if (cit != childIntervals.end() && (pit == parentIntervals.end() || cit->begin < pit->begin)) {
        intervals->emplace_back(cit->begin, cit->end);
        cit++;
    }
    else if (pit != parentIntervals.end()) {
        intervals->emplace_back(pit->begin, pit->end);
        pit++;
    }
    else {
        std::cerr << "ERROR!! Both parent an child intervals were empty.\n";
        std::exit(1);
    }
    std::vector<Interval>::iterator cur;
    while (pit != parentIntervals.end() || cit != childIntervals.end()) {
        // find the smallest interval between the heads of the two lists
        if (cit != childIntervals.end() && (pit == parentIntervals.end() || cit->begin < pit->begin) )  {
            cur = cit;
            cit++;
        }
        else if (pit != parentIntervals.end()) {
            cur = pit;
            pit++;  
        }
        else {
            std::cerr << "ERROR!! Shouldn't even enter the loop.\n";
            std::exit(1);
        }
        // merge the new interval
        // Note: since both lists are sorted
        // the new interval's begin is always >= the last inserted interval's
        if (fit->end >= cur->begin) { // if the new interval has an overlap with the last inserted one
            fit->end = std::max(cur->end, fit->end); // merge them
        } else { // insert the interval as a separate one and move fit forward
            intervals->emplace_back(cur->begin, cur->end);
            fit++;
        }
    }
}
/**
 * Sorts intervals
 * Merge intervals if possible
 * Calculates score
**/
void TaxaNode::cleanIntervals(bool isLeft) {
    // if we were writing intervals wrt read position we wouldn't need this part
    std::vector<Interval>* intervals;
    if (isLeft)
        intervals = &lintervals;
    else
        intervals = &rintervals;
    
    std::sort(intervals->begin(), intervals->end(), 
    [](Interval& i1, Interval& i2){
        return i1.begin != i2.begin?i1.begin < i2.begin:i1.end < i2.end;
    });
    // merge intervals if necessary
    // calculate score / coverage !! this whole process is repetition of mapping coverage calc!!
    for (auto it=intervals->begin(); it != intervals->end();) {
        // start from next item and merge (and erase) as much as possible
        bool merged = true;
        auto next = it+1;
        while (next != intervals->end() && merged) {
            // if they overlap, merge them
            if (it->end >= next->begin) {
                if (it->end < next->end)
                    it->end = next->end;
                intervals->erase(next); // erase next after merging with it
                next = it+1; // it pointer and its next are always valid
            }
            else {// no overlap anymore and no merging. Update score and leave the loop
                merged = false;
            }
        }
        // there is nothing to merge with the current it, so increase it
        it++;
    }
}

bool TaxaNode::addChild(TaxaNode* child) { 
    // if this child was not already added
    if (activeChildren.find(child->getId()) == activeChildren.end()) {
        activeChildren.insert(child->getId()); // add it
        notIncorporatedChildrenCounter++; // increase not incorprated children counter
        return true;
    }
    return false;
}

void TaxaNode::reset() {
     lintervals.clear();
     rintervals.clear();
     activeChildren.clear();
     notIncorporatedChildrenCounter = 0;
     score = 0;
}


KrakMap::KrakMap(std::string& taxonomyTree_filename, 
                 std::string& refId2TaxId_filename, 
                 std::string pruneLevelIn,
                 double filteringThresholdIn) {

    std::cerr << "KrakMap: Construct ..\n";
    // map rank string values to enum values
    filteringThreshold = filteringThresholdIn;
    initializeRanks();
    pruningLevel = str2Rank[pruneLevelIn];
    std::ifstream tfile;
    uint32_t id, pid;
    std::string rank, name;
    // load the reference id (name) to its corresponding taxonomy id
    tfile.open(refId2TaxId_filename);
    while(!tfile.eof()) {
        tfile >> name >> id;
        /* if (refId2taxId.find(name) != refId2taxId.end()) {
            std::cerr << "MULTIPLE REF_NAMES: " << name << " " << id << " --> ";
            std::cerr << refId2taxId[name] << " " << refId2taxId.size() << "\n";
        } */
        refId2taxId[name] = id;
    }
    tfile.close();

    // load the taxonomy child-parent tree and the rank of each node
    tfile.open(taxonomyTree_filename);
    std::string tmp;
    while (!tfile.eof()) {
        tfile >> id >> tmp >> pid >> tmp >> rank >> tmp;
        if (tmp != "|") {
            rank += " " + tmp;
        }
        /* if (taxaNodeMap.find(id) != taxaNodeMap.end()) {
            std::cerr << "MULTIPLE NODE IDs: " << id << " " << pid << " " << rank << " --> ";
            std::cerr << taxaNodeMap[id].getId() << " " << taxaNodeMap[id].getParentId() << " " << rankToStr(taxaNodeMap[id].getRank()) << " " << taxaNodeMap.size() << "\n";
            std::exit(1);
        } */
        taxaNodeMap[id] = TaxaNode(id, pid, str2Rank[rank]);
        if (taxaNodeMap[id].isRoot()) {
            rootId = id;
            std::cerr << "Root Id : " << id << "\n";
        }
        std::getline(tfile, tmp);
        
    }

    tfile.close();  
}

bool KrakMap::readHeader(std::ifstream& mfile) {
    std::string tmp, readType;
    mfile >> tmp >> readType;
    if (tmp != "#")
        return false;
    if (readType == "LT:S") 
        isPaired = false;
    else if (readType == "LT:P")
        isPaired = true;
    else
        return false;
    std::getline(mfile, tmp);
    return true;
}

void KrakMap::loadMappingInfo(std::ifstream& mfile) {
    std::string tname, tmp;
    uint64_t lcnt, rcnt, tid, ibeg, ilen;
    mfile >> tname;
    // first condition: Ignore those references that we don't have a taxaId for
    // secon condition: Ignore repeated exactly identical mappings (FIXME thing)
    if (refId2taxId.find(tname) != refId2taxId.end() &&
        activeTaxa.find(refId2taxId[tname]) == activeTaxa.end()) { 
        tid = refId2taxId[tname];
        activeTaxa.insert(tid);
        
        // fetch the taxon from the map
        TaxaNode* taxaPtr = &taxaNodeMap[tid];
        walk2theRoot(taxaPtr);
        mfile >> lcnt;
        if (isPaired)
            mfile >> rcnt;
        for (size_t i = 0; i < lcnt; ++i) {
            mfile >> ibeg >> ilen;
            taxaPtr->addInterval(ibeg, ilen, LEFT);
        }
        if (isPaired)
            for (size_t i = 0; i < rcnt; ++i) {
                mfile >> ibeg >> ilen;
                taxaPtr->addInterval(ibeg, ilen, RIGHT);
            }
        taxaPtr->cleanIntervals(LEFT);
        taxaPtr->cleanIntervals(RIGHT);
        taxaPtr->updateScore();
        hits.push_front(taxaPtr);
    }
    else { // otherwise we have to read till the end of the line and throw it away
        std::getline(mfile, tmp);
    }
}

bool KrakMap::classify(std::string& mapperOutput_filename) {
    std::cerr << "KrakMap: Classify ..\n";
    std::cerr << "\tMapping Output File: " << mapperOutput_filename << "\n";
    std::ifstream mfile(mapperOutput_filename);
    std::string rid, tname, tmp; // read id, taxa name temp
    uint64_t rlen, mcnt; // taxa id, read mapping count, # of interals, interval start, interval length
    uint64_t totalReadCnt = 0, totalUnmappedReads = 0, seqNotFound = 0;
    if (!readHeader(mfile)) {
        std::cerr << "ERROR: Invalid header for mapping output file.\n";
        std::exit(1);
    }
    std::cout<< "paired end? " << isPaired << "\n";
    while (!mfile.eof()) {
        mfile >> rid >> mcnt;
        totalReadCnt++;
        //std::cout << "r" << rid << " " << mcnt << "\n";
        if (mcnt != 0) {
            if (isPaired) {
                uint64_t rllen, rrlen;
                mfile >> rllen >> rrlen;
                rlen = rllen + rrlen;
            }
            else {
                mfile >> rlen;
            }
            std::set<uint64_t> seen;
            // reset everything we've done for previous read
            clearReadSubTree();
            // std::cerr << activeTaxa.size() << " ";
            // construct intervals for leaves
            for (size_t mappingCntr = 0; mappingCntr < mcnt; mappingCntr++) {
                loadMappingInfo(mfile);   
            } 
            if (activeTaxa.size() == 0) {
                seqNotFound++;
            }
            else {
                // propagate score and intervals to all internal nodes
                // std::cerr << "Update intervals and scores of internal nodes ..\n";
                propagateInfo();

                // find best path for this read
                // std::cerr << "Assign Read ..\n";
                assignRead(rlen);
            }
        } else {
                totalUnmappedReads++;
                std::getline(mfile, tmp);
        }
    }
    std::cout << "\nTotal Reads: " << totalReadCnt << "\n"
              << "Total Mapped Reads: " << totalReadCnt - totalUnmappedReads << "\n"
              << "Sequence Not Found: " << seqNotFound << "\n"
              << "Reads Classified: " << readCntr << "\n";
    return true;
}

void KrakMap::walk2theRoot(TaxaNode* child) {
    while (!child->isRoot()) {
        TaxaNode* parent = &taxaNodeMap[child->getParentId()]; // fetch parent
        activeTaxa.insert(parent->getId());
        parent->addChild(child); // add current node as its child
        child = parent; // walk to parent --> parent is the new child
    }
}

//TODO don't like it 
// It's not as beautiful as when the root is ripe, we are done, but it does the same thing :/
// will we pass being ripe and become overripe? (negative # of incorporated children .. )
// we shouldn't .. but it just doesn't happen because it doesn't .. :/
void KrakMap::propagateInfo() {
    while (!hits.empty()) {
        TaxaNode* taxaPtr = hits.back();
        hits.pop_back();
        // if the hit itself is not ripe, no need to push it back to the queue
        // when it's the turn for one of the other hits, it'll get ripe and updated
        while (!taxaPtr->isRoot() && taxaPtr->isRipe()) {
            TaxaNode* parentPtr = &taxaNodeMap[taxaPtr->getParentId()];
            parentPtr->updateIntervals(taxaPtr, LEFT);
            parentPtr->updateIntervals(taxaPtr, RIGHT);
            parentPtr->updateScore();
            parentPtr->setOneMoreChildAsProcessed();
            taxaPtr = parentPtr;
        }
    }
}

// why? why? why I don't understand this pointer/reference thing :mad: :dissapointed:
void KrakMap::assignRead(uint64_t readLen) {
	int32_t k{31};
    //double threshold = filteringThreshold >= 1 ? filteringThreshold : readLen*filteringThreshold;
    double threshold = filteringThreshold >= 1 ? filteringThreshold : ((readLen - k + 1)*filteringThreshold + (k-1));
    TaxaNode* walker = &taxaNodeMap[rootId];
    //std::cerr << walker->getScore() << " " << threshold << "\n";
    if (walker->getScore()<threshold) {
         if (mappedReadCntr.find(0) == mappedReadCntr.end())
            mappedReadCntr.insert(std::make_pair(0, TaxaInfo(1, Rank::UNCLASSIFIED)));
            //mappedReadCntr[0] = TaxaInfo(1, walker->getRank());
        else
            mappedReadCntr[0].increment();
        return;
    }
    while (walker->getRank() != pruningLevel) {
        uint64_t maxScore=0, maxId = -1, maxCntr = 0;
        for (auto childId : walker->getActiveChildren()) {
            
            TaxaNode& child = taxaNodeMap[childId];
            if (child.getScore() == maxScore) {
                maxCntr++;
            }
            else if (child.getScore() > maxScore) {
                maxId = childId;
                maxScore = child.getScore();
                maxCntr = 1;
            }
        }
        // zero --> no children (it's a leaf) || > 1 --> more than one child with max score
        if (maxCntr != 1 || taxaNodeMap[maxId].getScore()<threshold) { 
            break;
        }

        walker = &taxaNodeMap[maxId];
    }
    readCntr++;
    if (mappedReadCntr.find(walker->getId()) == mappedReadCntr.end()) {
        mappedReadCntr.insert(std::make_pair(walker->getId(), TaxaInfo(1, walker->getRank())));
    }
        //mappedReadCntr[walker->getId()] = TaxaInfo(1, walker->getRank());
    else
        mappedReadCntr[walker->getId()].increment();
    
    while (walker->getId() != rootId) {
        walker = &taxaNodeMap[walker->getParentId()];
        if (mappedReadCntr.find(walker->getId()) == mappedReadCntr.end())
            mappedReadCntr.insert(std::make_pair(walker->getId(), TaxaInfo(0, walker->getRank())));
        mappedReadCntr[walker->getId()].subTreeCnt++;     
    }
}

void KrakMap::serialize(std::string& output_filename) {
    std::cerr << "Write results in the file\n";
    std::ofstream ofile(output_filename);
    ofile << "taxaId\ttaxaRank\tcount\tsubTreeCount\n";
    for (auto& kv : mappedReadCntr) {
        ofile << kv.first << "\t" << rankToStr(kv.second.rank) << "\t" << kv.second.cnt << "\t" << kv.second.subTreeCnt << "\n";
    }
    ofile.close();
}
    
void KrakMap::clearReadSubTree() {
    for (auto & activeTaxum : activeTaxa) {
        taxaNodeMap[activeTaxum].reset();
    }
    activeTaxa.clear();
    hits.clear();
}


/**
 * "How to run" example:
 * make Pufferfish!
 * In the Pufferfish build directory run the following command:
 * /usr/bin/time src/krakmap 
 * -t /mnt/scratch2/avi/meta-map/kraken/KrakenDB/taxonomy/nodes.dmp  
 * -s /mnt/scratch2/avi/meta-map/kraken/KrakenDB/seqid2taxid.map 
 * -m /mnt/scratch2/avi/meta-map/kraken/puff/dmps/HC1.dmp 
 * -o HC1.out 
 * -l genus (optional)
 * -f 0.8 (optional)
 **/
int main(int argc, char* argv[]) {
  (void)argc;

  krakMapOpts kopts;
  CLI::App app{"krakMap : Taxonomy identification based on the output of Pufferfish mapper through the same process as Kraken."};
  app.add_option("-t,--taxtree", kopts.taxonomyTree_filename,
                 "path to the taxonomy tree file")
      ->required();
  app.add_option("-s,--seq2taxa", kopts.refId2TaxId_filename, "path to the refId 2 taxaId file")
      ->required();
  app.add_option("-m,--mapperout", kopts.mapperOutput_filename, "path to the pufferfish mapper output file")
      ->required();
  app.add_option("-o,--output", kopts.output_filename, "path to the output file to write results")
      ->required();
  app.add_option("-l,--level", kopts.level, "choose between (species, genus, family, order, class, phylum). Default:species")
      ->required(false);
  app.add_option("-f,--filter", kopts.filterThreshold, "choose the threshold (0-1) to filter out mappings with a score below that. Default: no filter")
      ->required(false);

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    return app.exit(e);
  }
  KrakMap krakMap(kopts.taxonomyTree_filename, kopts.refId2TaxId_filename, kopts.level, kopts.filterThreshold);
  krakMap.classify(kopts.mapperOutput_filename);
  krakMap.serialize(kopts.output_filename);
  return 0;
}
