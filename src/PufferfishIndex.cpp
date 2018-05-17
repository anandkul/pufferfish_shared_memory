#include <fstream>
#include <iostream>

#include "CLI/Timer.hpp"
#include "CanonicalKmerIterator.hpp"
#include "PufferFS.hpp"
#include "PufferfishIndex.hpp"
#include "cereal/archives/binary.hpp"
#include "cereal/archives/json.hpp"
#include <sys/ipc.h>
#include <sys/shm.h>
#include <errno.h>

#include "jellyfish/mer_dna.hpp"

static inline int
hash_string(const std::string& s) {
	int ret = 0;
	int a = 63689;
	int b = 378551;
	for(size_t i = 0; i < s.length(); i++) {
		ret = (ret * a) + (int)s[i];
		if(a == 0) {
			a += b;
		} else {
			a *= b;
		}
		if(a == 0) {
			a += b;
		}
	}
	return ret;
}

PufferfishIndex::PufferfishIndex() {}




PufferfishIndex::PufferfishIndex(const std::string& indexDir) {
  if (!puffer::fs::DirExists(indexDir.c_str())) {
    std::cerr << "The index directory " << indexDir << " does not exist!\n";
    std::exit(1);
  }
  

  {
    std::ifstream infoStream(indexDir + "/info.json");
    cereal::JSONInputArchive infoArchive(infoStream);
    infoArchive(cereal::make_nvp("k", k_));
    infoArchive(cereal::make_nvp("num_kmers", numKmers_));
    std::cerr << "k = " << k_ << '\n';
    std::cerr << "num kmers = " << numKmers_ << '\n';
    infoStream.close();
    twok_ = 2 * k_;
	/*
    	std::ifstream file(indexDir + "/reflengths.bin", std::ios::binary | std::ios::ate);
    	key_t key = (key_t)hash_string(indexDir + "/1111reflengths.bin" + "refLengths_");
	std::ifstream::pos_type size = file.tellg();
	size = size * 100;
   	int shmid = shmget(key, size, IPC_CREAT | IPC_EXCL);
	if(shmid != -1)
	{
		std::cerr << "1 ";
		shmctl(shmid, IPC_RMID, NULL);
		std::cerr << "2 ";
		shmid = shmget(key, size, IPC_CREAT | 0666);
		std::cerr << "3 ";
		int* ptr = (int*)shmat(shmid,NULL,0);
		std::cerr << "4 ";
		for(int i = 0; i < 50; i++)
		{
			std::cerr << "5 ";
			*ptr++ = 50-i+10;
		}
	}
	else
	{
		std::cerr << "6 ";
		shmid = shmget(key, size, IPC_CREAT | 0666);
		std::cerr << "7 ";
		int* ptr = (int*)shmat(shmid,NULL,0);
		std::cerr << "8 ";
		for(int i = 0; i < 50; i++)
		{
			std::cerr << "9 ";
			std::cerr << *ptr++ << " "; 
		}
	}
	*/

  } 

  {
	std::ifstream file(indexDir +  "/ctable.bin", std::ios::binary | std::ios::ate);
	key_t key = (key_t)hash_string(indexDir + "/qqqctable.bin" + "refLengths_");
    	std::ifstream::pos_type size = file.tellg();
    	size = size * 200;
    	//int shmid = shmget(key, size, IPC_CREAT | IPC_EXCL);
    	key_t key1 = (key_t)hash_string(indexDir + "/qqqctable.bin" + "refLengths_" + "vectorSize");
    	int shmid1 = shmget(key1, size/100, IPC_CREAT | IPC_EXCL);

	if(shmid1 != -1)
	{
	

		//shmctl(shmid, IPC_RMID, NULL);
		//shmid = shmget(key, size, IPC_CREAT | 0666);
		shmctl(shmid1, IPC_RMID, NULL);
		shmid1 = shmget(key1, size/100, IPC_CREAT | 0666);
		CLI::AutoTimer timer{"Loading contig table", CLI::Timer::Big};
    		std::ifstream contigTableStream(indexDir + "/ctable.bin");
   		cereal::BinaryInputArchive contigTableArchive(contigTableStream);
   		contigTableArchive(refNames_);
    		// contigTableArchive(cPosInfo_);
    		contigTableArchive(contigTable_);
    		contigTableStream.close();

		int* ptr1 = (int*)shmat(shmid1,NULL,0);
		*ptr1++ = contigTable_.size();
		key_t keyTemp = (key_t)hash_string(indexDir + "/qqqctable.bin" + "refLengths_"+ "vectorElements");
		int shmidTemp =  shmget(keyTemp, size, IPC_CREAT | 0666);
		util::Position* ptrTemp = (util::Position*)shmat(shmidTemp,NULL,0);

		for(int i = 0; i < contigTable_.size(); i++)
		{

			*ptr1++ = contigTable_[i].size();
			
			
			for(int k=0; k < contigTable_[i].size(); k++)
			{
				*ptrTemp++ = contigTable_[i][k];
			}
		}

		
	}
	else
	{
		std::cerr << "contigTable_ ELSE \n";

		CLI::AutoTimer timer{"Loading contig table", CLI::Timer::Big};
    		std::ifstream contigTableStream(indexDir + "/ctable.bin");
   		cereal::BinaryInputArchive contigTableArchive(contigTableStream);
   		contigTableArchive(refNames_);
    		// contigTableArchive(cPosInfo_);
    		contigTableArchive(contigTable_);
    		contigTableStream.close();
		shmid1 = shmget(key1, size/100, IPC_CREAT | 0666);
		int* ptr1 = (int*)shmat(shmid1,NULL,0);
		int ksize = *ptr1++;
		contigTable_.resize(ksize);
		key_t keyTemp = (key_t)hash_string(indexDir + "/qqqctable.bin" + "refLengths_"+ "vectorElements");
		int shmidTemp =  shmget(keyTemp, size, IPC_CREAT | 0666);
		util::Position* ptrTemp = (util::Position*)shmat(shmidTemp,NULL,0);
		for(int i =0 ; i < ksize; i++)
		{
			//std::cerr << *ptr1 << " ";
			int tempSize = *ptr1++;
			contigTable_[i].resize(tempSize);
			//std::cerr<<"CHK1 ";
			for(int k=0; k <tempSize; k++)
			{
				//std::cerr<<"CHK2 ";
				contigTable_[i][k] = (*ptrTemp++);
			}
			
		}
		
	}




    
  }
  numContigs_ = contigTable_.size();

  {
    {
	//std::cerr << "chk1 ";
    std::ifstream file2(indexDir + "/reflengths.bin", std::ios::binary | std::ios::ate);
    key_t key = (key_t)hash_string(indexDir + "/reflengths.bin" + "refLengths_");
    std::ifstream::pos_type size = file2.tellg();
    size = size * 10000;
    int shmid2 = shmget(key, size, IPC_CREAT | IPC_EXCL);
    key_t key1 = (key_t)hash_string(indexDir + "/reflengths.bin" + "refLengths_" + "vectorSize");
    int shmid3 = shmget(key1, size, IPC_CREAT | IPC_EXCL);
	//std::cerr << "chk2 ";
    if(shmid3 != -1)
    {
	//std::cerr << "chk3 ";
	shmctl(shmid2, IPC_RMID, NULL);
	shmid2 = shmget(key, size, IPC_CREAT | 0666);
	shmctl(shmid3, IPC_RMID, NULL);
	shmid3 = shmget(key1, size, IPC_CREAT | 0666);
  	uint32_t* ptr2 = (uint32_t*)shmat(shmid2,NULL,0);
	int* ptr3 = (int*)shmat(shmid3,NULL,0);
 	    std::string rlPath = indexDir + "/reflengths.bin";
	    if (puffer::fs::FileExists(rlPath.c_str())) {
	      CLI::AutoTimer timer{"Loading reference lengths", CLI::Timer::Big};
	      std::ifstream refLengthStream(rlPath);
	      cereal::BinaryInputArchive refLengthArchive(refLengthStream);
	      refLengthArchive(refLengths_);
	    } else {
	      refLengths_ = std::vector<uint32_t>(refNames_.size(), 1000);
	    }
	//std::cerr << "chk4 ";
	   refLengths_pointer = &refLengths_;
	//std::cerr << "chk5 ";
	*ptr3 = refLengths_.size();
	//std::cerr << "chk6 ";
	for(int i = 0; i < (int)refLengths_.size(); i++)
	{
		//std::cerr << "chk7 ";
		*ptr2++ = refLengths_[i];
	}
	

    }
    else
    {
	//std::cerr << "chk8 ";
	int shmid2 = shmget(key, size, IPC_CREAT | 0666);
	int shmid3 = shmget(key1, size, IPC_CREAT | 0666);
	int* ptr3 = (int*)shmat(shmid3,NULL,0);
	//std::cerr << "chk9 ";        
	uint32_t* ptr2 = (uint32_t*)shmat(shmid2,NULL,0);
	std::vector<uint32_t> refLengths_1(*ptr3);
	//std::cerr << "chk10 ";
	refLengths_.resize(*ptr3);
	for(int i = 0; i < (*ptr3); i++)
	{
		//std::cerr << "chk11 ";
		refLengths_[i] = *ptr2++;	
	}
	//std::cerr << "chk12 ";
    }
    
}
  }

  {
    CLI::AutoTimer timer{"Loading eq table", CLI::Timer::Big};
    std::ifstream eqTableStream(indexDir + "/eqtable.bin");
    cereal::BinaryInputArchive eqTableArchive(eqTableStream);
    eqTableArchive(eqClassIDs_);
    eqTableArchive(eqLabels_);
    eqTableStream.close();
  }

  {
    CLI::AutoTimer timer{"Loading mphf table", CLI::Timer::Big};
    std::string hfile = indexDir + "/mphf.bin";
    std::ifstream hstream(hfile);
    hash_.reset(new boophf_t);
    hash_->load(hstream);
    hstream.close();
    hash_raw_ = hash_.get();
  }

  {
    CLI::AutoTimer timer{"Loading contig boundaries", CLI::Timer::Big};
    std::string bfile = indexDir + "/rank.bin";
    sdsl::load_from_file(contigBoundary_, bfile);
    contigRank_ = decltype(contigBoundary_)::rank_1_type(&contigBoundary_);
    contigSelect_ = decltype(contigBoundary_)::select_1_type(&contigBoundary_);
  }
  /*
  selectPrecomp_.reserve(numContigs_+1);
  selectPrecomp_.push_back(0);
  for (size_t i = 1; i < numContigs_; ++i) {
    selectPrecomp_.push_back(contigSelect_(i));
  }
  selectPrecomp_.push_back(contigSelect_(numContigs_));
  */

  {
    CLI::AutoTimer timer{"Loading sequence", CLI::Timer::Big};
    std::string sfile = indexDir + "/seq.bin";
    sdsl::load_from_file(seq_, sfile);
    lastSeqPos_ = seq_.size() - k_;
  }

  {
    CLI::AutoTimer timer{"Loading positions", CLI::Timer::Big};
    std::string pfile = indexDir + "/pos.bin";
    sdsl::load_from_file(pos_, pfile);
  }

  {
    CLI::AutoTimer timer{"Loading edges", CLI::Timer::Big};
    std::string pfile = indexDir + "/edge.bin";
    sdsl::load_from_file(edge_, pfile);
  }
  /*
  {
    CLI::AutoTimer timer{"Loading edges", CLI::Timer::Big};
    std::string pfile = indexDir + "/revedge.bin";
    sdsl::load_from_file(revedge_, pfile);
  }
  */
}

PufferfishIndex::EqClassID PufferfishIndex::getEqClassID(uint32_t contigID) {
  return eqClassIDs_[contigID];
}

const PufferfishIndex::EqClassLabel&
PufferfishIndex::getEqClassLabel(uint32_t contigID) {
  return eqLabels_[getEqClassID(contigID)];
}

uint64_t PufferfishIndex::getRawPos(CanonicalKmer& mer) {
  auto km = mer.getCanonicalWord();
  size_t res = hash_raw_->lookup(km);
  if (res < numKmers_) {
    uint64_t pos = pos_[res];
    uint64_t fk = seq_.get_int(2 * pos, twok_);
    auto keq = mer.isEquivalent(fk);
    if (keq != KmerMatchType::NO_MATCH) {
      return pos;
    }
  }
  return std::numeric_limits<uint64_t>::max();
}

bool PufferfishIndex::contains(CanonicalKmer& mer) {
  return isValidPos(getRawPos(mer));
}

bool PufferfishIndex::isValidPos(uint64_t pos) {
  return pos != std::numeric_limits<uint64_t>::max();
}

uint32_t PufferfishIndex::contigID(CanonicalKmer& mer) {
  auto km = mer.getCanonicalWord();
  size_t res = hash_raw_->lookup(km);
  if (res < numKmers_) {
    uint64_t pos = pos_[res];
    uint64_t fk = seq_.get_int(2 * pos, twok_);
    auto keq = mer.isEquivalent(fk);
    if (keq != KmerMatchType::NO_MATCH) {
      auto rank = contigRank_(pos);
      return rank;
    }
  }
  return std::numeric_limits<uint32_t>::max();
}


std::string PufferfishIndex::getSeqStr(size_t globalPos, size_t length, bool isFw) {
	std::string outstr = "";
	uint64_t validLength = 0;
	uint64_t word = 0;
	uint8_t base = 0;
	while (length > 0) {
	validLength = std::min(length, (size_t)32);
	length -= validLength;
 	word = seq_.get_int(2*globalPos, 2*validLength);
	globalPos += validLength;
	if (isFw)
	  for(uint64_t i = 0; i < 2*validLength ;i+=2){
      base = (word >> i) & 0x03;
	    switch(base){
    	case 0:
        outstr += 'A';
        break ;
	    case 1:
        outstr += 'C';
	      break ;
    	case 2:
        outstr += 'G';
    	  break ;
	    case 3:
        outstr += 'T';
	      break ;
    	}
    }
	else
    for(uint64_t i = 0; i < 2*validLength ;i+=2){
      base = (word >> i) & 0x03;
      switch(base){
      case 0:
        outstr = 'T' + outstr;
        break ;
      case 1:
        outstr = 'G' + outstr;
        break ;
      case 2:
        outstr = 'C' + outstr;
        break ;
      case 3:
        outstr = 'A' + outstr;
        break ;
      }
    }
	}
  return outstr;
}



inline uint64_t our_read_int(uint64_t* word, uint8_t offset, const uint8_t len)
{
    uint64_t w1 = (*word)>>offset;
    if ((offset+len) > 64) { // if offset+len > 64
        return w1 |  // w1 or w2 adepted:
               ((*(word+1) & sdsl::bits::lo_set[(offset+len)&0x3F])   // set higher bits zero
                << (64-offset));  // move bits to the left
    } else {
        return w1 & sdsl::bits::lo_set[len];
    }
}

/**
 * Returns a ProjectedHits object containing all of the reference loci matching this
 * provided Canonical kmer (including the oritentation of the match).  The provided
 * QueryCache argument will be used to avoid redundant rank / select operations if feasible.
 */
auto PufferfishIndex::getRefPos(CanonicalKmer& mer, util::QueryCache& qc)
    -> util::ProjectedHits {
  using IterT = std::vector<util::Position>::iterator;
  auto km = mer.getCanonicalWord();
  size_t res = hash_raw_->lookup(km);
  if (res < numKmers_) {
    uint64_t pos = pos_[res];
    // if using quasi-dictionary idea (https://arxiv.org/pdf/1703.00667.pdf)
    /* 
    uint64_t hashbits = pos & 0xF;
    pos = pos >> 4;
    if ((km & 0xF) != hashbits) { 
        return {std::numeric_limits<uint32_t>::max(),
          std::numeric_limits<uint64_t>::max(),
          std::numeric_limits<uint32_t>::max(),
          true,
          0,
          k_,
          core::range<IterT>{}};
    }
    */
    uint64_t twopos = pos << 1;
    uint64_t fk = seq_.get_int(twopos, twok_);
    // say how the kmer fk matches mer; either
    // identity, twin (i.e. rev-comp), or no match
    auto keq = mer.isEquivalent(fk);
    if (keq != KmerMatchType::NO_MATCH) {
      // the index of this contig
      auto rank = contigRank_(pos);
      // the reference information in the contig table
      auto& pvec = contigTable_[rank];
      // start position of this contig
      uint64_t sp = 0;
      uint64_t contigEnd = 0;
      if (rank == qc.prevRank) {
        sp = qc.contigStart;
        contigEnd = qc.contigEnd;
      } else {
        sp = (rank == 0) ? 0 : static_cast<uint64_t>(contigSelect_(rank)) + 1;
        contigEnd = contigSelect_(rank + 1);
        qc.prevRank = rank;
        qc.contigStart = sp;
        qc.contigEnd = contigEnd;
      }

      // relative offset of this k-mer in the contig
      uint32_t relPos = static_cast<uint32_t>(pos - sp);

      // start position of the next contig - start position of this one
      auto clen = static_cast<uint64_t>(contigEnd + 1 - sp);
      // auto clen =
      // cPosInfo_[rank].length();//static_cast<uint64_t>(contigSelect_(rank +
      // 1) + 1 - sp);

      // how the k-mer hits the contig (true if k-mer in fwd orientation, false
      // otherwise)
      bool hitFW = (keq == KmerMatchType::IDENTITY_MATCH);
      return {static_cast<uint32_t>(rank),
              pos,
              relPos,
              hitFW,
              static_cast<uint32_t>(clen),
              k_,
              core::range<IterT>{pvec.begin(), pvec.end()}};
    } else {
      return {std::numeric_limits<uint32_t>::max(),
              std::numeric_limits<uint64_t>::max(),
              std::numeric_limits<uint32_t>::max(),
              true,
              0,
              k_,
              core::range<IterT>{}};
    }
  }

  return {std::numeric_limits<uint32_t>::max(),
          std::numeric_limits<uint64_t>::max(),
          std::numeric_limits<uint32_t>::max(),
          true,
          0,
          k_,
          core::range<IterT>{}};
}

auto PufferfishIndex::getRefPos(CanonicalKmer& mer) -> util::ProjectedHits {
  using IterT = std::vector<util::Position>::iterator;
  auto km = mer.getCanonicalWord();
  size_t res = hash_raw_->lookup(km);
  if (res < numKmers_) {
    uint64_t pos = pos_[res];
    uint64_t twopos = pos << 1;
    uint64_t fk = seq_.get_int(twopos, twok_);
    // say how the kmer fk matches mer; either
    // identity, twin (i.e. rev-comp), or no match
    auto keq = mer.isEquivalent(fk);
    if (keq != KmerMatchType::NO_MATCH) {
      // the index of this contig
      auto rank = contigRank_(pos);
      // the reference information in the contig table
      auto& pvec = contigTable_[rank];
      // start position of this contig
      uint64_t sp =
          (rank == 0) ? 0 : static_cast<uint64_t>(contigSelect_(rank)) + 1;
      uint64_t contigEnd = contigSelect_(rank + 1);

      // relative offset of this k-mer in the contig
      uint32_t relPos = static_cast<uint32_t>(pos - sp);

      // start position of the next contig - start position of this one
      auto clen = static_cast<uint64_t>(contigEnd + 1 - sp);
      // auto clen =
      // cPosInfo_[rank].length();//static_cast<uint64_t>(contigSelect_(rank +
      // 1) + 1 - sp);

      // how the k-mer hits the contig (true if k-mer in fwd orientation, false
      // otherwise)
      bool hitFW = (keq == KmerMatchType::IDENTITY_MATCH);
      return {static_cast<uint32_t>(rank),
              pos,
              relPos,
              hitFW,
              static_cast<uint32_t>(clen),
              k_,
              core::range<IterT>{pvec.begin(), pvec.end()}};
    } else {
      return {std::numeric_limits<uint32_t>::max(),
              std::numeric_limits<uint64_t>::max(),
              std::numeric_limits<uint32_t>::max(),
              true,
              0,
              k_,
              core::range<IterT>{}};
    }
  }

  return {std::numeric_limits<uint32_t>::max(),
          std::numeric_limits<uint64_t>::max(),
          std::numeric_limits<uint32_t>::max(),
          true,
          0,
          k_,
          core::range<IterT>{}};
}

uint32_t PufferfishIndex::k() { return k_; }

CanonicalKmer PufferfishIndex::getStartKmer(uint64_t rank){
  CanonicalKmer::k(k_) ;
  CanonicalKmer kb ;
  uint64_t sp = (rank == 0) ? 0 : static_cast<uint64_t>(contigSelect_(rank)) + 1;
  uint64_t fk = seq_.get_int(2*sp, 2*k_) ;
  kb.fromNum(fk) ;
  return kb ;

}
CanonicalKmer PufferfishIndex::getEndKmer(uint64_t rank){
  CanonicalKmer::k(k_) ;
  CanonicalKmer kb ;
  //uint64_t sp = (rank == 0) ? 0 : static_cast<uint64_t>(contigSelect_(rank)) + 1;
  uint64_t contigEnd = contigSelect_(rank + 1);

  uint64_t fk = seq_.get_int(2*(contigEnd - k_ + 1), 2*k_) ;
  kb.fromNum(fk) ;
  return kb ;
}


std::vector<CanonicalKmer> PufferfishIndex::getNextKmerOnGraph(uint64_t rank, util::Direction dir, bool isCurContigFwd){
  //get the edge vec
  std::vector<CanonicalKmer> nextKmers ;
  uint8_t edgeVec = edge_[rank] ;
  uint8_t mask = 1 ;
  std::vector<char> nuclmap = {'C','G','T','A','C','G','T','A'} ;
  std::map<char, char> cMap = {{'A','T'}, {'T','A'}, {'C','G'}, {'G','C'}} ;

  if(dir == util::Direction::FORWARD){
    // We need to append so let's concentrate on the lower 4 bits
    auto ke = getEndKmer(rank) ;
    auto ktmp = ke ;
    for(uint8_t i=0; i < 4; ++i){
      ktmp = ke ;
      if(edgeVec & (mask << i)){
        char c = nuclmap[i] ;
        char charToAdd = (isCurContigFwd) ? c : cMap[c] ;
        ktmp.shiftFw(charToAdd) ;
        nextKmers.push_back(ktmp) ;
      }
    }
  }else{
    auto kb = getStartKmer(rank) ;
    auto ktmp = kb ;
    for(uint8_t i=4; i < 8; ++i){
      ktmp = kb ;
      if(edgeVec & (mask << i)){
        char c = nuclmap[i] ;
        char charToAdd = (isCurContigFwd) ? c : cMap[c] ;
        ktmp.shiftBw(charToAdd) ;
        nextKmers.push_back(ktmp) ;
      }
    }
  }
  return nextKmers ;
}

uint32_t PufferfishIndex::getContigLen(uint64_t rank){
  uint64_t sp = (rank == 0) ? 0 : static_cast<uint64_t>(contigSelect_(rank)) + 1;
  uint64_t contigEnd = contigSelect_(rank + 1);
  return (static_cast<uint32_t>(contigEnd - sp + 1)) ;
}

uint64_t PufferfishIndex::getGlobalPos(uint64_t rank){
  uint64_t sp = (rank == 0) ? 0 : static_cast<uint64_t>(contigSelect_(rank)) + 1;
  return sp ;
}

auto  PufferfishIndex::getContigBlock(uint64_t rank)->util::ContigBlock{
  CanonicalKmer::k(k_) ;
  CanonicalKmer kb;
  CanonicalKmer ke;
  uint64_t sp = (rank == 0) ? 0 : static_cast<uint64_t>(contigSelect_(rank)) + 1;
  uint64_t contigEnd = contigSelect_(rank+1) ;

  uint32_t clen = static_cast<uint32_t>(contigEnd - sp + 1) ;
  uint64_t fk = seq_.get_int(2*sp, 2*k_) ;
  kb.fromNum(fk) ;

  fk = seq_.get_int(2*(contigEnd - k_ + 1), 2*k_) ;
  ke.fromNum(fk) ;

  std::string seq = getSeqStr(sp,clen) ;

  return {rank,sp,clen,seq};
}

/**
 * Return the position list (ref_id, pos) corresponding to a contig.
 */
const std::vector<util::Position>&
PufferfishIndex::refList(uint64_t contigRank) {
  return contigTable_[contigRank];
}

const std::string& PufferfishIndex::refName(uint64_t refRank) {
  return refNames_[refRank];
}

uint32_t PufferfishIndex::refLength(uint64_t refRank) const {
	//std::cerr << "chkpt 15 ";
  return refLengths_[refRank];
}

const std::vector<std::string>& PufferfishIndex::getRefNames() {
  return refNames_;
}

const std::vector<uint32_t>& PufferfishIndex::getRefLengths() const {
	//std::cerr <<"chkpt 16 ";
  return refLengths_;
}

