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


class membuf : public std::basic_streambuf<char> {
public:
  membuf(const char *p, size_t l) {
    setg((char*)p, (char*)p, (char*)p + l);
  }
};


class memstream : public std::istream {
public:
  memstream(const char *p, size_t l) :
    std::istream(&_buffer),
    _buffer(p, l) {
    rdbuf(&_buffer);
  }

private:
  membuf _buffer;
};



std::vector<std::vector<util::Position>>* ptr1 = NULL;
int shmid1;
int shmid_2;
std::vector<uint32_t>* ptr_2 = NULL;


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
  int shmid;
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
  } 

  {
    std::ifstream file1(indexDir + "/ctable.bin", std::ios::binary | std::ios::ate);
    std::string tempFile = indexDir + "/ctable.bin";
    FILE * pFile = fopen (tempFile.c_str(), "rb");
    
    key_t key = (key_t)hash_string(indexDir + "/ctable.bin" + "contigTable_");
    std::ifstream::pos_type size = file1.tellg();
    shmid1 = shmget(key, size, IPC_CREAT | IPC_EXCL);
    std::cerr << "KEY :" << key << "\n";
    std::cerr << "SHMID : " << shmid1 << "\n";
    if(shmid1 != -1)
    {
	std::cerr << "SHMID: " << shmid1 << "\n";
	std::cerr << "Step If 1";
        CLI::AutoTimer timer{"Loading contig table", CLI::Timer::Big};
	std::cerr << "Bookmark 1 \n";
        std::ifstream contigTableStream(indexDir + "/ctable.bin");
	std::cerr << "Bookmark 2\n";
        cereal::BinaryInputArchive contigTableArchive(contigTableStream);
	std::cerr << "Bookmark 3\n";        
	contigTableArchive(refNames_);
	std::cerr << "Bookmark 4\n";
        // contigTableArchive(cPosInfo_);
        contigTableArchive(contigTable_);
	std::cerr << "Bookmark 5\n";
        contigTableStream.close();
	std::cerr << "Bookmark 6\n";
        //std::vector<std::vector<util::Position>>* ptr = NULL;
	shmctl(shmid1, IPC_RMID, NULL);
	shmid1 = shmget(key, size, IPC_CREAT | 0666);
	//ptr1 = (std::vector<std::vector<util::Position>>*)shmat(shmid1,NULL,0);
	char* ptr_ch = (char*)shmat(shmid1,NULL,0);
	char * memblock =  new char[size];
	file1.seekg (0, std::ios::beg);
	file1.read (memblock, size);
	
        
	if(ptr_ch == NULL)
	{
		std::cerr << "POINTER IS NULL\n";		
	}
	else
	{
		std::cerr << "POINTER IS NOT NULL\n";	
	}        
	std::cerr << "Bookmark 7\n";
	std::cerr << "COUNT : " << size << "\n";
	/*for(int i = 0; i < file1.tellg(); i++)
	{
		*ptr_ch++ = memblock[i];


	}*/
	fread(ptr_ch, 1, file1.tellg(), pFile);
	std::cerr << "SHMIDIDID111: " << shmid1 << "\n";	
	std::cerr << "KEY1: " << key << "\n";
	//std::cerr << "POINTER1: " << ptr1 << "\n";
	//*ptr1 = contigTable_;
	//std::cerr << "POINTER11: " << ptr1 << "\n";
	std::cerr << "Bookmark 8\n";
	std::cerr << "Step If End 1\n";
	
		
    }
    else
    {
	std::cerr << "Entered ELSE 1\n";
	std::cerr << "Step Else 1\n";
	if(errno == ENOMEM) 
	{
		std::cerr << "ENOMEM" << "\n";
	}
	else if(errno == EACCES) 
	{
		std::cerr << "EACCES" << "\n";
	}
	else if(errno == EEXIST) 
	{
		std::cerr << "EEXIST" << "\n";
	}
	else if(errno == EINVAL) 
	{
		std::cerr << "EINVAL" << "\n";
	}
	else if(errno == ENFILE) 
	{
		std::cerr << "ENFILE" << "\n";
	}
	else if(errno == ENOENT) 
	{
		std::cerr << "ENOENT" << "\n";
	}
	else if(errno == ENOSPC) 
	{
		std::cerr << "ENOSPC" << "\n";
	}
	else if(errno == EPERM) 
	{
		std::cerr << "EPERM" << "\n";
	}
	else  
	{
		std::cerr << "EA SOME OTHER ERROR" << "\n";
	}
	shmid1 = shmget(key, size, IPC_CREAT | 0666);


        char* ptr_ch = (char*)shmat(shmid1,NULL,0);




	std::cerr << "KEY11 :" << key << "\n";
        std::cerr << "SHMID11 : " << shmid1 << "\n";
	std::cerr << "POINTER1 : " << ptr_ch;
	if(ptr_ch == NULL)
	{
		std::cerr << "POINTER IS NULL\n";		
	}
	else
	{
		std::cerr << "POINTER IS NOT NULL\n";	
	}      
	//std::vector<std::vector<util::Position>> *temp = new std::vector<std::vector<util::Position>>();
        //temp  = ptr1;
	//unsigned int i,j,k;
	//std::cerr << "Check1 : \n";
	//contigTable_ = new std::vector<std::vector<util::Position>>(*ptr1);
	/*for(i = 0; i < ((*ptr1).size()); i++)
	{
		//std::cerr << "Check2 : \n";
		std::vector<util::Position> row;
		contigTable_.push_back(row);	
	}
	for(j = 0; j < ((*ptr1).size());  j++)
	{
		std::cerr << "Check3 : \n";
		for(k = 0; k < (*ptr1)[j].size(); k++)
		{
			std::cerr << "Check4 : \n";
			contigTable_[j].push_back((*ptr1)[j][k]);	
		}	
	}*/


        //CLI::AutoTimer timer{"Loading contig table", CLI::Timer::Big};
	
	std::cerr << "ELSE 1 : " << "memstream start \n";

        memstream contigTableStream(ptr_ch, size);
	
	std::cerr << "ELSE 1 : " << "memstream end \n";

	std::cerr << "ELSE 1 : " << "BinaryInputArchive start \n";
	
        cereal::BinaryInputArchive contigTableArchive(contigTableStream);

	std::cerr << "ELSE 1 : " << "BinaryInputArchive end \n";

	
	contigTableArchive(refNames_);

	std::cerr << "ELSE 1 : " << "refNames_ \n";


        // contigTableArchive(cPosInfo_);
        contigTableArchive(contigTable_);
	
	std::cerr << "ELSE 1 : " << "contigTableArchive \n";

        //contigTableStream.close();
	
        std::cerr << "Aabra ka dabra \n";
    }
  }
  //numContigs_ = ptr1->size();
  //numContigs_ = contigTable_.size();
  std::cerr << "CONTIG TABLE SIZE : " << contigTable_.size() << "\n";
  std::cerr << "Aabra ka dabra11111 \n";
  //std::cerr << "SIZE : " <<  ptr1->size() << "\n";

  {
 	std::ifstream file2(indexDir + "/reflengths.bin", std::ios::binary | std::ios::ate);
	key_t key = (key_t)hash_string(indexDir + "/reflengths.bin" + "refLengths_");
	std::ifstream::pos_type size = file2.tellg();
	int shmid2 = shmget(key, size, IPC_CREAT | IPC_EXCL);


	  if(shmid2 != -1)
	  {

	    std::cerr << "IF2 \n";
	    std::string rlPath = indexDir + "/reflengths.bin";
	    if (puffer::fs::FileExists(rlPath.c_str())) {
	      CLI::AutoTimer timer{"Loading reference lengths", CLI::Timer::Big};
	      std::ifstream refLengthStream(rlPath);
	      cereal::BinaryInputArchive refLengthArchive(refLengthStream);
	      refLengthArchive(refLengths_);
	    } else {
	      refLengths_ = std::vector<uint32_t>(refNames_.size(), 1000);
	    }
		shmctl(shmid2, IPC_RMID, NULL);
		shmid2 = shmget(key, size, IPC_CREAT | 0666);
		char* ptr_ch = (char*)shmat(shmid2,NULL,0);
		char * memblock =  new char[size];
		file2.seekg (0, std::ios::beg);
		file2.read (memblock, size);
	
		if(ptr_ch == NULL)
		{
			std::cerr << "POINTER IS NULL\n";		
		}
		else
		{
			std::cerr << "POINTER IS NOT NULL\n";	
		}        
		for(int i = 0; i < file2.tellg(); i++)
		{
			*ptr_ch++ = memblock[i];
		}
	
	    }



	  else
	  {
		std::cerr << "ELSE2 \n";
		shmid2 = shmget(key, size, IPC_CREAT | 0666);
		char* ptr_ch = (char*)shmat(shmid2,NULL,0);
		if(ptr_ch == NULL)
		{
			std::cerr << "POINTER IS NULL\n";		
		}
		else
		{
			std::cerr << "POINTER IS NOT NULL\n";	
		}      
		memstream refLengthStream(ptr_ch, size);
		cereal::BinaryInputArchive refLengthArchive(refLengthStream);
		refLengthArchive(refLengths_);
	
	  }
  }

  {
	std::ifstream file3(indexDir + "/eqtable.bin", std::ios::binary | std::ios::ate);
	key_t key = (key_t)hash_string(indexDir + "/eqtable.bin" + "eqTableStream");
	std::ifstream::pos_type size = file3.tellg();
	int shmid3 = shmget(key, size, IPC_CREAT | IPC_EXCL);


	if(shmid3 != -1)
	{
		std::cerr << "IF 3 \n";
	    	CLI::AutoTimer timer{"Loading eq table", CLI::Timer::Big};
	    	std::ifstream eqTableStream(indexDir + "/eqtable.bin");
	    	cereal::BinaryInputArchive eqTableArchive(eqTableStream);
	    	eqTableArchive(eqClassIDs_);
	    	eqTableArchive(eqLabels_);
	    	eqTableStream.close();
		shmctl(shmid3, IPC_RMID, NULL);
		shmid3 = shmget(key, size, IPC_CREAT | 0666);
		char* ptr_ch = (char*)shmat(shmid3,NULL,0);
		char * memblock =  new char[size];
		file3.seekg (0, std::ios::beg);
		file3.read (memblock, size);
	
		if(ptr_ch == NULL)
		{
			std::cerr << "POINTER IS NULL\n";		
		}
		else
		{
			std::cerr << "POINTER IS NOT NULL\n";	
		}        
		for(int i = 0; i < file3.tellg(); i++)
		{
			*ptr_ch++ = memblock[i];
		}
	
	}



	else
	{
		std::cerr << "ELSE 3 \n";
		shmid3 = shmget(key, size, IPC_CREAT | 0666);
		char* ptr_ch = (char*)shmat(shmid3,NULL,0);
		if(ptr_ch == NULL)
		{
			std::cerr << "POINTER IS NULL\n";		
		}
		else
		{
			std::cerr << "POINTER IS NOT NULL\n";	
		}      
		memstream eqTableStream(ptr_ch, size);
		cereal::BinaryInputArchive eqTableArchive(eqTableStream);
	    	eqTableArchive(eqClassIDs_);
	    	eqTableArchive(eqLabels_);
	
	}
  }

  {
    	std::ifstream file4(indexDir + "/mphf.bin", std::ios::binary | std::ios::ate);
	key_t key = (key_t)hash_string(indexDir + "/mphf.bin" + "hstream");
	std::ifstream::pos_type size = file4.tellg();
	int shmid4 = shmget(key, size, IPC_CREAT | IPC_EXCL);


	if(shmid4 != -1)
	{
	    std::cerr << "IF 4 \n";
	    CLI::AutoTimer timer{"Loading mphf table", CLI::Timer::Big};
	    std::string hfile = indexDir + "/mphf.bin";
	    std::ifstream hstream(hfile);
	    hash_.reset(new boophf_t);
	    hash_->load(hstream);
	    hstream.close();
	    hash_raw_ = hash_.get();
		shmctl(shmid4, IPC_RMID, NULL);
		shmid4 = shmget(key, size, IPC_CREAT | 0666);
		char* ptr_ch = (char*)shmat(shmid4,NULL,0);
		char * memblock =  new char[size];
		file4.seekg (0, std::ios::beg);
		file4.read (memblock, size);
		
		if(ptr_ch == NULL)
		{
			std::cerr << "POINTER IS NULL\n";		
		}
		else
		{
			std::cerr << "POINTER IS NOT NULL\n";	
		}        
		for(int i = 0; i < file4.tellg(); i++)
		{
			*ptr_ch++ = memblock[i];
		}
		
	}



	else
	{
		std::cerr << "ELSE 4 \n";
		shmid4 = shmget(key, size, IPC_CREAT | 0666);
		char* ptr_ch = (char*)shmat(shmid4,NULL,0);
		if(ptr_ch == NULL)
		{
			std::cerr << "POINTER IS NULL\n";		
		}
		else
		{
			std::cerr << "POINTER IS NOT NULL\n";	
		}      
		memstream hstream(ptr_ch, size);
	    hash_.reset(new boophf_t);
	    hash_->load(hstream);
	    hash_raw_ = hash_.get();
		
	}
  }

    {
    		std::ifstream file5(indexDir + "/rank.bin", std::ios::binary | std::ios::ate);
		key_t key = (key_t)hash_string(indexDir + "/rank.bin" + "contigBoundary_");
		std::ifstream::pos_type size = file5.tellg();
		int shmid5 = shmget(key, size, IPC_CREAT | IPC_EXCL);


	if(shmid5 != -1)
	{
	    std::cerr << "IF 5 \n";
	    CLI::AutoTimer timer{"Loading contig boundaries", CLI::Timer::Big};
    	    std::string bfile = indexDir + "/rank.bin";
    	    sdsl::load_from_file(contigBoundary_, bfile);
    	    contigRank_ = decltype(contigBoundary_)::rank_1_type(&contigBoundary_);
    	    contigSelect_ = decltype(contigBoundary_)::select_1_type(&contigBoundary_);
		shmctl(shmid5, IPC_RMID, NULL);
		shmid5 = shmget(key, size, IPC_CREAT | 0666);
		char* ptr_ch = (char*)shmat(shmid5,NULL,0);
		char * memblock =  new char[size];
		file5.seekg (0, std::ios::beg);
		file5.read (memblock, size);
		
		if(ptr_ch == NULL)
		{
			std::cerr << "POINTER IS NULL\n";		
		}
		else
		{
			std::cerr << "POINTER IS NOT NULL\n";	
		}        
		for(int i = 0; i < file5.tellg(); i++)
		{
			*ptr_ch++ = memblock[i];
		}
		std::cerr << "SHMID5 : " << shmid5 << "\n";
		
	}



	else
	{
		std::cerr << "ELSE 5 \n";
		shmid5 = shmget(key, size, IPC_CREAT | 0666);
		char* ptr_ch = (char*)shmat(shmid5,NULL,0);
		if(ptr_ch == NULL)
		{
			std::cerr << "POINTER IS NULL\n";		
		}
		else
		{
			std::cerr << "POINTER IS NOT NULL\n";	
		}      
		memstream rstream(ptr_ch, size);
		contigBoundary_.load(rstream);
    	contigRank_ = decltype(contigBoundary_)::rank_1_type(&contigBoundary_);
    	contigSelect_ = decltype(contigBoundary_)::select_1_type(&contigBoundary_);
	std::cerr << "SHMID5 : " << shmid5 << "\n";
		
	}
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
    	std::ifstream file6(indexDir + "/seq.bin", std::ios::binary | std::ios::ate);
		key_t key = (key_t)hash_string(indexDir + "/seq.bin" + "seq_");
		std::ifstream::pos_type size = file6.tellg();
		int shmid6 = shmget(key, size, IPC_CREAT | IPC_EXCL);


	if(shmid6 != -1)
	{
	    std::cerr << "IF 6 \n";
	    CLI::AutoTimer timer{"Loading sequence", CLI::Timer::Big};
    	std::string sfile = indexDir + "/seq.bin";
    	sdsl::load_from_file(seq_, sfile);
    	lastSeqPos_ = seq_.size() - k_;
		shmctl(shmid6, IPC_RMID, NULL);
		shmid6 = shmget(key, size, IPC_CREAT | 0666);
		char* ptr_ch = (char*)shmat(shmid6,NULL,0);
		char * memblock =  new char[size];
		file6.seekg (0, std::ios::beg);
		file6.read (memblock, size);
		
		if(ptr_ch == NULL)
		{
			std::cerr << "POINTER IS NULL\n";		
		}
		else
		{
			std::cerr << "POINTER IS NOT NULL\n";	
		}        
		for(int i = 0; i < file6.tellg(); i++)
		{
			*ptr_ch++ = memblock[i];
		}
		
	}



	else
	{
		std::cerr << "ELSE 6 \n";
		shmid6 = shmget(key, size, IPC_CREAT | 0666);
		char* ptr_ch = (char*)shmat(shmid6,NULL,0);
		if(ptr_ch == NULL)
		{
			std::cerr << "POINTER IS NULL\n";		
		}
		else
		{
			std::cerr << "POINTER IS NOT NULL\n";	
		}      
		memstream seqstream(ptr_ch, size);
		seq_.load(seqstream);
    	lastSeqPos_ = seq_.size() - k_;
		
	}
  }

    {
    		std::ifstream file7(indexDir + "/pos.bin", std::ios::binary | std::ios::ate);
		key_t key = (key_t)hash_string(indexDir + "/pos.bin" + "pos_");
		std::ifstream::pos_type size = file7.tellg();
		int shmid7 = shmget(key, size, IPC_CREAT | IPC_EXCL);


	if(shmid7 != -1)
	{
	    std::cerr << "IF 7 \n";
	    CLI::AutoTimer timer{"Loading positions", CLI::Timer::Big};
    	std::string pfile = indexDir + "/pos.bin";
    	sdsl::load_from_file(pos_, pfile);
		shmctl(shmid7, IPC_RMID, NULL);
		shmid7 = shmget(key, size, IPC_CREAT | 0666);
		char* ptr_ch = (char*)shmat(shmid7,NULL,0);
		char * memblock =  new char[size];
		file7.seekg (0, std::ios::beg);
		file7.read (memblock, size);
		
		if(ptr_ch == NULL)
		{
			std::cerr << "POINTER IS NULL\n";		
		}
		else
		{
			std::cerr << "POINTER IS NOT NULL\n";	
		}
		//fread(ptr_ch, 1, file7.tellg(), file7);        
		for(int i = 0; i < file7.tellg(); i++)
		{
			*ptr_ch++ = memblock[i];
		}
		
	}



	else
	{
		std::cerr << "ELSE 7 \n";
		shmid7 = shmget(key, size, IPC_CREAT | 0666);
		char* ptr_ch = (char*)shmat(shmid7,NULL,0);
		if(ptr_ch == NULL)
		{
			std::cerr << "POINTER IS NULL\n";		
		}
		else
		{
			std::cerr << "POINTER IS NOT NULL\n";	
		}      
		memstream posstream(ptr_ch, size);
		pos_.load(posstream);   	
	}
  }

    {
    	std::ifstream file8(indexDir + "/edge.bin", std::ios::binary | std::ios::ate);
		key_t key = (key_t)hash_string(indexDir + "/edge.bin" + "edge_");
		std::ifstream::pos_type size = file8.tellg();
		int shmid8 = shmget(key, size, IPC_CREAT | IPC_EXCL);


	if(shmid8 != -1)
	{
	    std::cerr << "IF 8 \n";
	    CLI::AutoTimer timer{"Loading edges", CLI::Timer::Big};
    	std::string pfile = indexDir + "/edge.bin";
    	sdsl::load_from_file(edge_, pfile);
		shmctl(shmid8, IPC_RMID, NULL);
		shmid8 = shmget(key, size, IPC_CREAT | 0666);
		char* ptr_ch = (char*)shmat(shmid8,NULL,0);
		char * memblock =  new char[size];
		file8.seekg (0, std::ios::beg);
		file8.read (memblock, size);
		
		if(ptr_ch == NULL)
		{
			std::cerr << "POINTER IS NULL\n";		
		}
		else
		{
			std::cerr << "POINTER IS NOT NULL\n";	
		}        
		for(int i = 0; i < file8.tellg(); i++)
		{
			*ptr_ch++ = memblock[i];
		}
		
	}



	else
	{
		std::cerr << "ELSE 8 \n";
		shmid8 = shmget(key, size, IPC_CREAT | 0666);
		char* ptr_ch = (char*)shmat(shmid8,NULL,0);
		if(ptr_ch == NULL)
		{
			std::cerr << "POINTER IS NULL\n";		
		}
		else
		{
			std::cerr << "POINTER IS NOT NULL\n";	
		}      
		memstream edgestream(ptr_ch, size);
		edge_.load(edgestream);   	
	}
  }
  /*
  {
    CLI::AutoTimer timer{"Loading edges", CLI::Timer::Big};
    std::string pfile = indexDir + "/revedge.bin";
    sdsl::load_from_file(revedge_, pfile);
  }
  */
    std::cerr << "SLEEEEPPP : \n";
    //sleep(10);
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
      //std::cerr << "SEGFAULT 2 start : \n";
      auto& pvec = contigTable_[rank];
      //std::cerr << "SEGFAULT 2 end : \n";
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
      //std::cerr << "BK1 : \n";
      // relative offset of this k-mer in the contig
      uint32_t relPos = static_cast<uint32_t>(pos - sp);

      // start position of the next contig - start position of this one
      auto clen = static_cast<uint64_t>(contigEnd + 1 - sp);
      // auto clen =
      // cPosInfo_[rank].length();//static_cast<uint64_t>(contigSelect_(rank +
      // 1) + 1 - sp);

      // how the k-mer hits the contig (true if k-mer in fwd orientation, false
      // otherwise)
      //std::cerr << "BK2 : \n";
      bool hitFW = (keq == KmerMatchType::IDENTITY_MATCH);
      //std::cerr << "BK3 : \n";
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
      std::cerr << "SEGFAULT 3 start : \n";
      auto& pvec = contigTable_[rank];
      std::cerr << "SEGFAULT 3 end : \n";
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
  std::cerr << "SEGFAULT 1 start : \n";
  return contigTable_[contigRank];
  std::cerr << "SEGFAULT 1 end : \n";
}

const std::string& PufferfishIndex::refName(uint64_t refRank) {
  return refNames_[refRank];
}

uint32_t PufferfishIndex::refLength(uint64_t refRank) const {
  return refLengths_[refRank];
}

const std::vector<std::string>& PufferfishIndex::getRefNames() {
  return refNames_;
}

const std::vector<uint32_t>& PufferfishIndex::getRefLengths() const {
  return refLengths_;
}

