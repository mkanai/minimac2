#include "HaplotypeSet.h"
#include "StringArray.h"
#include "MemoryAllocators.h"
#include "MemoryInfo.h"
#include "Error.h"
#include <string>

#define square(x) ((x) * (x))

const char * HaplotypeSet::bases[8] = {"", "A", "C", "G", "T", "D", "I", "R"};
bool HaplotypeSet::hapmapFormat = false;
bool HaplotypeSet::vcfReference = false;
bool HaplotypeSet::TransposeReference = false;
bool HaplotypeSet::ImputeReference = false;
double HaplotypeSet::startposition = 0.0;
double HaplotypeSet::endposition = 300000000.0; //300Mb
double HaplotypeSet::corestartposition = 0.0;
double HaplotypeSet::coreendposition = 300000000.0; //300Mb
String HaplotypeSet::schr ="20";
String HaplotypeSet::mimicArray ="";

double HaplotypeSet::chr = 0.0;
double HaplotypeSet::window = 0.0;
int HaplotypeSet::vcfFields = 0;



HaplotypeSet::HaplotypeSet()
   {
   haplotypes = NULL;
   freq = NULL;
   translate = true;
   markerCount = 0;
   count = 0;
   flipped = false;
   }

HaplotypeSet::~HaplotypeSet()
   {
   if (haplotypes != NULL)
      FreeCharMatrix(haplotypes, flipped ? markerCount : count);


   if (freq != NULL)
      FreeFloatMatrix(freq, 5);

   if (major != NULL)
      delete [] major;
   }


bool HaplotypeSet::setupFile(String& inputVcf, VcfFileReader& reader, VcfHeader& header)
{

    bool overlap = false;
    bool params = false;

    // Open the file.
    reader.open(inputVcf, header);

    // If a region was specified, read the index.
    if(HaplotypeSet::window > 0)
    {
            try
            {
                reader.readVcfIndex();
                printf("Using VCF index\n");
            }
            catch(std::exception& e)
            {
                std::cerr
                    << "Warning index file not found, not using an index\n";
            }

            printf("%s\n",HaplotypeSet::schr.c_str());

            reader.set1BasedReadSection(HaplotypeSet::schr.c_str(), HaplotypeSet::startposition, HaplotypeSet::endposition, overlap);

    }
    return(true);
}


void HaplotypeSet::processKeepSample(VcfRecord& record, int sampleIndex)
{
    const char* alleles;
    for(int i = 0; i < record.getNumGTs(sampleIndex); i++)
    {
        alleles = record.getAlleles(record.getGT(sampleIndex, i));
        // Alleles are the alleles for this GT.
        // TODO
    }
}


void HaplotypeSet::processExcludeSample(VcfRecord& record, int sampleIndex)
{
    // Check to see if we should keep this position.
    // TODO

    const char* alleles;
    for(int i = 0; i < record.getNumGTs(sampleIndex); i++)
    {
        alleles = record.getAlleles(record.getGT(sampleIndex, i));
        // Alleles are the alleles for this GT.
        // TODO
    }
}


int HaplotypeSet::LoadReferenceSNPsFromVcf (const char * filename, StringArray & markerList, bool rs)
{
  int count = 0;
  IFILE file = ifopen(filename, "rb");
  if (file == NULL)
    error("File [%s] with phased haplotypes could not be opened\n", filename);

  String      buffer;
  StringArray tokens;

  // load markers, and count missingness or unphased
  // load positions ($chr:$pos for all markers ignoring rsIDs, added :$alt for non-SNPs because some have the same positions:e.g., 14:50032210; 55811665)
  //printf("Loading markers from VCF-format input file\n");
  HaplotypeSet::vcfFields = 0;
  int numDiploid = 0;
  int linenum = 0;
  while (!ifeof(file))
  {
    linenum ++;
    buffer.ReadLine(file);
    tokens.ReplaceTokens(buffer);
    if (tokens.Length() == 0) continue;

    //skip comment lines
    if (tokens[0].SubStr(0,2) == "##") continue;    // count number of individuals from the header line

    if (tokens[0].SubStr(0,6) == "#CHROM")
    {
		// HaplotypeSet::vcfFields = tokens.Length();
		// numDiploid = HaplotypeSet::vcfFields - HaplotypeSet::VCF_HEADING_FIELDS;
		// count = numDiploid * 2;
      continue;
    }


    if (tokens[1].AsDouble() < HaplotypeSet::startposition) continue;
    if (tokens[1].AsDouble() > HaplotypeSet::endposition) break;

    StringArray altAlleles;
    altAlleles.ReplaceTokens(tokens[4], ",");
    if (altAlleles.Length() > 1) continue; // skip markers with >1 ALT alleles

    String MarkerPosition;

    if (count==0)
    {
		   HaplotypeSet::vcfFields = tokens.Length();
		   int diplo = 0, haplo = 0;
    	   for (int i = VCF_HEADING_FIELDS; i < vcfFields; i++)
    	    {
    	      StringArray ThisPersonInformation;
    	      ThisPersonInformation.ReplaceTokens(tokens[i], ":");
    	      int genotypeLength = ThisPersonInformation[0].Length();
    	      if (ThisPersonInformation[0].Length() == 3)
                  diplo++;
        	  if (ThisPersonInformation[0].Length() == 1)
                  haplo++;
        	//  printf(" %i Haplotypes identified");
    	    }
    	   count = diplo * 2 + haplo;
    }

    if (tokens.Length() != HaplotypeSet::vcfFields)
      error("Each line in vcf input file should contain exact %d fields (%d header fields + %d individuals)\n"
            "while line %d does not:\n\n", HaplotypeSet::vcfFields, HaplotypeSet::VCF_HEADING_FIELDS, numDiploid, linenum);


    // rs numbers?
    if (rs && tokens[2] != ".")
      {
        	 MarkerPosition = tokens[2].Trim();
      }
    else
      {
		MarkerPosition = tokens[0]+":"+tokens[1];
		if (tokens[3].Length() > 1 || tokens[4].Length() > 1)
		{
		  String extraIndelAllele = tokens[3]+"_"+tokens[4];
		  int extraLength = extraIndelAllele.Length();
		  if (extraLength > 5) extraLength = 5;
		  MarkerPosition = MarkerPosition+":"+extraIndelAllele.SubStr(0,extraLength);
		}

      }
    //printf("%s\n",MarkerPosition.c_str());
    markerList.Add(MarkerPosition);
  }
  ifclose(file);

  return count;
}

int HaplotypeSet::LoadReferenceSNPsFromVcf (String& vcf, StringArray & markerList, bool rs, String exclude)
{
  int count = 0;

  VcfFileReader reader;
  VcfHeader header;

  // The first pass through the file is site only.
  reader.setSiteOnly(true);

  if(!setupFile(vcf, reader, header))
  {
      std::cerr << "Failed to open/setup the input vcf file for reading\n";
      std::cerr << "Exiting.\n";
      return(-1);
  }

  count = header.getNumSamples();

  // Read the variant ids into a StringArray.
  VcfRecord record;
  StringArray chromPosArray,id;
  String chromPosString,ref,alt;

  while(reader.readRecord(record))
  {
      //ONLY PASS?
	  if (record.get1BasedPosition() < HaplotypeSet::startposition) continue;
	  if (record.get1BasedPosition() > HaplotypeSet::endposition) break;

	  if(!chr) chr=atof(record.getChromStr());

      id.ReplaceTokens(record.getIDStr(), ":");
      if(!rs || (id[0] == record.getChromStr())) {
          chromPosString = record.getChromStr();
          chromPosString += ":";
          chromPosString += record.get1BasedPosition();
          ref=record.getRefStr();
          alt=record.getAltStr();

          //Indel
          if (ref.Length() > 1 || alt.Length() > 1)
          {
              String extraIndelAllele = ref+"_"+alt;
              int extraLength = extraIndelAllele.Length();
              if (extraLength > 5) extraLength = 5;
              chromPosString = chromPosString+":"+extraIndelAllele.SubStr(0,extraLength);
    	  }
      } else {
          chromPosString = id[0];
      }

      markerList.Add(chromPosString);

      if (record.get1BasedPosition() >= (HaplotypeSet::startposition+HaplotypeSet::window) && HaplotypeSet::startposition > 0 && HaplotypeSet::corestartposition == 0.0)
    	  HaplotypeSet::corestartposition = record.get1BasedPosition();
      if (record.get1BasedPosition() <= (HaplotypeSet::endposition-HaplotypeSet::window))
    	  HaplotypeSet::coreendposition = record.get1BasedPosition();
  }

  // Close the file and reopen in to process the haplotypes.
  reader.close();
  return count;
}

void HaplotypeSet::LoadLegend(const char * filename)
{
  IFILE file = ifopen(filename, "rb");
  if (file == NULL)
    error("File [%s] with phased haplotypes could not be opened\n", filename);

  AlleleOneLabelString.Clear();
  AlleleTwoLabelString.Clear();

  String      buffer;
  StringArray tokens;

  printf("Loading alleles from HapMap-style legend file ...\n");
  //header
  buffer.ReadLine(file);
  while (!ifeof(file))
  {
    buffer.ReadLine(file);
    tokens.ReplaceTokens(buffer);

    if (tokens.Length() == 0) continue;

    if (tokens.Length() != 4)
      error("Each line should list the marker name, position and alleles,\n"
            "but the following line includes %d items (instead of 4):\n\n"
            "%s\n", tokens.Length(), (const char *) buffer);

    if (tokens[1].AsDouble() < startposition) continue;
    if (tokens[1].AsDouble() > endposition) break;

    AlleleOneLabelString.Add(tokens[2]);
    AlleleTwoLabelString.Add(tokens[3]);
  }
  ifclose(file);
}


void HaplotypeSet::LoadSamples(const char * filename)
{
  IFILE file = ifopen(filename, "rb");
  if (file == NULL)
    error("File [%s] with sample information could not be opened\n", filename);

  String      buffer;
  StringArray tokens;

  printf("Loading ShapeIT sample information ...\n");

  // skip two header lines
  buffer.ReadLine(file);
  buffer.ReadLine(file);
  while (!ifeof(file))
  {
    buffer.ReadLine(file);
    tokens.ReplaceTokens(buffer);

    if (tokens.Length() == 0) continue;


      
    labels.Add(tokens[0]+"->"+tokens[1]);
    labels.Add(tokens[0]+"->"+tokens[1]);

  }
  ifclose(file);

}


void HaplotypeSet::LoadHaplotypes(const char * filename, bool allowMissing, bool refHap, String exclude)
   {

	 if (refHap && vcfReference)
	 {
	    LoadVcf(filename, exclude);
	    return;
	 }

   IFILE f = ifopen(filename, "rb");

   if (f == NULL)
      {
      error("File [%s] with phased haplotypes could not be opened\n", filename);
      return;
      }

   LoadHaplotypes(f, allowMissing, refHap);
   ifclose(f);
   }


void HaplotypeSet::LoadImputeHaplotypes (IFILE file)
{
  // Don't load haplotypes unless we have a marker list
  if (markerCount == 0)
  {
    printf("  WARNING -- Since no legend file was provided, haplotype file will be ignored\n\n");
    return;
  }

  count = 0;

  printf("Loading IMPUTE-style phased haplotypes ...\n");

  String      buffer;
  StringArray tokens;

   // first load IDfile if IDs2Remove is non-empty
   // ylwtx: not copied
   int numberIDs = -9;

  // In the first pass, we simply count if the number of markers (lines)
  // to see if agree with markerCount
  int numHapLines = 0;
  while (!ifeof(file))
  {
    buffer.ReadLine(file);
    if (buffer.Length() == 0) continue;
    numHapLines++;
  }
  if (numHapLines != markerCount)
    error("According to the legend file, there should be %d markers but \n"
          "the IMPUTE-style haplotype file contains %d lines\n\n", markerCount, numHapLines);

  // In the second pass, we simply count the number of haplotypes, by counting
  // the number of characters in the first line, which is for the first marker
  ifrewind(file);
  buffer.ReadLine(file);
  tokens.ReplaceTokens(buffer);
  count = tokens.Length();
  int numberExpIDsInHap = count;
  if (count != numberIDs && numberIDs != -9)
	  error("Number of IDs in ID file %d does not agree with what is in the reference haplotype %d\n\n",
		  numberIDs, count);
  //if (IDs2Remove.Entries() > 0) count = IDs2Keep.Length();

  // Check if we got some valid input
  if (count == 0 || markerCount == 0) return;

  // Then, we allocate memory for storing the phased haplotypes
  haplotypes = AllocateCharMatrix(count, markerCount);
  major = new char [markerCount];
  if (haplotypes == NULL) error("Failed to allocate memory for reference haplotypes\n\n");

  // And finally, we load the data in a second pass
  ifrewind(file);

  int line = 0, markerindex = 0;
  while (!ifeof(file))
  {
    line++;
    buffer.ReadLine(file);
    buffer.Trim();
    if (buffer.Length() == 0) continue;

    //if (buffer.Length() != 2 * markerCount - 1) continue;
    if (buffer.Length() != 2 * numberExpIDsInHap - 1)
      error("Expect %d alleles for each marker/line for IMPUTE-style haplotype file, while line %d is invalid.\n\n",
        		numberExpIDsInHap, line);

    bool badchar = false;
    // ylwtx: if (IDs2Keep.Length() > 0)
    if (numberExpIDsInHap == count)
    {
        for (int i = 0; i < count; i++)
        {
          if (buffer[i * 2] == '0')
             haplotypes[i][markerindex] = 1;
          else if (buffer[i * 2] == '1')
             haplotypes[i][markerindex] = 2;
          else
             badchar = true;

          if (badchar || (buffer[i * 2 + 1] != ' ' && buffer[i * 2 + 1] != 0))
             error("IMPUTE-style haplotype file should include a series of '0's and '1's,\n"
                   "separated by spaces. However, an unexpected character was\n"
                   "encountered in line %d.\n", line);
        }
    } // if (numberExpIDsInHap == count)
    markerindex++;
  }
}


void HaplotypeSet::LoadShapeITHaplotypes (const char * filename, StringArray & markerList, bool rs, String chr)
{

  IFILE file = ifopen(filename, "rb");
  if (file == NULL)
	 error("File [%s] with sample information could not be opened\n", filename);

  count = 0;

  String      buffer;
  StringArray tokens;

  // In the first pass, we simply count if the number of markers (lines)
  // to see if agree with markerCount
  markerCount = 0;
  while (!ifeof(file))
  {
    buffer.ReadLine(file);
    if (buffer.Length() == 0) continue;
    markerCount++;
  }

  ifrewind(file);
  buffer.ReadLine(file);
  tokens.ReplaceTokens(buffer);

  int column_offset = 5;

  count = tokens.Length() - column_offset;
  int numberExpIDsInHap = labels.Length();

  if (count != numberExpIDsInHap )
	  error("Number of IDs in sample file (%d) does not agree with what is in the haplotype file (%d)\n\n",
			  numberExpIDsInHap, count);

  // Check if we got some valid input
  if (count == 0 || markerCount == 0) return;

  // Then, we allocate memory for storing the phased haplotypes
  haplotypes = AllocateCharMatrix(count, markerCount);
  major = new char [markerCount];
  if (haplotypes == NULL) error("Failed to allocate memory for framework haplotypes\n\n");

  // And finally, we load the data in a second pass
  ifrewind(file);

  int line = 0, markerindex = 0;
  while (!ifeof(file))
  {
    line++;
    buffer.ReadLine(file);
    buffer.Trim();
    tokens.ReplaceTokens(buffer);

    if (tokens.Length() == 0) continue;

    if (tokens.Length() - column_offset != numberExpIDsInHap)
      error("Expect %d alleles for each marker/line for IMPUTE-style haplotype file, while line %d is invalid.\n\n",
        		numberExpIDsInHap, line);

    bool badchar = false;

    if(rs)
       markerList.Add(tokens[1].Trim());
    else
       markerList.Add(chr+":"+tokens[2].Trim());

    int al0, al1;
    if (tokens[3] == "A" || tokens[3] == "a") al0 = 1;
    else if (tokens[3] == "C" || tokens[3] == "c") al0 = 2;
    else if (tokens[3] == "G" || tokens[3] == "g") al0 = 3;
    else if (tokens[3] == "T" || tokens[3] == "t") al0 = 4;
    else error("SNP alleles can only contain alleles A (A or a),\n"
                "C (C or c), G (G or g), and T (T or t).\n");

    if (tokens[4] == "A" || tokens[4] == "a") al1 = 1;
    else if (tokens[4] == "C" || tokens[4] == "c") al1 = 2;
    else if (tokens[4] == "G" || tokens[4] == "g") al1 = 3;
    else if (tokens[4] == "T" || tokens[4] == "t") al1 = 4;
    else error("SNP alleles can only contain alleles A (A or a),\n"
                "C (C or c), G (G or g), and T (T or t).\n");

	for (int i = 0; i < count; i++)
	{
	  if (tokens[i+ column_offset] == "0")
		 haplotypes[i][markerindex] = al0;
	  else if (tokens[i+ column_offset] == "1")
		 haplotypes[i][markerindex] = al1;
	  else
		 error("IMPUTE-style haplotype file should include a series of '0's and '1's,\n"
			   "separated by spaces. However, an unexpected character was\n"
			   "encountered in line %d.\n", line);
	}
    markerindex++;
  }
}


void HaplotypeSet::LoadTargetHaplotypes(const char * filename, HaplotypeSet & haps, IntArray & markerIndex, StringArray & markerLabels, bool allowMissing)
{
  IFILE file = ifopen(filename, "rb");
  if (file == NULL)
    error("File [%s] with phased haplotypes could not be opened\n", filename);

  if (markerCount == 0)
  {
    printf("  WARNING -- Since no marker list was provided, haplotype file will be ignored\n\n");
    return;
  }

  String      buffer;
  StringArray tokens;

  // In the first pass, we simply count the number of non-blank lines
  while (!ifeof(file))
  {
    buffer.ReadLine(file);
    tokens.ReplaceTokens(buffer);
    if (!tokens.Length()) continue;
    count++;
  }

  // Check if we got some valid input
  if (count == 0 || markerCount == 0) return;

  // Then, we allocate memory for storing the phased haplotypes
  haplotypes = AllocateCharMatrix(count, markerCount);
  major = new char [markerCount];
  if (haplotypes == NULL) error("Failed to allocate memory for target haplotypes\n\n");

  labels.Dimension(count);

  // And finally, we load the data in a second pass
  ifrewind(file);

  int line = 0, index = 0;
  while (!ifeof(file))
  {
    line++;
    buffer.ReadLine(file);
    tokens.ReplaceTokens(buffer);

    if (tokens.Length() == 0) continue;

    labels[index] = tokens[0];

    int hapstart = tokens.Length() - 1;
    int offset = markerCount;

    while ((offset -= tokens[hapstart].Length()) > 0 && hapstart > 0)
      hapstart--;

    if (offset != 0)
      error("The haplotype file format was not recognized\n"
            "(Problem occured reading haplotype #%d in line #%d)\n\n"
            "Check that the number of markers matches the SNPs list\n",
            ++line, ++index);

    for (int i = 0; i < markerCount; i++)
    {
      if (offset == tokens[hapstart].Length())
      {
        offset = 0;
        hapstart++;
      }

      char allele = tokens[hapstart][offset++];
      char al;

      if (translate)
        switch (allele)
        {
          case '1' : allele = 'A'; break;
          case '2' : allele = 'C'; break;
          case '3' : allele = 'G'; break;
          case '4' : allele = 'T'; break;
          case '5' : allele = 'D'; break;
          case '6' : allele = 'I'; break;
          case '7' : allele = 'R'; break;
        }

      switch (allele)
      {
        case 'A' : case 'a' : al = 1; break;
        case 'C' : case 'c' : al = 2; break;
        case 'G' : case 'g' : al = 3; break;
        case 'T' : case 't' : al = 4; break;
        case 'D' : case 'd' : al = 5; break;
        case 'I' : case 'i' : al = 6; break;
        case 'R' : case 'r' : al = 7; break;
        case '0' : case '.' : case 'N' : case 'n' :
          if (allowMissing) { al = 0; break; }
        default  :
          error("Haplotypes can only contain alleles A ('A', 'a' or '1'),\n"
                "C ('C', 'c' or 2), G ('G', 'g' or 3), T ('T', 't' or '4'),\n"
                "D ('D', 'd' or 5), I ('I', 'i' or 6) and R ('R', 'r' or 7).\n");
      }

      String alString;
      switch (al)
      {
        case 1 : alString = "A"; break;
        case 2 : alString = "C"; break;
        case 3 : alString = "G"; break;
        case 4 : alString = "T"; break;
        case 5 : alString = "D"; break;
        case 6 : alString = "I"; break;
        case 7 : alString = "R"; break;
        case 0 : alString = "0"; break;
        default  :
          error("Haplotypes can only contain alleles A ('A', 'a' or '1'),\n"
                "C ('C', 'c' or 2), G ('G', 'g' or 3), T ('T', 't' or '4'),\n"
                "D ('D', 'd' or 5), I ('I', 'i' or 6) and R ('R', 'r' or 7).\n");
      }

      // markerIndex[i] >= 0 <=> target SNP in reference
      if (markerIndex[i] >= 0 && alString != "0")
      {
        if (haps.AlleleOneLabelString[markerIndex[i]] == alString) al = 1;
        else if (haps.AlleleTwoLabelString[markerIndex[i]] == alString) al = 2;
        else error("For marker '%s': target allele [%s] matches neither in reference: %s or %s\n\n", (const char *) markerLabels[i], (const char * )alString, (const char *) haps.AlleleOneLabelString[markerIndex[i]], (const char *) haps.AlleleTwoLabelString[markerIndex[i]]);
      }

      haplotypes[index][i] = al;
    }
    index++;
  }

  ifclose(file);
}

void HaplotypeSet::LoadHaplotypes(IFILE file, bool allowMissing, bool refHap)
{

  if (refHap && ImputeReference)
  {
    LoadImputeHaplotypes(file);
    return;
  }

   // Don't load haplotypes unless we have a marker list
   if (markerCount == 0)
      {
      printf("  WARNING -- Since no marker list was provided, haplotype file will be ignored\n\n");
      return;
      }

   String      buffer;
   StringArray tokens;

   // In the first pass, we simply count the number of non-blank lines
   while (!ifeof(file))
      {
      buffer.ReadLine(file);
      tokens.ReplaceTokens(buffer);

      if (!tokens.Length()) continue;

      count++;
      }

   // Check if we got some valid input
   if (count == 0 || markerCount == 0)
      return;

   // Then, we allocate memory for storing the phased haplotypes
   haplotypes = AllocateCharMatrix(count, markerCount);
   major = new char [markerCount];
   if (haplotypes == NULL) error("Failed to allocate memory for reference haplotypes\n\n");

   labels.Dimension(count);

   // And finally, we load the data in a second pass
   ifrewind(file);

   int line = 0, index = 0;
   while (!ifeof(file))
      {
      line++;
      buffer.ReadLine(file);
      tokens.ReplaceTokens(buffer);

      if (tokens.Length() == 0) continue;

      labels[index] = tokens[0];

      int hapstart = tokens.Length() - 1;
      int offset = markerCount;

      while ((offset -= tokens[hapstart].Length()) > 0 && hapstart > 0)
         hapstart--;

      if (offset != 0)
         error("The haplotype file format was not recognized\n"
               "(Problem occured reading haplotype #%d in line #%d)\n\n"
               "Check that the number of markers matches the SNPs list\n"
                "Offset:%d Token.length:%d - %d, MarkerCount:%d, Hapstart:%d",
               ++index, line, offset,tokens.Length(),tokens[hapstart].Length(), markerCount,hapstart);

      for (int i = 0; i < markerCount; i++)
         {
         if (offset == tokens[hapstart].Length())
            {
            offset = 0;
            hapstart++;
            }

         char allele = tokens[hapstart][offset++];
         char al;

         if (translate)
            switch (allele) {
               case '1' : allele = 'A'; break;
               case '2' : allele = 'C'; break;
               case '3' : allele = 'G'; break;
               case '4' : allele = 'T'; break;
               case '5' : allele = 'D'; break;
               case '6' : allele = 'I'; break;
               case '7' : allele = 'R'; break;
            }

         switch (allele)
               {
               case 'A' : case 'a' : al = 1; break;
               case 'C' : case 'c' : al = 2; break;
               case 'G' : case 'g' : al = 3; break;
               case 'T' : case 't' : al = 4; break;
               case 'D' : case 'd' : al = 5; break;
               case 'I' : case 'i' : al = 6; break;
               case 'R' : case 'r' : al = 7; break;
               case '0' : case '.' : case 'N' : case 'n' :
                  if (allowMissing) { al = 0; break; }
               default  :
          error("Haplotypes can only contain alleles A ('A', 'a' or '1'),\n"
                "C ('C', 'c' or 2), G ('G', 'g' or 3), T ('T', 't' or '4'),\n"
                "D ('D', 'd' or 5), I ('I', 'i' or 6) and R ('R', 'r' or 7).\n");
               }

         haplotypes[index][i] = al;
         }
      index++;
      }
   }

void HaplotypeSet::LoadVcf(IFILE file)
{
  String      buffer;
  StringArray tokens;

  // Check if we got some valid input
  if (markerCount == 0)
  {
    printf("  WARNING -- Since we read ZERO marker from vcf input file, the vcf input will be ignored\n\n");
    return;
  }

  // Check if we got some valid input
  if (count == 0)
  {
    printf("  WARNING -- Since we read ZERO individuals from vcf input file, the vcf input will be ignored\n\n");
    return;
  }

  // We first allocate memory for storing the phased haplotypes
  haplotypes = AllocateCharMatrix(count, markerCount);
  major = new char [markerCount];
  if (haplotypes == NULL) error("Failed to allocate memory for reference haplotypes\n\n");

  // load haplotypes as well as reset allele labels
  printf("Loading haplotypes from VCF-format input file (WARNING: genotypes would be treated as phased regardless of delimiter) ...\n");

  int linenum = 0;
  int markerIndex = 0;

  while (!ifeof(file))
  {
    linenum ++;
    buffer.ReadLine(file);
    tokens.ReplaceTokens(buffer);

    if (tokens.Length() == 0) continue;

    //skip comment lines and header line
    if (tokens[0].SubStr(0,1) == "#") continue;

    // already checked #fields in the first pass

    if (tokens[1].AsDouble() < startposition) continue;
    if (tokens[1].AsDouble() > endposition) break;

    StringArray altAlleles;
    altAlleles.ReplaceTokens(tokens[4], ",");
    if (altAlleles.Length() > 1) continue; // skip markers with >1 ALT alleles


    // set chromosome
    if(!chr) chr=tokens[0].AsDouble();

    if (tokens[1].AsDouble() >= (startposition+window) && startposition > 0 && corestartposition == 0.0)
       corestartposition = tokens[1].AsDouble();
    if (tokens[1].AsDouble() <= (endposition-window))
       coreendposition = tokens[1].AsDouble();

    // does not check for multi-allelic
    // store allele label
    AlleleOneLabelString.Add(tokens[3]);
    AlleleTwoLabelString.Add(tokens[4]);

    int al1, al2;
    if (tokens[3].Length() == 1 && tokens[4].Length() == 1)
    //SNP
    {
      if (tokens[3] == "A" || tokens[3] == "a") al1 = 1;
      else if (tokens[3] == "C" || tokens[3] == "c") al1 = 2;
      else if (tokens[3] == "G" || tokens[3] == "g") al1 = 3;
      else if (tokens[3] == "T" || tokens[3] == "t") al1 = 4;
      else error("SNP alleles can only contain alleles A (A or a),\n"
                "C (C or c), G (G or g), and T (T or t).\n"
                "but chr%s:%s has invalid REF allele %s.\n\n", (const char *) tokens[0], (const char *) tokens[1], (const char *) tokens[3]);

      if (tokens[4] == "A" || tokens[4] == "a") al2 = 1;
      else if (tokens[4] == "C" || tokens[4] == "c") al2 = 2;
      else if (tokens[4] == "G" || tokens[4] == "g") al2 = 3;
      else if (tokens[4] == "T" || tokens[4] == "t") al2 = 4;
      else error("SNP alleles can only contain alleles A (A or a),\n"
                "C (C or c), G (G or g), and T (T or t).\n"
                "but chr%s:%s has invalid REF allele %s.\n\n", (const char *) tokens[0], (const char *) tokens[1], (const char *) tokens[4]);
    }
    else
    //indel
    {
      al1 = 7;
      if (tokens[3].Length() > tokens[4].Length()) al2 = 5;
      else al2 = 6;

    }

    // only 0 and 1 allowed
    bool badchar = false;
    int haplotypeIndex = 0;

    for (int i = VCF_HEADING_FIELDS; i < vcfFields; i++)
    {

      int personIndex = i - VCF_HEADING_FIELDS;

      StringArray ThisPersonInformation;
      ThisPersonInformation.ReplaceTokens(tokens[i], ":");
      int genotypeLength = ThisPersonInformation[0].Length();

      if (!(ThisPersonInformation[0].Length() == 3 || ThisPersonInformation[0].Length() == 1))
        error("Each genotype has to be a string of three characters: allele1, delimiter, and allele2. \n"
              "The following genotype in line %d does not: %s\n\n", linenum, (const char *) ThisPersonInformation[0]);

      //load first allele
      if (ThisPersonInformation[0].SubStr(0,1) == "0") haplotypes[haplotypeIndex][markerIndex] = al1;
      else if (ThisPersonInformation[0].SubStr(0,1) == "1") haplotypes[haplotypeIndex][markerIndex] = al2;
      else badchar = true;

      	  haplotypeIndex ++;

      if (badchar) error("Allele 1 in vcf can only be '0' or '1'\n"
                         "The following Allele(s) in line %d does not: %s\n\n", linenum, (const char *) ThisPersonInformation[0]);

      if (ThisPersonInformation[0].Length() == 3){
      //load second allele
		  if (ThisPersonInformation[0].SubStr(2,1) == "0") haplotypes[haplotypeIndex][markerIndex] = al1;
		  else if (ThisPersonInformation[0].SubStr(2,1) == "1") haplotypes[haplotypeIndex][markerIndex] = al2;
		  else badchar = true;

		  if (badchar) error("Allele 2 in vcf can only be '0' or '1', separated by any one-character delimiter.\n"
							 "The following Allele(s) in line %d does not: %s\n\n", linenum, (const char *) ThisPersonInformation[0]);
		  haplotypeIndex ++;
      }
    }
    markerIndex ++;
  }
}


void HaplotypeSet::LoadVcf(String file, String exclude)
{

	VcfFileReader reader;
	VcfHeader header;
	VcfRecord record;

	int numSamples;

	// Now we want the genotype information.
    reader.setSiteOnly(false);
    if(!setupFile(file, reader, header))
    {
        std::cerr << "Failed to open/setup the input vcf file for reading\n";
        std::cerr << "Exiting.\n";
    }

    // Determine the position of the sample to exclude.
    int excludeIndex = -1;
    numSamples = header.getNumSamples();

    if(exclude != "")
    {

    	excludeIndex = header.getSampleIndex(exclude);
        if(excludeIndex >= numSamples || excludeIndex == -1)
        {
            // Index is out of range of the samples.
            excludeIndex = -1;
            count = numSamples;
        }
        else
        	count = numSamples -1;
    }
    else{
        count = numSamples;
    }

    count = count * 2;
    // We first allocate memory for storing the phased haplotypes
    printf("%d - %d\n\n",count,markerCount);

  	haplotypes = AllocateCharMatrix(count, markerCount);
  	major = new char [markerCount];
  	if (haplotypes == NULL) error("Failed to allocate memory for reference haplotypes\n\n");

  	// load haplotypes as well as reset allele labels
  	printf("Loading haplotypes from VCF-format input file (WARNING: genotypes would be treated as phased regardless of delimiter) ...\n");

  	// Only store the GT field.
    VcfRecordGenotype::addStoreField("GT");

     // Loop through the file and read the haplotypes.
    const char* alleles;
    int haplotypeIndex;
    int marker = 0;

    while(reader.readRecord(record))
    {

      bool indel = false;
      int reflen = -1;
      int alttype = 7;

      //ONLY PASS?
  	  if (record.get1BasedPosition() < HaplotypeSet::startposition) continue;
  	  if (record.get1BasedPosition() > HaplotypeSet::endposition) break;

        // Loop through and get the genotype info for
        // each sample.

    	//Indels
     	if(strlen(record.getRefStr())>1 || strlen(record.getAltStr())>1){
     	   indel = true;
     	   reflen = strlen(record.getRefStr());

    	   AlleleOneLabelString.Add("R");
    	   if (strlen(record.getRefStr()) > strlen(record.getAltStr())){
    		   AlleleTwoLabelString.Add("D");
    		   alttype = 5;
    	   }
           else{
    		   AlleleTwoLabelString.Add("I");
               alttype = 6;
           }
    	}else{
        	AlleleOneLabelString.Add(record.getRefStr());
            AlleleTwoLabelString.Add(record.getAltStr());
    	}


    	haplotypeIndex = 0;

        for(int i = 0; i < excludeIndex; i++)
        {

        	alleles = record.getAlleles(record.getGT(i, 0));
        	if(indel){
        	   if(strlen(alleles)==reflen)
        		   haplotypes[haplotypeIndex][marker] = 7;
        	   else
        		   haplotypes[haplotypeIndex][marker] = alttype;
        	}
        	else
        	  haplotypes[haplotypeIndex][marker] = Base2Int(alleles[0],0);

        	haplotypeIndex ++;

   		    alleles = record.getAlleles(record.getGT(i, 1));
        	if(indel){
        	   if(strlen(alleles)==reflen)
        		   haplotypes[haplotypeIndex][marker] = 7;
        	   else
        		   haplotypes[haplotypeIndex][marker] = alttype;
        	}
        	else
   		    haplotypes[haplotypeIndex][marker] = Base2Int(alleles[0],1);
   		    haplotypeIndex ++;

        }
        // process the exclude index if not -1
        if(excludeIndex != -1)
        {
            //processExcludeSample(record, excludeIndex);

        }
        // process the rest of the samples.
        for(int i = (excludeIndex+1); i < numSamples; i++)
        {

        	alleles = record.getAlleles(record.getGT(i, 0));
        	if(indel){
        	   if(strlen(alleles)==reflen)
        		   haplotypes[haplotypeIndex][marker] = 7;
        	   else
        		   haplotypes[haplotypeIndex][marker] = alttype;
        	}
        	else
        	  haplotypes[haplotypeIndex][marker] = Base2Int(alleles[0],0);

        	haplotypeIndex ++;

   		    alleles = record.getAlleles(record.getGT(i, 1));
        	if(indel){
        	   if(strlen(alleles)==reflen)
        		   haplotypes[haplotypeIndex][marker] = 7;
        	   else
        		   haplotypes[haplotypeIndex][marker] = alttype;
        	}
        	else
   		     haplotypes[haplotypeIndex][marker] = Base2Int(alleles[0],1);

        	haplotypeIndex ++;


        }
        marker++;
       // printf("%5d\b\b\b\b\b",marker); fflush(stdout);
    }
}


void HaplotypeSet::ClipHaplotypes(int & firstMarker, int & lastMarker)

   {

   if (firstMarker < 0)

      firstMarker = 0;


   if (lastMarker < 0 || lastMarker >= markerCount - 1)

      lastMarker = markerCount - 1;


   if (firstMarker > lastMarker)

      firstMarker = lastMarker;


   int newMarkerCount = lastMarker - firstMarker + 1;


   char ** newHaplotypes = AllocateCharMatrix(count, newMarkerCount);


   for (int i = 0; i < count; i++)

      for (int j = 0; j < newMarkerCount; j++)

         newHaplotypes[i][j] = haplotypes[i][j + firstMarker];


   FreeCharMatrix(haplotypes, count);


   haplotypes = newHaplotypes;

   markerCount = newMarkerCount;

   }

void HaplotypeSet::ListMajorAlleleStrings()
// could use ListMajorAlleles but HapMap for sure will be bi-allelic, no need to calculcate for FOUR nucleotides
{
  // HapMap or vcf format: only two alleles, 0 or 1 (loaded as 1 or 2)
  for (int i = 0; i < markerCount; i++)
  {
    int freqs[3];
    for (int j = 0; j < 3; j++) freqs[j] = 0;

    for (int j = 0; j < count; j++) freqs[haplotypes[j][i]]++;

    major[i] = 1;
    if (freqs[2] >= freqs[1]) major[i] = 2;
  }
}

void HaplotypeSet::ListMajorAlleles()
   {
   for (int i = 0; i < markerCount; i++)
      {
      int freqs[8];

      for (int j = 0; j < 8; j++)
         freqs[j] = 0;

      for (int j = 0; j < count; j++)
         freqs[haplotypes[j][i]]++;

      major[i] = 1;

      for (int j = 2; j < 8; j++)
         if (freqs[j] >= freqs[major[i]])
            major[i] = j;
      }
   }

void HaplotypeSet::CalculateFrequencies()
   {
   if (freq == NULL)
      freq = AllocateFloatMatrix(8, markerCount);
//   freq--;

   for (int i = 1; i < 8; i++)
      for (int j = 0; j < markerCount; j++)
         freq[i][j] = 0;

   for (int i = 0; i < count; i++)
      for (int j = 0; j < markerCount; j++)
         if (haplotypes[i][j] != 0)
            freq[haplotypes[i][j]][j]++;

   for (int j = 0; j < markerCount; j++)
      {
      double sum = freq[1][j] + freq[2][j] + freq[3][j] + freq[4][j] + freq[5][j] + freq[6][j] + freq[7][j];

      if (sum == 0.0) continue;

      double scale = 1.0 / sum;

      for (int i = 1; i <= 7; i++)
         freq[i][j] *= scale;
      }
   }

void HaplotypeSet::CompareFrequencies(HaplotypeSet & haps, IntArray & index, StringArray & labels, int startIndex, int stopIndex)
   {
   int problems = 0;

   for (int i = 0; i < markerCount; i++)
//     for (int i = startIndex; i <= stopIndex; i++)

	   if (index[i] >= 0)

	   //if (index[i] >= startIndex && index[i] <= stopIndex)
         {
         int knownCount = 0;
         for (int j = 0; j < count; j++)
            if (haplotypes[j][i] != 0)
               knownCount++;

         int knownCountHaps = 0;
         for (int j = 0; j < haps.count; j++)
            if (haps.haplotypes[j][index[i]] != 0)
               knownCountHaps++;

         double chisq = 0.0;
         for (int j = 1; j <= 7; j++)
            if (freq[j][i] + haps.freq[j][index[i]] > 0)
               {
               double total = freq[j][i] * knownCount + haps.freq[j][index[i]] * knownCountHaps;
               double expected = total / (knownCount + knownCountHaps) * knownCount;

               double delta = freq[j][i] * knownCount - expected;

               chisq += square(delta) / expected +
                        square(delta) / (total - expected);
               }

         char a[8] = {0, 7, 6, 5, 4, 3, 2, 1};

         double chisq_after_strand_flip = 0.0;
         for (int j = 1; j <= 7; j++)
            if (freq[a[j]][i] + haps.freq[j][index[i]] > 0)
               {
               double total = freq[a[j]][i] * knownCount + haps.freq[j][index[i]] * knownCountHaps;
               double expected = total / (knownCount + knownCountHaps) * knownCount;

               double delta = freq[a[j]][i] * knownCount - expected;

               chisq_after_strand_flip += square(delta) / expected +
                                          square(delta) / (total - expected);
               }

         if (chisq > 300)
            {
            String alleles, freq1, freq2;

            for (int j = 1; j <= 7; j++)
               if (freq[j][i] + haps.freq[j][index[i]] > 0)
                  {
                  if (alleles.Length()) alleles += ",";
                  alleles += bases[j];

                  if (freq1.Length()) freq1 += ",";
                  freq1.catprintf("%.2f", freq[j][i]);

                  if (freq2.Length()) freq2 += ",";
                  freq2.catprintf("%.2f", haps.freq[j][index[i]]);
                  }

            if (hapmapFormat || vcfReference || TransposeReference || ImputeReference)
              printf("  %s for '%s': "
                     "f[%s] = [%s] vs [%s], chisq %.1f\n",
                     chisq_after_strand_flip < chisq ? "Possible strand flip" : "Mismatched frequencies",
                     (const char *) labels[i],
#					 //(const char *) haps.AlleleOneLabelString[index[i]], (const char *) haps.AlleleTwoLabelString[index[i]],
                     (const char *) alleles,
                     (const char *) freq1, (const char *) freq2,
                     chisq);
            else
              printf("  %s for '%s': "
                     "f[%s] = [%s] vs [%s], chisq %.1f\n",
                     chisq_after_strand_flip < chisq ? "Possible strand flip" : "Mismatched frequencies",
                     (const char *) labels[i],
                     (const char *) alleles, (const char *) freq1, (const char *) freq2,
                     chisq);

            problems++;
            }
         }

   if (problems)
      printf("  %d markers with potential frequency mismatches\n", problems);
   }

const char * HaplotypeSet::MajorAlleleLabel(int marker)
   {
   int hi = 1;
   for (int i = 2; i <= 7; i++)
      if (freq[i][marker] >= freq[hi][marker])
         hi = i;

   return bases[hi];
   }

const char * HaplotypeSet::MinorAlleleLabel(int marker)
   {
    int lo = 1;

	   while (freq[lo][marker] == 0 && lo < 7)
                 lo++;

		  for (int i = lo + 1; i <= 7; i++)
                    if (freq[i][marker] > 0 && freq[i][marker] < freq[lo][marker])
                                 lo = i;

             return bases[lo];
   }


void HaplotypeSet::Transpose()
{
	printf("Transposing haplotype matrix...\n");
	char **x = AllocateCharMatrix(markerCount, count);
	for (int i = 0; i < markerCount; i++)
	{
		for (int j = 0; j < count-1; j += 2)
		{
            x[i][(j>>1)] = haplotypes[j][i];
            x[i][(j>>1)+count/2] = haplotypes[j+1][i];
		}
		if (count & 1)
			x[i][count-1] = haplotypes[count-1][i];
	}
	FreeCharMatrix(haplotypes, count);
	haplotypes = x;
	flipped = !flipped;
	printf("Done.\n");
}




int HaplotypeSet::Base2Int(char allele, int pos)
{
	switch(allele)
	{
		case 'A':
			return(1);
			break;
		case 'C':
			return(2);
			break;
		case 'G':
			return(3);
			break;
		case 'T':
			return(4);
			break;
		case 'D':
			return(5);
			break;
		case 'I':
			return(6);
			break;
		case 'R':
			return(7);
			break;
		default:
			if(pos==0){
			  return(1);
			  break;
			}
			if(pos==1){
			  return(2);
			  break;
			}else
			  std::cerr << "VcfRecord::getIntAllele, unknown allele, "
					  << allele << std::endl;
	}
}
