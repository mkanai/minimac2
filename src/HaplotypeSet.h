#ifndef __HAPLOTYPESET_H__
#define __HAPLOTYPESET_H__

#include "IntArray.h"
#include "StringArray.h"
#include "InputFile.h"
#include "VcfFileReader.h"

#include <map>

using namespace std;


class HaplotypeSet
   {
   public:
      int         count;
      int         markerCount;
      StringArray labels;
      char **     haplotypes;
      char *      major;
      float **    freq;
      bool        translate;
      bool        flipped;

      HaplotypeSet();
	  ~HaplotypeSet();

	  static const int VCF_HEADING_FIELDS = 9;
      static bool hapmapFormat;
      static bool vcfReference;
      static bool TransposeReference;
      static bool ImputeReference;

      static int  vcfFields;

      static double startposition;
      static double endposition;
      static double chr;
      static String schr;
      static String mimicArray;


      static double corestartposition;
      static double coreendposition;
      static double window;

	  static int LoadReferenceSNPsFromVcf (const char * filename, StringArray & markerList, bool rs);
	  static int LoadReferenceSNPsFromVcf (String& vcf, StringArray & markerList, bool rs, String exclude);
	  static int Base2Int(char allele, int pos);

      void LoadImputeHaplotypes(IFILE file);

      static bool setupFile(String& inputVcf, VcfFileReader& reader, VcfHeader& header);
      void processKeepSample(VcfRecord& record, int sampleIndex);
      void processExcludeSample(VcfRecord& record, int sampleIndex);


      //load both positions (chr$chr:$pos for all, ignoring rsIDs) and haplotypes
      void LoadVcf(IFILE file);
      void LoadVcf(String file, String exclude);

      void LoadLegend(const char * filename);
      void LoadSamples(const char * filename);
      void LoadShapeITHaplotypes(const char * filename, StringArray & markerLabels, bool rs = false, String chr="");


      void LoadHaplotypes(const char * filename, bool allowMissing = false, bool refHap = false, String exclude = "");
      void LoadHaplotypes(IFILE file, bool allowMissing = false, bool refHap = false);
      // the function below is for loading target haplotypes when reference haplotypes are in vcf (0/1 allele) format
      void LoadTargetHaplotypes(const char * filename, HaplotypeSet & haps, IntArray & markerIndex, StringArray & markerLabels, bool allowMissing = false);

	  void ClipHaplotypes(int & firstMarker, int & lastMarker);

	  void ListMajorAlleles();
	  void ListMajorAlleleStrings();
	  void Transpose();

	  void CalculateFrequencies();
	  void CompareFrequencies(HaplotypeSet & sets, IntArray & index, StringArray & names,int startIndex, int stopIndex);

	  const char * MajorAlleleLabel(int marker);
	  const char * MinorAlleleLabel(int marker);

	  StringArray AlleleOneLabelString;
	  StringArray AlleleTwoLabelString;

   private:
	  static const char * bases[8];
   };

#endif

