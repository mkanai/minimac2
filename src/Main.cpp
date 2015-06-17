#include "Error.h"
#include "IntArray.h"
#include "Parameters.h"
#include "StringArray.h"
#include "StringHash.h"
#include "StringAlias.h"
#include "HaplotypeSet.h"
#include "ImputationStatistics.h"
#include "HaplotypeClipper.h"
#include "MarkovModel.h"

#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define square(x)    ((x)*(x))

time_t log_now(char *what, time_t start)
   {
   time_t stop = time(NULL);
   int seconds = stop - start;
   printf("\n%s in %d:%02d:%02d on %s\n\n", what,
          seconds / 3600, (seconds % 3600) / 60, seconds % 60,
          ctime(&stop));
   return(stop);
   }

int main(int argc, char ** argv)
   {
   setbuf(stdout, NULL);
   setbuf(stderr, NULL);

   time_t start = time(NULL);

   printf("MiniMac2 - Imputation into phased haplotypes\n"
		  "(c) 2014 Christian Fuchsberger, Goncalo Abecasis, David Hinds\n");
#ifdef VERSION
   printf("RELEASE STAMP " VERSION "\n");
#else
   printf("UNDOCUMENTED RELEASE\n");
#endif

   int rounds = 5, states = 200, cpus = 0;
   bool em = false, gzip = false, info_only = false;
   bool phased = false, alleles = false, probs = false;
   bool rs = false;
   bool mach = true;

   String referenceHaplotypes, referenceSnps, snpAliases;
   String haplotypes, shape_haplotypes, snps, sample, chr, exclude;
   String prefix("pre");
   String firstMarker, lastMarker;
   String clipFile;

   String recombinationRates, errorRates;

   BEGIN_LONG_PARAMETERS(longParameters)
	  LONG_PARAMETER_GROUP("Reference Haplotypes")
		 LONG_STRINGPARAMETER("refHaps", &referenceHaplotypes)
		 LONG_STRINGPARAMETER("refSnps", &referenceSnps)
		 LONG_STRINGPARAMETER("snpAliases", &snpAliases)
		 LONG_PARAMETER("vcfReference", &HaplotypeSet::vcfReference)
	  LONG_PARAMETER_GROUP("Target Haplotypes (MaCH)")
		 LONG_STRINGPARAMETER("haps", &haplotypes)
		 LONG_STRINGPARAMETER("snps", &snps)
		 LONG_PARAMETER("rs", &rs)
	  LONG_PARAMETER_GROUP("Target Haplotypes (ShapeIT)")
		 LONG_STRINGPARAMETER("shape_haps", &shape_haplotypes)
		 LONG_STRINGPARAMETER("sample", &sample)
		 LONG_STRINGPARAMETER("chr", &chr)
		 LONG_PARAMETER("rs", &rs)
	  LONG_PARAMETER_GROUP("Starting Parameters")
		 LONG_STRINGPARAMETER("rec", &recombinationRates)
		 LONG_STRINGPARAMETER("erate", &errorRates)
	  LONG_PARAMETER_GROUP("Parameter Fitting")
		 LONG_INTPARAMETER("rounds", &rounds)
		 LONG_INTPARAMETER("states", &states)
		 LONG_PARAMETER("em", &em)
	  LONG_PARAMETER_GROUP("Output Files")
		 LONG_STRINGPARAMETER("prefix", &prefix)
		 LONG_PARAMETER("phased", &phased)
		 LONG_PARAMETER("alleles", &alleles)
		 LONG_PARAMETER("probs", &probs)
		 LONG_PARAMETER("gzip", &gzip)
	  LONG_PARAMETER_GROUP("Regional Options")
		 LONG_DOUBLEPARAMETER("vcfstart", &HaplotypeSet::startposition)
		 LONG_DOUBLEPARAMETER("vcfend", &HaplotypeSet::endposition)
         LONG_DOUBLEPARAMETER("vcfwindow", &HaplotypeSet::window)
	  LONG_PARAMETER_GROUP("Clipping Window")
		 LONG_STRINGPARAMETER("start", &firstMarker)
		 LONG_STRINGPARAMETER("stop", &lastMarker)
		 LONG_STRINGPARAMETER("autoClip", &clipFile)
	  LONG_PARAMETER_GROUP("Experimental")
		 LONG_PARAMETER("info_only", &info_only)
 	     LONG_STRINGPARAMETER("exclude", &exclude)
	     LONG_STRINGPARAMETER("mimicArray", &HaplotypeSet::mimicArray)
	     LONG_STRINGPARAMETER("vcfchr", &HaplotypeSet::schr)
	   
#ifdef _OPENMP
	  LONG_PARAMETER_GROUP("Multi-Threading")
		 LONG_INTPARAMETER("cpus", &cpus)
#endif
   END_LONG_PARAMETERS();

   ParameterList pl;

   pl.Add(new LongParameters("Command Line Options", longParameters));
   pl.Read(argc, argv);
   pl.Status();

#ifdef _OPENMP
   if (cpus > 0)
	  omp_set_num_threads(cpus);
#endif

   HaplotypeSet target;
   StringArray markerList;

  // Load framework haplotypes
   if(HaplotypeSet::window > 0 && rs)
	{
		   error("Framework chunking does not support rs numbers!");
	}
	   
   // Load framework haplotypes
   if(!haplotypes.IsEmpty() && !snps.IsEmpty())
    {
	   printf("Loading MaCH target haplotypes ...\n");

	   markerList.Read(snps);
	   if (markerList.Length() == 0)
	     error("Framework marker list is empty - please verify filenames and contents!");

	   target.markerCount = markerList.Length();
	   target.LoadHaplotypes(haplotypes, true, false);
    }
   else
    {
	   if(!shape_haplotypes.IsEmpty() && !sample.IsEmpty())
	    {
	    	if(chr.IsEmpty() && rs==false)
	    		error("Please provide chr!");
		    mach = false;
	    	printf("Loading ShapeIT target haplotypes ...\n");
	    	target.LoadSamples(sample);
            target.LoadShapeITHaplotypes(shape_haplotypes,markerList,rs,chr);
	    }
	    else
	    {
		    error("Missing target haplotypes - please verify filenames and contents!");
	    }
    }

   printf("  %d Target Haplotypes Loaded ...\n\n", target.count);

   // regional imputation only works with vcf input haplotypes, can easily extend to HapMap format: HaplotypeSet::LoadHapMapHaplotypes()
   if (!HaplotypeSet::vcfReference && (HaplotypeSet::startposition != 0.0 || HaplotypeSet::endposition != 300000000.0))
	 error("--vcfstart and --vcfend options for regional imputation only work vcf input reference haplotypes.\n\n");
   if(HaplotypeSet::startposition < 0.0 || HaplotypeSet::endposition > 300000000.0)
	 error("--vcfstart or --vcfend are out of range.\n\n");

	StringAlias	* aliases = new StringAlias;

	if (snpAliases.Length())
		{
		// Read alternate names for SNPs
		printf("Loading SNP Aliases List ...\n");

		aliases->ReadFromFile(snpAliases);
		}

     // Read marker list
  printf("Reading Reference Marker List ...\n");

  StringArray refMarkerList;
  int countFromVcf = 0;

  if (HaplotypeSet::vcfReference)
   {

	  //set new start/stop based on window
	   if (HaplotypeSet::window > 0)
	      {

	 		  if (HaplotypeSet::startposition-HaplotypeSet::window < 0)
	 			 HaplotypeSet::startposition = 0;
	 		  else
	 			 HaplotypeSet::startposition -= HaplotypeSet::window;

	 		 HaplotypeSet::endposition += HaplotypeSet::window;
	      }
	    countFromVcf = HaplotypeSet::LoadReferenceSNPsFromVcf(referenceHaplotypes, refMarkerList, rs, exclude);
		printf("Read %d reference haplotypes from vcf file\n", countFromVcf);
   }
  else
	 refMarkerList.Read(referenceSnps);

  int changedLabels = aliases->GetAliases(refMarkerList);

  if (changedLabels)
	printf("   Changed %d Marker Labels Based on Aliases List\n", changedLabels);

  // Index markers
  StringIntHash referenceHash;
  for (int i = 0; i < refMarkerList.Length(); i++)
	referenceHash.Add(refMarkerList[i].Trim(), i);

  printf("  %d Markers in Reference Haplotypes...\n\n", refMarkerList.Length());

   if (refMarkerList.Length() == 0)
      error("Reference panel marker list is empty - please verify filenames and contents!");

   // Load reference haplotypes
   printf("Loading reference haplotypes ...\n");

   HaplotypeSet reference;
   reference.markerCount = refMarkerList.Length();

   if (HaplotypeSet::vcfReference)
	   reference.count = countFromVcf;

   reference.LoadHaplotypes(referenceHaplotypes,false,true,exclude);

   printf("  %d Reference Haplotypes Loaded ...\n\n", reference.count);

   changedLabels = aliases->GetAliases(markerList);

   if (changedLabels)
      printf("   Changed %d Marker Labels Based on Aliases List\n", changedLabels);

   StringArray clipInfo;
   clipInfo.Read(clipFile);
   String prelastMarker="";

   if (clipInfo.Length() && markerList.Length())
      {
      bool match = false, data = false;
      printf("   Reading clipping instructions from %s ...\n", (const char *) clipFile);

      StringArray tokens;
      for (int i = 0; i < clipInfo.Length(); i++)
         {
         tokens.ReplaceTokens(clipInfo[i]);

         if (tokens.Length() != 4) continue;

         tokens[0] = aliases->GetAlias(tokens[0]);
         tokens[1] = aliases->GetAlias(tokens[1]);

         if(tokens[2]=="start")
        	 data = true;

         if (tokens[0] != markerList[0].Trim() ||
             tokens[1] != markerList.Last().Trim())
         {
             if(data)
        	 prelastMarker = aliases->GetAlias(tokens[3]);
             //printf("%s\n",prelastMarker.c_str());
	         continue;
         }

         firstMarker = aliases->GetAlias(tokens[2]); // core start
         lastMarker = aliases->GetAlias(tokens[3]);  // core end

         match = true;
         break;
         }

      if(prelastMarker != "")
    	 firstMarker = prelastMarker;

      if (!match)
         printf("      No match found, automatic clipping disabled.\n");
      else
         printf("      Clipping will start at %s, end at %s.\n",
                (const char *) firstMarker, (const char *) lastMarker);
      }

   //ClipReference(reference, refMarkerList, referenceHash, markerList,
   //              firstMarker, lastMarker);

   // Done with aliases list, reclaim memory
   delete aliases;

   // Crossref Marker Names to Reference Panel Positions
   IntArray markerIndex;
   markerIndex.Dimension(markerList.Length());

   int matches = 0;

   for (int i = 0; i < markerList.Length(); i++)
      {
      markerIndex[i] = referenceHash.Integer(markerList[i].Trim());

      if (markerIndex[i] >= 0) matches++;
      }

   printf("  %d Markers in Framework Haplotypes Overlap Reference ...\n", matches);

   if (matches == 0)
      error("No markers overlap between target and reference\n"
            "Please check correct reference is being used and markers are named consistently");

   printf("  %d Other Markers in Framework Haplotypes Discarded ...\n\n", markerList.Length() - matches);

   // Check for flips in reference vs. target haplotypes
   int flips = 0;
   int previous = -1;
   for (int i = 0; i < markerIndex.Length(); i++)
      if (markerIndex[i] >= 0)
         if (markerIndex[i] < previous)
            {
            if (flips++ < 10)
               printf("  -> Marker %s precedes %s in reference, but follows it in target\n",
                     (const char *) refMarkerList[previous],
                     (const char *) markerList[i]);
            previous = markerIndex[i];
            }
   if (flips > 10)
      printf("  -> %d Additional Marker Order Changes Not Listed\n", flips - 10);
   if (flips)
      printf("  %d Marker Pairs Change Order in Target vs Framework Haplotypes\n", flips);


   //
   int startIndex = firstMarker.IsEmpty() ? 0 : referenceHash.Integer(firstMarker);
   int stopIndex = lastMarker.IsEmpty() ? reference.markerCount - 1 : referenceHash.Integer(lastMarker);

   // tweak index
   if(!lastMarker.IsEmpty())
     {
 	      stopIndex = stopIndex-1;
     }

   // Calculate frequencies mismatch
   reference.CalculateFrequencies();
   target.CalculateFrequencies();
   target.CompareFrequencies(reference, markerIndex, markerList, startIndex, stopIndex);

   if (startIndex < 0 || stopIndex < 0)
      error("Clipping requested, but no position available for one of the endpoints");

   // List the major allele at each location
   if (HaplotypeSet::hapmapFormat)
	reference.ListMajorAlleleStrings();
   else
   	reference.ListMajorAlleles();


   // define output indexes if vcfWindow is definied
   int outStartIndex = startIndex;
   int outStopIndex = stopIndex;

   if(HaplotypeSet::window > 0)
    {
	   String marker;
	   if(HaplotypeSet::corestartposition > 0)
	    {
		   marker.printf("%0.0f:%0.0f",HaplotypeSet::chr,HaplotypeSet::corestartposition);
		   //printf("Core-Start: %s\n",marker.c_str());
		   outStartIndex = referenceHash.Integer(marker);
	    }
	   marker.printf("%0.0f:%0.0f",HaplotypeSet::chr,HaplotypeSet::coreendposition);
	   outStopIndex = referenceHash.Integer(marker);
	   if(outStopIndex == -1 )
		   outStopIndex = stopIndex;
	   //printf("Core-Stop: %s\n",marker.c_str());
    }

   printf("Generating Draft .info File ...\n\n");
   // Output some basic information
   IFILE info = ifopen(prefix + ".info.draft", "wb");

   ifprintf(info, "SNP\tAl1\tAl2\tFreq1\tGenotyped\n");

  if (HaplotypeSet::hapmapFormat || HaplotypeSet::TransposeReference || HaplotypeSet::ImputeReference)
  {
    for (int i = 0, j = 0; i <= outStopIndex; i++)
      if (i >= outStartIndex)
      {
        String majorAllele, minorAllele;
        if (reference.major[i] == 1)
        {
          majorAllele = reference.AlleleOneLabelString[i];
          minorAllele = reference.AlleleTwoLabelString[i];
        }
        else
        {
          majorAllele = reference.AlleleTwoLabelString[i];
          minorAllele = reference.AlleleOneLabelString[i];
        }
        ifprintf(info, "%s\t%s\t%s\t%.4f\t%s\n",
          (const char *) refMarkerList[i],
          (const char *) majorAllele, (const char *) minorAllele,
          reference.freq[reference.major[i]][i],
          j < markerIndex.Length() && i == markerIndex[j] ? (j++, "Genotyped") : "-");
      }
      else
        if (j < markerIndex.Length() && i == markerIndex[j])
          j++;
  }
  else
  {
    for (int i = 0, j = 0; i <= outStopIndex; i++)
      if (i >= outStartIndex)
        ifprintf(info, "%s\t%s\t%s\t%.4f\t%s\n",
          (const char *) refMarkerList[i],
          reference.MajorAlleleLabel(i), reference.MinorAlleleLabel(i),
          reference.freq[reference.major[i]][i],
          j < markerIndex.Length() && i == markerIndex[j] ? (j++, "Genotyped") : "-");
      else
        if (j < markerIndex.Length() && i == markerIndex[j])
          j++;
  }

  ifclose(info);

  if(info_only) exit(1);

   
   printf("Setting up Markov Model...\n\n");

   // Setup Markov Model
   MarkovParameters mp;

   mp.Allocate(reference.markerCount);

   if (rounds > 0)
      printf("Initializing Model Parameters (using %s and up to %d haplotypes)\n",
             em ? "E-M" : "MCMC", states);

   // Simple initial estimates of error and recombination rate
   for (int i = 0; i < reference.markerCount; i++)
      mp.E[i] = 0.01;

   for (int i = 0; i < reference.markerCount - 1; i++)
      mp.R[i] = 0.001;

	   if(!errorRates.IsEmpty()){
		   if (mp.ReadErrorRates(errorRates))
			   printf("  Updated error rates using data in %s ...\n", (const char *) errorRates);
		   else
			   error("Cannot update error rate, please check file!");
	   }
	   
	   if(!recombinationRates.IsEmpty()){
		   if (mp.ReadCrossoverRates(recombinationRates))
			   printf("  Updated recombination rates using %s ...\n", (const char *) recombinationRates);
		   else
			   error("Cannot update recombination map, please check file!");
	   }
	   
   // Parameter estimation loop
   for (int round = 0; round < rounds; round++)
      {
      printf("  Round %d of Parameter Refinement ...\n", round + 1);

      int iterations = states < reference.count ? states : reference.count;

      MarkovModel original;
      original.CopyParameters(mp);

      #pragma omp parallel for
      for (int i = 0; i < iterations; i++)
         {
         MarkovModel mm;

         mm.Allocate(reference.markerCount, reference.count - 1);
         mm.CopyParameters(original);

         // Reference leave one out (loo) panel
         char ** reference_loo = new char * [reference.count - 1];
         for (int in = 0, out = 0; in < reference.count; in++)
            if (in != i)
               reference_loo[out++] = reference.haplotypes[in];

         mm.WalkLeft(reference.haplotypes[i], reference_loo, reference.freq, false);

         if (em)
            mm.CountExpected(reference.haplotypes[i], reference_loo, reference.freq);
         else
            {
            #pragma omp critical
            { mm.ProfileModel(reference.haplotypes[i], reference_loo, reference.freq); }
            }

         delete [] reference_loo;

         #pragma omp critical
         mp += mm;
         }

      if (round >= rounds / 2)
         {
         int iterations = states < target.count ? states : target.count;

         #pragma omp parallel for
         for (int i = 0; i < iterations; i++)
            {
            MarkovModel mm;

            mm.Allocate(reference.markerCount, reference.count);
            mm.CopyParameters(original);

            // Padded version of target haplotype, including missing sites
            char * padded = new char [reference.markerCount];
            for (int k = 0; k < reference.markerCount; k++)
               padded[k] = 0;

            // Copy current haplotype into padded vector
            for (int j = 0; j < target.markerCount; j++)
               if (markerIndex[j] >= 0)
                  padded[markerIndex[j]] = target.haplotypes[i][j];

            mm.WalkLeft(padded, reference.haplotypes, reference.freq, false);

            if (em)
               mm.CountExpected(padded, reference.haplotypes, reference.freq);
            else
               {
               #pragma omp critical
               { mm.ProfileModel(padded, reference.haplotypes, reference.freq); }
               }

            delete [] padded;

            #pragma omp critical
            mp += mm;
            }
         }

      mp.UpdateModel();

      double crossovers = 0;
      for (int i = 0; i < reference.markerCount - 1; i++)
         crossovers += mp.R[i];

	  double errors = 0;
	  for (int i = 0; i < reference.markerCount; i++)
		 {
		 double heterozygosity = 1.0 - square(reference.freq[1][i])
									 - square(reference.freq[2][i])
									 - square(reference.freq[3][i])
									 - square(reference.freq[4][i])
									 - square(reference.freq[5][i])
									 - square(reference.freq[6][i])
									 - square(reference.freq[7][i]);


		 errors += mp.E[i] * heterozygosity;
		 }
	  errors /= reference.markerCount + 1e-30;

	  printf("      %.1f mosaic crossovers expected per haplotype\n", crossovers);
	  printf("      %.1f%% of crossovers are due to reference flips\n", mp.empiricalFlipRate * 100.);
	  printf("      %.3g errors in mosaic expected per marker\n", errors);
	  }

   if (rounds > 0)
	  {
	  printf("  Saving estimated parameters for future use ...\n");
	  mp.WriteParameters(refMarkerList, prefix, gzip);
	  }

   printf("\n");

   reference.Transpose();

   time_t istart = log_now("Setup completed", start);
   printf("Imputing Genotypes ...\n");

   IFILE dosages, probabilities, hapdose, haps;

   if (!phased)
      {
      dosages = ifopen(prefix + ".dose" + (gzip ? ".gz" : ""), "wb");
      
      if (dosages == NULL)
	 error("Error opening output file for storing dosage information.");
      }

   if (phased)
      {
      hapdose = ifopen(prefix + ".hapDose" + (gzip ? ".gz" : ""), "wb");

      if (hapdose == NULL)
         error("Error opening output files for storing phased dosage information.");
      }


   if (alleles)
      {
      haps = ifopen(prefix + ".haps" + (gzip ? ".gz" : ""), "wb");

      if (haps == NULL)
         error("Error opening output files for storing haplotype information.");
      }


   if (probs)
      {
      probabilities = ifopen(prefix + ".prob" + (gzip ? ".gz" : ""), "wb");

      if (probabilities == NULL)
         error("Error opening output file for storing genotype probabilities.");
      }

   ImputationStatistics stats(reference.markerCount);

   // Impute each haplotype
   #pragma omp parallel for
   for (int i = 0; i < target.count; i++)
      {
      if (i != 0 && target.labels[i] == target.labels[i-1])
         continue;

      MarkovModel mm;

      mm.Allocate(reference.markerCount, reference.count);
      mm.ClearImputedDose();
      mm.CopyParameters(mp);

      // Padded version of target haplotype, including missing sites
      char * padded = new char [reference.markerCount];
      for (int j = 0; j < reference.markerCount; j++)
         padded[j] = 0;

      int k = i;

      do {
         printf("  Processing Haplotype %d of %d ...\n", k + 1, target.count);

         // Copy current haplotype into padded vector
         for (int j = 0; j < target.markerCount; j++)
            if (markerIndex[j] >= 0)
               padded[markerIndex[j]] = target.haplotypes[k][j];

         mm.WalkLeft(padded, reference.haplotypes, reference.freq, true);
         mm.Impute(reference.major, padded, reference.haplotypes, reference.freq);

         #pragma omp critical
         { stats.Update(mm.imputedHap, mm.leaveOneOut, padded, reference.major); }

         #pragma omp critical
	    {
	    if (phased)
	       {
	       ifprintf(hapdose, "%s\tHAPLO%d", (const char *) target.labels[i], k - i + 1);
	       for (int j = outStartIndex; j <= outStopIndex; j++)
		  ifprintf(hapdose, "\t%.3f", mm.imputedHap[j]);
	       ifprintf(hapdose, "\n");
	       }
	    if (alleles)
	       {
	       ifprintf(haps, "%s\tHAPLO%d\t", (const char *) target.labels[i], k - i + 1);
	       for (int j = outStartIndex; j <= outStopIndex; j++)
		  //ifprintf(haps, "%s%c", j % 8 == 0 ? " " : "", mm.imputedAlleles[j]);
		  ifprintf(haps, "%c", mm.imputedAlleles[j]);
	       ifprintf(haps, "\n");
	       }
	    }
	 
         k++;
      } while (k < target.count && target.labels[k] == target.labels[i]);

      printf("    Outputting Individual %s ...\n", (const char *) target.labels[i]);

      #pragma omp critical
	 {
	 if (!phased)
	    {
	    ifprintf(dosages, "%s\tDOSE", (const char *) target.labels[i]);
	    for (int j = outStartIndex; j <= outStopIndex; j++)
	       ifprintf(dosages, "\t%.3f", mm.imputedDose[j]);
	    ifprintf(dosages, "\n");
	    }
	 if (probs)
            {
            ifprintf(probabilities, "%s\tPROB", (const char *) target.labels[i]);
            for (int j = outStartIndex; j <= outStopIndex; j++)
               {
               double p1 = mm.imputedDose[j] - mm.imputedHap[j];
               double p2 = mm.imputedHap[j];

               if (p1 < 0.) p1 = 0.;

               double pHom = p1 * p2;
               double pHet = p1 * (1. - p2) + (1. - p1) * p2;

               ifprintf(probabilities, "\t%.3f\t%.3f", pHom, pHet);
               }
            ifprintf(probabilities, "\n");
            }
         }

      delete [] padded;
      }

   if (!phased)
      ifclose(dosages);
   if (phased)
      ifclose(hapdose);
   if (alleles)
      ifclose(haps);
   if (probs)
      ifclose(probabilities);

   // Output some basic information
   info = ifopen(prefix + ".info" + (gzip ? ".gz" : ""), "wb");

   ifprintf(info, "SNP\tAl1\tAl2\tFreq1\tMAF\tAvgCall\tRsq\tGenotyped\tLooRsq\tEmpR\tEmpRsq\tDose1\tDose2\n");

   // Padded version of target haplotype, including missing sites
   char * padded = new char [reference.markerCount];
   for (int k = 0; k < reference.markerCount; k++)
      padded[k] = 0;

   // Mark genotyped SNPs in padded vector
   for (int j = 0; j < target.markerCount; j++)
      if (markerIndex[j] >= 0)
          padded[markerIndex[j]] = 1;

   for (int i = outStartIndex; i <= outStopIndex; i++)
      if (reference.MajorAlleleLabel(i)[0] == reference.MinorAlleleLabel(i)[0])
         {
         reference.MinorAlleleLabel(i)[0];
         }


   if (HaplotypeSet::hapmapFormat || HaplotypeSet::TransposeReference || HaplotypeSet::ImputeReference)
     for (int i = outStartIndex; i <= outStopIndex; i++)
      {
      String majorAllele, minorAllele;
      if (reference.major[i] == 1)
      {
        majorAllele = reference.AlleleOneLabelString[i];
        minorAllele = reference.AlleleTwoLabelString[i];
      }
      else
      {
        majorAllele = reference.AlleleTwoLabelString[i];
        minorAllele = reference.AlleleOneLabelString[i];
      }
      ifprintf(info, "%s\t%s\t%s\t%.5f\t%.5f\t%.5f\t%.5f\t",
            (const char *) refMarkerList[i],
            (const char *) majorAllele,
            (const char *) minorAllele,
            stats.AlleleFrequency(i),
            stats.AlleleFrequency(i) > 0.5 ? 1.0 - stats.AlleleFrequency(i) : stats.AlleleFrequency(i),
            stats.AverageCallScore(i),
            stats.Rsq(i));

      if (padded[i])
         ifprintf(info, "Genotyped\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n",
                  stats.LooRsq(i), stats.EmpiricalR(i), stats.EmpiricalRsq(i),
                  stats.LooMajorDose(i), stats.LooMinorDose(i));
      else
         ifprintf(info, "-\t-\t-\t-\t-\t-\n");
	  }
   else
	   for (int i = 0; i <= outStopIndex; i++)
	   {
		 if (i >= outStartIndex){
		  ifprintf(info, "%s\t%s\t%s\t%.5f\t%.5f\t%.5f\t%.5f\t",
				(const char *) refMarkerList[i],
				reference.MajorAlleleLabel(i),
				reference.MinorAlleleLabel(i),
				stats.AlleleFrequency(i),
				stats.AlleleFrequency(i) > 0.5 ? 1.0 - stats.AlleleFrequency(i) : stats.AlleleFrequency(i),
				stats.AverageCallScore(i),
				stats.Rsq(i));

		  if (padded[i])
			 ifprintf(info, "Genotyped\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n",
					  stats.LooRsq(i), stats.EmpiricalR(i), stats.EmpiricalRsq(i),
					  stats.LooMajorDose(i), stats.LooMinorDose(i));
		  else
			 ifprintf(info, "-\t-\t-\t-\t-\t-\n");
		   }
	   }

   ifclose(info);

   delete [] padded;

   log_now("Imputation completed", istart);
   log_now("Run completed", start);
   }


/* Local variables: */
/* c-file-style: "minimac" */
/* End: */
