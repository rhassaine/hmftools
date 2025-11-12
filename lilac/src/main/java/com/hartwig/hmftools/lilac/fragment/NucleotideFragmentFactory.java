package com.hartwig.hmftools.lilac.fragment;

import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.lilac.LilacUtils.calcNucelotideLocus;
import static com.hartwig.hmftools.lilac.LilacUtils.formRange;
import static com.hartwig.hmftools.lilac.ReferenceData.GENE_CACHE;
import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.calcAminoAcidIndices;
import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.expandIndices;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.DEL_STR;

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.Codons;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.utils.SuffixTree;
import com.hartwig.hmftools.lilac.ReferenceData;
import com.hartwig.hmftools.lilac.evidence.Nucleotide;
import com.hartwig.hmftools.lilac.hla.HlaGene;
import com.hartwig.hmftools.lilac.read.Indel;
import com.hartwig.hmftools.lilac.read.Read;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

public class NucleotideFragmentFactory
{
    private final ReferenceData mReferenceData;
    private final LinkedHashMap<HlaSequenceLoci, SuffixTree> mInsertSuffixTrees;
    private final LinkedHashMap<HlaSequenceLoci, SuffixTree> mDeleteSuffixTrees_;
    private final byte mMinBaseQuality;

    public NucleotideFragmentFactory(final ReferenceData referenceData)
    {
        mMinBaseQuality = LOW_BASE_QUAL_THRESHOLD;
        mReferenceData = referenceData;

        mInsertSuffixTrees = Maps.newLinkedHashMap();
        mDeleteSuffixTrees_ = Maps.newLinkedHashMap();

        mReferenceData.AminoAcidSequencesWithInserts.forEach(x -> mInsertSuffixTrees.put(x, new SuffixTree(x.sequence())));
        mReferenceData.AminoAcidSequencesWithDeletes_.forEach(x -> mDeleteSuffixTrees_.put(x, new SuffixTree(x.sequence())));
    }

    public Fragment createFragment(final Read read, final HlaGene geneName, final byte geneStrand)
    {
        if(read.ReadIndexStart < 0 || read.ReadIndexEnd < read.ReadIndexStart || read.PositionEnd < read.PositionStart)
            return null;

        boolean reverseStrand = geneStrand == NEG_STRAND;

        int codingPositionStartLoci = calcNucelotideLocus(GENE_CACHE.Transcripts, read.PositionStart);
        int codingPositionEndLoci = calcNucelotideLocus(GENE_CACHE.Transcripts, read.PositionEnd);

        int samCodingStartLoci = !reverseStrand ? codingPositionStartLoci : codingPositionEndLoci;
        int samCodingEndLoci = !reverseStrand ? codingPositionEndLoci : codingPositionStartLoci;

        int readLength = read.ReadIndexEnd - read.ReadIndexStart + 1;
        final char[] codingRegionReadBases = new char[readLength];
        final byte[] codingRegionQualities = new byte[readLength];

        read.populateCodingRegion(codingRegionReadBases, codingRegionQualities, reverseStrand);

        if(read.containsValidIndel_() || read.containsSoftClip())
        {
            List<Integer> aminoAcidIndices = calcAminoAcidIndices(samCodingStartLoci, samCodingEndLoci);
            int firstAAIndex = aminoAcidIndices.get(0);
            int nucleotideStartLoci = firstAAIndex * 3;
            String sequence = String.valueOf(codingRegionReadBases);
            int startLoci = nucleotideStartLoci - samCodingStartLoci;

            if(startLoci < 0 || startLoci >= sequence.length())
            {
                // likely due to a delete in this region
                LL_LOGGER.trace("invalid startLoci({}) requested: read({}) gene({})", startLoci, read.Id, geneName);
                return null;
            }
            String aminoAcids = Codons.aminoAcidFromBases(sequence.substring(startLoci));

            if(!aminoAcids.isEmpty())
            {
                int matchRangeAllowedStart_ = firstAAIndex - read.SoftClippedStart / 3 - read.maxValidIndelSize_();
                int matchRangeAllowedEnd_ = firstAAIndex + read.maxValidIndelSize_() + read.SoftClippedEnd / 3;

                Fragment matchedFragment = checkMatchedInsertDeleteSequence__(
                        read, geneName, aminoAcids, matchRangeAllowedStart_, matchRangeAllowedEnd_, mInsertSuffixTrees);

                if(matchedFragment != null)
                    return matchedFragment;

                matchedFragment = checkMatchedInsertDeleteSequence__(
                        read, geneName, aminoAcids, matchRangeAllowedStart_, matchRangeAllowedEnd_, mDeleteSuffixTrees_);

                if(matchedFragment != null)
                    return matchedFragment;
            }

            if(read.containsValidIndel_())
            {
                // TODO: remove this
                StringJoiner logMsgBuilder = new StringJoiner(",");
                logMsgBuilder.add("*** NON-IGNORED");
                logMsgBuilder.add(read.bamRecord().getReadName());
                logMsgBuilder.add(read.bamRecord().getContig());
                logMsgBuilder.add(String.valueOf(read.bamRecord().getAlignmentStart()));
                for(Indel indel : read.getValidIndels_())
                    logMsgBuilder.add(indel.toString());

                System.out.println(logMsgBuilder.toString());

                return null;
            }
        }

        if(samCodingStartLoci < 0 || samCodingEndLoci < 0)
            return null;

        int rangeLength = samCodingEndLoci - samCodingStartLoci + 1;
        List<Nucleotide> nucleotides = Lists.newArrayListWithCapacity(rangeLength);

        final int inc = reverseStrand ? -1 : 1;
        final int iStart = reverseStrand ? readLength - 1 : 0;
        int samCodingLoci = reverseStrand ? samCodingEndLoci : samCodingStartLoci;
        int readIndex = read.ReadIndexStart;
        for(int i = iStart; i >= 0 && i < readLength;)
        {
            // TODO: clean this up.
            Indel currentIndel = null;
            for(Indel indel : read.getIgnoredIndels_())
            {
                if(indel.ReadIndex == readIndex)
                {
                    currentIndel = indel;
                    break;
                }
            }

            if(currentIndel == null)
            {
                nucleotides.add(new Nucleotide(
                        samCodingLoci, codingRegionQualities[i], String.valueOf(codingRegionReadBases[i])));

                samCodingLoci += inc;
                readIndex++;
                i += inc;
                continue;
            }

            if(currentIndel.IsInsert)
            {
                nucleotides.add(new Nucleotide(
                        samCodingLoci, codingRegionQualities[i], String.valueOf(codingRegionReadBases[i])));

                samCodingLoci += inc;
                readIndex += 1 + currentIndel.Length;
                i += inc * (1 + currentIndel.Length);
                continue;
            }

            // currentIndel is del
            nucleotides.add(new Nucleotide(
                    samCodingLoci, codingRegionQualities[i], String.valueOf(codingRegionReadBases[i])));
            samCodingLoci += inc;
            readIndex++;
            i += inc;

            for(int j = 1; j < currentIndel.Ref.length(); j++)
            {
                char base = currentIndel.Ref.charAt(j);
                if(reverseStrand)
                    base = Nucleotides.swapDnaBase(base);

                nucleotides.add(new Nucleotide(
                        samCodingLoci, LOW_BASE_QUAL_THRESHOLD, String.valueOf(base)));

                samCodingLoci += inc;
            }
        }

        if(reverseStrand)
            Collections.reverse(nucleotides);

        return new Fragment(read, geneName, Sets.newHashSet(geneName), nucleotides);
    }

    private Fragment checkMatchedInsertDeleteSequence__(
            final Read record_, final HlaGene geneName,
            final String aminoAcids_, int matchRangeAllowedStart, int matchRangeAllowedEnd,
            final LinkedHashMap<HlaSequenceLoci, SuffixTree> sequenceMap_)
    {
        List<List<Integer>> matchedIndicesList_ = Lists.newArrayList();
        List<HlaSequenceLoci> matchedSeqLoci_ = Lists.newArrayList();

        for(Map.Entry<HlaSequenceLoci, SuffixTree> entry_ : sequenceMap_.entrySet())
        {
            HlaSequenceLoci seqLoci_ = entry_.getKey();
            List<Integer> aaIndices_ = entry_.getValue().indices(aminoAcids_);

            List<Integer> filteredAaIndices_ = aaIndices_.stream()
                    .filter(x -> x >= matchRangeAllowedStart && x <= matchRangeAllowedEnd)
                    .collect(Collectors.toList());

            if(!filteredAaIndices_.isEmpty())
            {
                matchedSeqLoci_.add(seqLoci_);
                matchedIndicesList_.add(filteredAaIndices_);
            }
        }

        for(int i = 0; i < matchedSeqLoci_.size(); ++i)
        {
            HlaSequenceLoci seqLoci_ = matchedSeqLoci_.get(i);
            List<Integer> filteredAaIndices_ = matchedIndicesList_.get(i);
            Fragment fragment = createIndelFragment__(record_, geneName, filteredAaIndices_.get(0), aminoAcids_, seqLoci_);
            if(!fragment.nucleotidesByLoci().isEmpty())
                return fragment;
        }

        return null;
    }

    private Fragment createIndelFragment__(
            final Read record_, final HlaGene geneName, final int startLoci, final String bamSequence__, final HlaSequenceLoci hlaSequence)
    {
        int endLoci = endLoci(startLoci, bamSequence__, hlaSequence);
        List<Integer> aminoAcidLoci = formRange(startLoci, endLoci);
        List<Integer> nucleotideLoci = expandIndices(aminoAcidLoci);

        List<String> nucleotides = Lists.newArrayList();

        aminoAcidLoci.stream()
                .map(hlaSequence::sequence)
                .map(NucleotideFragmentFactory::createNucleotidesFromAminoAcid)
                .forEach(nucleotides::addAll);

        List<Byte> qualities = nucleotideLoci.stream().map(x -> mMinBaseQuality).collect(Collectors.toList());

        return new Fragment(record_, geneName, Sets.newHashSet(geneName), nucleotideLoci, qualities, nucleotides);
    }

    private static int endLoci(int startLoci, final String bamSequence, final HlaSequenceLoci hlaSequence)
    {
        String tmpSequence = "";

        for(int locus = startLoci; locus < hlaSequence.length(); ++locus)
        {
            tmpSequence += hlaSequence.getSequences().get(locus);

            if(!bamSequence.startsWith(tmpSequence))
                return locus - 1;
        }

        return hlaSequence.length() - 1;
    }

    public static List<String> createNucleotidesFromAminoAcid(final String aminoAcid)
    {
        if(aminoAcid.equals(DEL_STR))
            return Lists.newArrayList(DEL_STR, DEL_STR, DEL_STR);

        String codons = Codons.aminoAcidsToCodons(aminoAcid);
        return Lists.newArrayList(
                String.valueOf(codons.charAt(0)), String.valueOf(codons.charAt(1)), codons.substring(2));
    }

    public Fragment createAlignmentFragments(final Read record, final HlaGene geneName, final byte geneStrand)
    {
        List<Fragment> fragments = record.alignmentsOnly().stream()
                .filter(Objects::nonNull)
                .map(x -> createFragment(record, geneName, geneStrand))
                .filter(Objects::nonNull)
                .collect(Collectors.toList());

        if(fragments.isEmpty())
            return null;

        if(fragments.size() == 1)
            return fragments.get(0);

        return FragmentUtils.mergeFragmentsById(fragments).get(0);
    }

    public byte calculatePercentileBaseQuality(final Iterable<Fragment> fragments, double percentile)
    {
        // calculates the nth percentile base quality for each fragment's nucleotides
        int maxBaseQual = mMinBaseQuality * 2;
        int[] baseQualFrequeny = new int[maxBaseQual + 1];
        long totalBases = 0;

        for(Fragment fragment : fragments)
        {
            for(byte baseQual : Nucleotide.qualities(fragment.rawNucleotidesByLoci().values()))
            {
                ++totalBases;
                ++baseQualFrequeny[min(baseQual, maxBaseQual)];
            }
        }

        long percentileEntry = round(totalBases * percentile);
        long cumulativeTotal = 0;

        for(int i = 0; i < baseQualFrequeny.length; ++i)
        {
            cumulativeTotal += baseQualFrequeny[i];

            if(cumulativeTotal >= percentileEntry)
                return (byte) i;
        }

        return mMinBaseQuality;
    }

    public static Map<HlaGene, int[]> calculateGeneCoverage(final Iterable<Fragment> fragments)
    {
        final Map<HlaGene, int[]> geneBaseDepth = Maps.newHashMap();

        for(HlaGene geneName : GENE_CACHE.GeneNames)
        {
            int geneNucleotideCount = GENE_CACHE.NucleotideLengths.get(geneName);
            geneBaseDepth.put(geneName, new int[geneNucleotideCount]);
        }

        for(Fragment fragment : fragments)
        {
            int[] baseDepth = geneBaseDepth.get(fragment.readGene());

            for(int locus : fragment.rawNucleotidesByLoci().keySet())
            {
                if(locus < baseDepth.length)
                    ++baseDepth[locus];
            }
        }

        return geneBaseDepth;
    }

}
