package com.hartwig.hmftools.bamtools.checker;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.MISMATCHES_AND_DELETIONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.XS_ATTRIBUTE;

import java.util.List;
import java.util.Set;
import java.util.concurrent.atomic.LongAdder;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

// supports writing the output BAM with the ref genome's sequence dictionary, eg for a BAM aligned to GRCh38.d1.vd1 with viral contigs
// the BAM writer resolves each record's contig and mate contig by name against its header, so records must not refer to any contig
// absent from the ref genome: primaries on such contigs are unmapped (placed at their mate if it is on a ref contig), supplementaries
// and secondaries are dropped, mate fields pointing at such contigs are set to mate-unmapped, and SA entries on them are removed
public class NonRefContigFilter
{
    private final Set<String> mRefContigs;

    private final LongAdder mDropped;
    private final LongAdder mUnmapped;
    private final LongAdder mMateUnmapped;
    private final LongAdder mSuppDataTrimmed;

    public NonRefContigFilter(final SAMSequenceDictionary refDictionary)
    {
        mRefContigs = Sets.newHashSet();
        refDictionary.getSequences().forEach(x -> mRefContigs.add(x.getSequenceName()));

        mDropped = new LongAdder();
        mUnmapped = new LongAdder();
        mMateUnmapped = new LongAdder();
        mSuppDataTrimmed = new LongAdder();
    }

    public static SAMFileHeader formOutputHeader(final SAMFileHeader inputHeader, final SAMSequenceDictionary refDictionary)
    {
        // keep read groups, programs and comments, and take contigs from the ref genome
        SAMFileHeader outputHeader = inputHeader.clone();
        outputHeader.setSequenceDictionary(refDictionary);
        return outputHeader;
    }

    public static boolean checkContigLengths(final SAMFileHeader inputHeader, final SAMSequenceDictionary refDictionary)
    {
        int bamOnlyCount = 0;
        boolean valid = true;

        for(SAMSequenceRecord bamSequence : inputHeader.getSequenceDictionary().getSequences())
        {
            SAMSequenceRecord refSequence = refDictionary.getSequence(bamSequence.getSequenceName());

            if(refSequence == null)
            {
                ++bamOnlyCount;
            }
            else if(refSequence.getSequenceLength() != bamSequence.getSequenceLength())
            {
                BT_LOGGER.error("contig({}) length differs: BAM({}) ref genome({})",
                        bamSequence.getSequenceName(), bamSequence.getSequenceLength(), refSequence.getSequenceLength());
                valid = false;
            }
        }

        BT_LOGGER.info("writing ref genome sequence dictionary: BAM contigs({}) ref genome contigs({}) BAM contigs not in ref genome({})",
                inputHeader.getSequenceDictionary().size(), refDictionary.size(), bamOnlyCount);

        return valid;
    }

    // returns false if the record should not be written
    public boolean apply(final SAMRecord record)
    {
        String refName = record.getReferenceName();

        if(!isRefContig(refName))
        {
            if(record.getSupplementaryAlignmentFlag() || record.isSecondaryAlignment())
            {
                mDropped.increment();
                return false;
            }

            unmapRecord(record);
        }

        if(record.getReadPairedFlag() && !isRefContig(record.getMateReferenceName()))
            setMateUnmapped(record);

        removeNonRefSuppData(record);

        return true;
    }

    private boolean isRefContig(final String contig)
    {
        return contig == null || NO_CHROMOSOME_NAME.equals(contig) || mRefContigs.contains(contig);
    }

    private void unmapRecord(final SAMRecord record)
    {
        mUnmapped.increment();

        record.setReadUnmappedFlag(true);
        record.setProperPairFlag(false);
        record.setInferredInsertSize(0);
        record.setMappingQuality(0);
        record.setCigarString(NO_CIGAR);
        record.setAttribute(SUPPLEMENTARY_ATTRIBUTE, null);
        record.setAttribute(NUM_MUTATONS_ATTRIBUTE, null);
        record.setAttribute(MISMATCHES_AND_DELETIONS_ATTRIBUTE, null);
        record.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, null);
        record.setAttribute(XS_ATTRIBUTE, null);

        String mateRefName = record.getMateReferenceName();

        if(record.getReadPairedFlag() && !record.getMateUnmappedFlag() && !NO_CHROMOSOME_NAME.equals(mateRefName)
        && isRefContig(mateRefName))
        {
            // place an unmapped read at its mapped mate
            int mateStart = record.getMateAlignmentStart();
            record.setReferenceName(mateRefName);
            record.setAlignmentStart(mateStart);
        }
        else
        {
            record.setReferenceName(NO_CHROMOSOME_NAME);
            record.setAlignmentStart(0);
        }
    }

    private void setMateUnmapped(final SAMRecord record)
    {
        // the mate is on a contig absent from the ref genome, so it has been or will be unmapped
        mMateUnmapped.increment();

        record.setMateUnmappedFlag(true);
        record.setProperPairFlag(false);
        record.setInferredInsertSize(0);
        record.setAttribute(MATE_CIGAR_ATTRIBUTE, null);

        String refName = record.getReferenceName();

        if(NO_CHROMOSOME_NAME.equals(refName))
        {
            record.setMateReferenceName(NO_CHROMOSOME_NAME);
            record.setMateAlignmentStart(0);
            return;
        }

        // an unmapped mate is placed at the primary of this read, which for a supplementary is its first SA entry
        String mateRefName = refName;
        int mateStart = record.getAlignmentStart();

        if(record.getSupplementaryAlignmentFlag())
        {
            SupplementaryReadData primaryData = SupplementaryReadData.extractAlignment(record);

            if(primaryData != null && isRefContig(primaryData.Chromosome))
            {
                mateRefName = primaryData.Chromosome;
                mateStart = primaryData.Position;
            }
        }

        record.setMateReferenceName(mateRefName);
        record.setMateAlignmentStart(mateStart);
    }

    private void removeNonRefSuppData(final SAMRecord record)
    {
        List<SupplementaryReadData> suppDataList = SupplementaryReadData.extractAlignments(record);

        if(suppDataList == null)
            return;

        List<SupplementaryReadData> refSuppDataList = suppDataList.stream()
                .filter(x -> isRefContig(x.Chromosome)).collect(Collectors.toList());

        if(refSuppDataList.size() == suppDataList.size())
            return;

        mSuppDataTrimmed.increment();

        record.setAttribute(
                SUPPLEMENTARY_ATTRIBUTE, refSuppDataList.isEmpty() ? null : SupplementaryReadData.alignmentsToSamTag(refSuppDataList));
    }

    public void logStats()
    {
        BT_LOGGER.info("non-ref contig records: dropped({}) unmapped({}) mateUnmapped({}) suppDataTrimmed({})",
                mDropped.sum(), mUnmapped.sum(), mMateUnmapped.sum(), mSuppDataTrimmed.sum());
    }
}
