package com.hartwig.hmftools.bamtools.checker;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MISMATCHES_AND_DELETIONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;

import java.util.concurrent.atomic.LongAdder;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

// recalculates NM, and MD if present, for mapped reads against the ref genome, since a BAM aligned to a different build of the
// assembly has these tags relative to that build, eg GRC-masked regions which are N in the ref genome
// one instance per thread since the ref genome reader isn't thread-safe; ref bases are cached in a window around the current reads
public class AlignmentTagCalculator
{
    private final RefGenomeInterface mRefGenome;
    private final int mWindowSize;

    private String mWindowChromosome;
    private int mWindowStart; // 1-based position of the first cached base
    private byte[] mWindowBases;

    private final LongAdder mRecalculated;
    private final LongAdder mNmChanged;

    private static final int DEFAULT_WINDOW_SIZE = 100_000;

    public AlignmentTagCalculator(final RefGenomeInterface refGenome)
    {
        this(refGenome, DEFAULT_WINDOW_SIZE);
    }

    @VisibleForTesting
    public AlignmentTagCalculator(final RefGenomeInterface refGenome, int windowSize)
    {
        mRefGenome = refGenome;
        mWindowSize = windowSize;
        mWindowChromosome = null;
        mWindowStart = 0;
        mWindowBases = null;
        mRecalculated = new LongAdder();
        mNmChanged = new LongAdder();
    }

    public long recalculatedCount() { return mRecalculated.sum(); }
    public long nmChangedCount() { return mNmChanged.sum(); }

    public void recalculate(final SAMRecord record)
    {
        if(record.getReadUnmappedFlag())
            return;

        String chromosome = record.getReferenceName();
        int alignmentStart = record.getAlignmentStart();
        int alignmentEnd = record.getAlignmentEnd();

        if(alignmentEnd > mRefGenome.getChromosomeLength(chromosome))
            return; // alignment beyond the contig end, so leave the tags as they are

        if(record.getCigar().getReadLength() != record.getReadLength())
            return; // no or truncated bases (eg SEQ *), so the tags cannot be recalculated

        if(!loadWindow(chromosome, alignmentStart, alignmentEnd))
            return;

        int numMutations = 0;
        int matchCount = 0;
        StringBuilder mdString = new StringBuilder();

        byte[] readBases = record.getReadBases();
        int readIndex = 0;
        int refIndex = alignmentStart - mWindowStart;

        for(CigarElement element : record.getCigar().getCigarElements())
        {
            int length = element.getLength();

            switch(element.getOperator())
            {
                case M:
                case EQ:
                case X:
                    for(int i = 0; i < length; ++i)
                    {
                        byte refBase = mWindowBases[refIndex + i];

                        if(SequenceUtil.basesEqual(readBases[readIndex + i], refBase))
                        {
                            ++matchCount;
                        }
                        else
                        {
                            mdString.append(matchCount);
                            mdString.append((char)SequenceUtil.upperCase(refBase));
                            matchCount = 0;
                            ++numMutations;
                        }
                    }

                    readIndex += length;
                    refIndex += length;
                    break;

                case D:
                    mdString.append(matchCount).append('^');

                    for(int i = 0; i < length; ++i)
                    {
                        mdString.append((char)SequenceUtil.upperCase(mWindowBases[refIndex + i]));
                    }

                    matchCount = 0;
                    numMutations += length;
                    refIndex += length;
                    break;

                case I:
                    numMutations += length;
                    readIndex += length;
                    break;

                case S:
                    readIndex += length;
                    break;

                case N:
                    refIndex += length;
                    break;

                default: // H and P consume neither read bases nor ref bases
                    break;
            }
        }

        mdString.append(matchCount);

        mRecalculated.increment();

        Integer existingNm = record.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE);

        if(existingNm == null || existingNm != numMutations)
            mNmChanged.increment();

        record.setAttribute(NUM_MUTATONS_ATTRIBUTE, numMutations);

        // MD is optional, so only rewritten if it was present
        if(record.hasAttribute(MISMATCHES_AND_DELETIONS_ATTRIBUTE))
            record.setAttribute(MISMATCHES_AND_DELETIONS_ATTRIBUTE, mdString.toString());
    }

    private boolean loadWindow(final String chromosome, int alignmentStart, int alignmentEnd)
    {
        if(mWindowBases != null && chromosome.equals(mWindowChromosome)
        && alignmentStart >= mWindowStart && alignmentEnd < mWindowStart + mWindowBases.length)
        {
            return true;
        }

        int windowEnd = Math.min(mRefGenome.getChromosomeLength(chromosome), Math.max(alignmentEnd, alignmentStart + mWindowSize - 1));

        byte[] bases = mRefGenome.getBases(chromosome, alignmentStart, windowEnd);

        if(bases == null || bases.length < alignmentEnd - alignmentStart + 1)
            return false;

        mWindowChromosome = chromosome;
        mWindowStart = alignmentStart;
        mWindowBases = bases;
        return true;
    }
}
