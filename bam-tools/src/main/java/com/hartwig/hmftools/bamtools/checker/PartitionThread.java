package com.hartwig.hmftools.bamtools.checker;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_EXTENSION;

import java.io.File;
import java.util.List;
import java.util.NoSuchElementException;

import com.hartwig.hmftools.bamtools.common.PartitionTask;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.perf.TaskQueue;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class PartitionThread extends Thread
{
    private final SAMFileWriter mBamWriter;
    private final SamReader mSamReader;
    private final TaskQueue mPartitions;

    private final PartitionChecker mPartitionChecker;
    private final String mBamFilename;
    private final boolean mReverseBqsr;
    private final boolean mConvertRefGenome;
    private final AlignmentTagCalculator mAlignmentTagCalculator;
    private long mNonRefContigDropped;

    protected static final String UNSORTED_BAM_ID = "unsorted";
    protected static final String SORTED_BAM_ID = "sorted";

    public PartitionThread(
            final CheckConfig config, final FragmentCache fragmentCache, final TaskQueue partitions, final int threadId)
    {
        mConvertRefGenome = config.ConvertRefGenome;

        // NM and MD are recalculated against the ref genome, with a reader per thread
        mAlignmentTagCalculator = mConvertRefGenome && config.writeBam()
                ? new AlignmentTagCalculator(RefGenomeSource.loadRefGenome(config.RefGenomeFile)) : null;
        mNonRefContigDropped = 0;
        mPartitions = partitions;
        mReverseBqsr = config.ReverseBqsr;

        mSamReader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(config.RefGenomeFile))
                .validationStringency(config.BamStringency)
                .open(new File(config.BamFile));

        if(config.writeBam())
        {
            // create a BAM writer per thread
            mBamFilename = config.formFilename(format("%s_%02d", UNSORTED_BAM_ID, threadId), BAM_EXTENSION);

            SAMFileHeader fileHeader = mSamReader.getFileHeader().clone();

            // the writer resolves each record's contigs by name against the ref genome's sequence dictionary
            if(mConvertRefGenome)
                fileHeader.setSequenceDictionary(config.RefDictionary);

            fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

            mBamWriter = new SAMFileWriterFactory().makeBAMWriter(fileHeader, false, new File(mBamFilename));
        }
        else
        {
            mBamWriter = null;
            mBamFilename = "";
        }

        mPartitionChecker = new PartitionChecker(
                config, fragmentCache, mSamReader, mBamWriter, mAlignmentTagCalculator);
    }

    public static List<ChrBaseRegion> splitRegionsIntoPartitions(final CheckConfig config)
    {
        return PartitionTask.splitRegionsIntoPartitions(
                config.RefGenomeFile, config.Threads, config.PartitionSize, config.BamFile, config.SpecificChrRegions);
    }

    public FragmentStats stats() { return mPartitionChecker.stats(); }
    public String bamFilename() { return mBamFilename; }
    public AlignmentTagCalculator alignmentTagCalculator() { return mAlignmentTagCalculator; }
    public long nonRefContigDropped() { return mNonRefContigDropped + mPartitionChecker.nonRefContigDropped(); }

    @Override
    public void run()
    {
        while(true)
        {
            try
            {
                ChrBaseRegion region = (ChrBaseRegion)mPartitions.removeItem();

                mPartitionChecker.processPartition(region);
            }
            catch(NoSuchElementException e)
            {
                BT_LOGGER.trace("all tasks complete");
                break;
            }
            catch(Exception e)
            {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }

    public void writeUnmappedReads()
    {
        SAMRecordIterator iterator = mSamReader.queryUnmapped();

        while(iterator.hasNext())
        {
            SAMRecord record = iterator.next();

            if(mReverseBqsr)
                PartitionChecker.reverseBqsrQuals(record);

            writeRecord(record);
        }
    }

    public void writeIncompleteReads(final List<SAMRecord> reads)
    {
        reads.forEach(x -> writeRecord(x));
    }

    private void writeRecord(final SAMRecord record)
    {
        if(mConvertRefGenome && !PartitionChecker.hasRefGenomeContigs(record))
        {
            ++mNonRefContigDropped;
            return;
        }

        if(mAlignmentTagCalculator != null)
            mAlignmentTagCalculator.recalculate(record);

        mBamWriter.addAlignment(record);
    }

    public void close()
    {
        if(mBamWriter != null)
            mBamWriter.close();
    }
}
