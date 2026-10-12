package com.hartwig.hmftools.bamtools.checker;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

import org.junit.Test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class ConvertRefGenomeTest
{
    private static final String VIRAL_CONTIG = "HPV16";

    private final SAMFileHeader mInputHeader;
    private final SAMSequenceDictionary mRefDictionary;

    public ConvertRefGenomeTest()
    {
        // input BAM: chr1, chr2, chrEBV, HPV16; ref genome: chr1, chr2, chr1_alt, chrEBV, so chrEBV changes index from 2 to 3
        mInputHeader = new SAMFileHeader();
        mInputHeader.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 10000), new SAMSequenceRecord("chr2", 10000),
                new SAMSequenceRecord("chrEBV", 5000), new SAMSequenceRecord(VIRAL_CONTIG, 8000))));

        mRefDictionary = new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 10000), new SAMSequenceRecord("chr2", 10000),
                new SAMSequenceRecord("chr1_alt", 2000), new SAMSequenceRecord("chrEBV", 5000)));

        mRefDictionary.getSequences().forEach(x -> CheckConfig.Params.ValidContigs.add(x.getSequenceName()));
    }

    private SAMRecord createRead(final String contig, int position, final String mateContig, int matePosition)
    {
        SAMRecord read = new SAMRecord(mInputHeader);
        read.setReadName("READ_01");
        read.setReadPairedFlag(true);
        read.setFirstOfPairFlag(true);
        read.setReferenceName(contig);
        read.setAlignmentStart(position);
        read.setCigarString("20M");
        read.setReadBases("ACGTACGTACGTACGTACGT".getBytes());
        read.setBaseQualities(new byte[20]);
        read.setMateReferenceName(mateContig);
        read.setMateAlignmentStart(matePosition);
        return read;
    }

    @Test
    public void testReadsWithNonRefContigs()
    {
        assertTrue(PartitionChecker.hasRefGenomeContigs(createRead("chrEBV", 100, "chr2", 500)));
        assertFalse(PartitionChecker.hasRefGenomeContigs(createRead(VIRAL_CONTIG, 1, "chr1", 500)));
        assertFalse(PartitionChecker.hasRefGenomeContigs(createRead("chr1", 500, VIRAL_CONTIG, 1)));

        SAMRecord unpairedRead = createRead("chr1", 500, VIRAL_CONTIG, 1);
        unpairedRead.setReadPairedFlag(false);
        assertTrue(PartitionChecker.hasRefGenomeContigs(unpairedRead));
    }

    @Test
    public void testWriterResolvesContigsByName() throws IOException
    {
        SAMRecord read = createRead("chrEBV", 100, "chr2", 500);
        assertEquals(2, read.getReferenceIndex().intValue());

        File bamFile = Files.createTempFile("convert_ref_genome", ".bam").toFile();
        bamFile.deleteOnExit();

        // as for the partition BAM writers
        SAMFileHeader outputHeader = mInputHeader.clone();
        outputHeader.setSequenceDictionary(mRefDictionary);
        outputHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        try(SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(outputHeader, false, bamFile))
        {
            writer.addAlignment(read);
        }

        try(SamReader reader = SamReaderFactory.makeDefault().open(bamFile))
        {
            assertTrue(reader.getFileHeader().getSequenceDictionary().isSameDictionary(mRefDictionary));

            List<SAMRecord> written = Lists.newArrayList(reader.iterator());
            assertEquals(1, written.size());
            assertEquals("chrEBV", written.get(0).getReferenceName());
            assertEquals(3, written.get(0).getReferenceIndex().intValue());
            assertEquals(1, written.get(0).getMateReferenceIndex().intValue());
        }
    }

    @Test
    public void testContigLengthCheck()
    {
        assertTrue(BamChecker.checkContigLengths(mInputHeader, mRefDictionary));

        SAMFileHeader mismatchHeader = new SAMFileHeader();
        mismatchHeader.setSequenceDictionary(new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", 9999))));
        assertFalse(BamChecker.checkContigLengths(mismatchHeader, mRefDictionary));
    }

    @Test
    public void testRefDictionaryFallsBackToFastaIndex() throws IOException
    {
        File refDir = Files.createTempDirectory("ref_genome").toFile();
        refDir.deleteOnExit();

        File fastaFile = new File(refDir, "ref.fasta");
        Files.writeString(fastaFile.toPath(), ">chr1\nACGTACGT\n>chr2\nACGT\n");
        File indexFile = new File(refDir, "ref.fasta.fai");
        Files.writeString(indexFile.toPath(), "chr1\t8\t6\t8\t9\nchr2\t4\t21\t4\t5\n");
        fastaFile.deleteOnExit();
        indexFile.deleteOnExit();

        RefGenomeSource refGenome = RefGenomeSource.loadRefGenome(fastaFile.getAbsolutePath());
        assertNull(refGenome.refGenomeFile().getSequenceDictionary());

        SAMSequenceDictionary dictionary = CheckConfig.loadRefDictionary(
                fastaFile.getAbsolutePath(), refGenome.refGenomeFile().getSequenceDictionary());

        assertEquals(2, dictionary.size());
        assertEquals("chr1", dictionary.getSequence(0).getSequenceName());
        assertEquals(8, dictionary.getSequence(0).getSequenceLength());
        assertEquals("chr2", dictionary.getSequence(1).getSequenceName());
        assertEquals(4, dictionary.getSequence(1).getSequenceLength());

        // a dictionary when present is used as is
        assertEquals(mRefDictionary, CheckConfig.loadRefDictionary(fastaFile.getAbsolutePath(), mRefDictionary));
    }
}
