package com.hartwig.hmftools.bamtools.checker;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;

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

public class NonRefContigFilterTest
{
    private static final String VIRAL_CONTIG = "HPV16";
    private static final String ALT_CONTIG = "chr1_alt";

    private final SAMFileHeader mInputHeader;
    private final SAMSequenceDictionary mRefDictionary;
    private final NonRefContigFilter mFilter;

    public NonRefContigFilterTest()
    {
        // input BAM: chr1, chr2, chrEBV, HPV16; ref genome: chr1, chr2, chr1_alt, chrEBV, so chrEBV changes index from 2 to 3
        SAMSequenceDictionary inputDictionary = new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 10000), new SAMSequenceRecord("chr2", 10000),
                new SAMSequenceRecord("chrEBV", 5000), new SAMSequenceRecord(VIRAL_CONTIG, 8000)));

        mInputHeader = new SAMFileHeader();
        mInputHeader.setSequenceDictionary(inputDictionary);

        mRefDictionary = new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 10000), new SAMSequenceRecord("chr2", 10000),
                new SAMSequenceRecord(ALT_CONTIG, 2000), new SAMSequenceRecord("chrEBV", 5000)));

        mFilter = new NonRefContigFilter(mRefDictionary);
    }

    private SAMRecord createRead(final String readId, final String contig, int position, final String mateContig, int matePosition)
    {
        SAMRecord read = new SAMRecord(mInputHeader);
        read.setReadName(readId);
        read.setReadPairedFlag(true);
        read.setFirstOfPairFlag(true);
        read.setReferenceName(contig);
        read.setAlignmentStart(position);
        read.setCigarString("20M");
        read.setMappingQuality(60);
        read.setReadBases("ACGTACGTACGTACGTACGT".getBytes());
        read.setBaseQualities(new byte[20]);
        read.setMateReferenceName(mateContig);
        read.setMateAlignmentStart(matePosition);
        return read;
    }

    @Test
    public void testRefContigsUnchanged()
    {
        SAMRecord read = createRead("READ_01", "chrEBV", 100, "chr2", 500);
        read.setProperPairFlag(true);
        read.setAttribute(SUPPLEMENTARY_ATTRIBUTE, "chr1,600,+,10S10M,60,0;");

        assertTrue(mFilter.apply(read));
        assertEquals("chrEBV", read.getReferenceName());
        assertEquals("chr2", read.getMateReferenceName());
        assertTrue(read.getProperPairFlag());
        assertEquals("chr1,600,+,10S10M,60,0;", read.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE));
    }

    @Test
    public void testViralPrimaryUnmappedAtMate()
    {
        SAMRecord read = createRead("READ_01", VIRAL_CONTIG, 1, "chr1", 500);
        read.setAttribute(SUPPLEMENTARY_ATTRIBUTE, "chr1,600,+,10S10M,60,0;");
        read.setAttribute(NUM_MUTATONS_ATTRIBUTE, 0);
        read.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, 20);

        assertTrue(mFilter.apply(read));
        assertTrue(read.getReadUnmappedFlag());
        assertEquals("chr1", read.getReferenceName());
        assertEquals(500, read.getAlignmentStart());
        assertEquals(0, read.getMappingQuality());
        assertEquals("*", read.getCigarString());
        assertNull(read.getAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertNull(read.getAttribute(NUM_MUTATONS_ATTRIBUTE));
        assertNull(read.getAttribute(ALIGNMENT_SCORE_ATTRIBUTE));
        assertEquals("ACGTACGTACGTACGTACGT", read.getReadString());
    }

    @Test
    public void testViralPairFullyUnmapped()
    {
        SAMRecord read = createRead("READ_01", VIRAL_CONTIG, 1, VIRAL_CONTIG, 300);

        assertTrue(mFilter.apply(read));
        assertTrue(read.getReadUnmappedFlag());
        assertTrue(read.getMateUnmappedFlag());
        assertEquals(NO_CHROMOSOME_NAME, read.getReferenceName());
        assertEquals(NO_CHROMOSOME_NAME, read.getMateReferenceName());
    }

    @Test
    public void testViralSupplementaryAndSecondaryDropped()
    {
        SAMRecord supp = createRead("READ_01", VIRAL_CONTIG, 1, "chr1", 500);
        supp.setSupplementaryAlignmentFlag(true);
        assertFalse(mFilter.apply(supp));

        SAMRecord secondary = createRead("READ_02", VIRAL_CONTIG, 1, "chr1", 500);
        secondary.setSecondaryAlignment(true);
        assertFalse(mFilter.apply(secondary));
    }

    @Test
    public void testHumanReadWithViralMate()
    {
        SAMRecord read = createRead("READ_01", "chr1", 500, VIRAL_CONTIG, 1);
        read.setProperPairFlag(true);
        read.setAttribute(MATE_CIGAR_ATTRIBUTE, "20M");

        assertTrue(mFilter.apply(read));
        assertFalse(read.getReadUnmappedFlag());
        assertTrue(read.getMateUnmappedFlag());
        assertFalse(read.getProperPairFlag());
        assertEquals("chr1", read.getMateReferenceName());
        assertEquals(500, read.getMateAlignmentStart());
        assertNull(read.getAttribute(MATE_CIGAR_ATTRIBUTE));
    }

    @Test
    public void testSupplementaryWithViralMateAndSuppData()
    {
        SAMRecord read = createRead("READ_01", "chr2", 2000, VIRAL_CONTIG, 1);
        read.setSupplementaryAlignmentFlag(true);
        read.setAttribute(SUPPLEMENTARY_ATTRIBUTE, "chr1,500,+,20M,60,0;" + VIRAL_CONTIG + ",50,+,5S15M,60,0;");

        assertTrue(mFilter.apply(read));
        assertTrue(read.getMateUnmappedFlag());
        assertEquals("chr1", read.getMateReferenceName());
        assertEquals(500, read.getMateAlignmentStart());
        assertEquals("chr1,500,+,20M,60,0", read.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE));
    }

    @Test
    public void testWriterResolvesContigsByName() throws IOException
    {
        // the BAM writer with the ref genome dictionary re-indexes records by contig name
        List<SAMRecord> reads = Lists.newArrayList(
                createRead("READ_01", "chrEBV", 100, "chr2", 500),
                createRead("READ_02", VIRAL_CONTIG, 1, "chr1", 300),
                createRead("READ_03", "chr1", 300, VIRAL_CONTIG, 1));

        assertEquals(2, reads.get(0).getReferenceIndex().intValue());

        File bamFile = Files.createTempFile("non_ref_contig_filter", ".bam").toFile();
        bamFile.deleteOnExit();

        SAMFileHeader outputHeader = NonRefContigFilter.formOutputHeader(mInputHeader, mRefDictionary);
        outputHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        try(SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(outputHeader, false, bamFile))
        {
            for(SAMRecord read : reads)
            {
                if(mFilter.apply(read))
                    writer.addAlignment(read);
            }
        }

        try(SamReader reader = SamReaderFactory.makeDefault().open(bamFile))
        {
            assertTrue(reader.getFileHeader().getSequenceDictionary().isSameDictionary(mRefDictionary));

            List<SAMRecord> written = Lists.newArrayList(reader.iterator());
            assertEquals(3, written.size());

            assertEquals("chrEBV", written.get(0).getReferenceName());
            assertEquals(3, written.get(0).getReferenceIndex().intValue());
            assertEquals(1, written.get(0).getMateReferenceIndex().intValue());

            assertTrue(written.get(1).getReadUnmappedFlag());
            assertEquals("chr1", written.get(1).getReferenceName());

            assertTrue(written.get(2).getMateUnmappedFlag());
            assertEquals("chr1", written.get(2).getMateReferenceName());
        }
    }

    @Test
    public void testContigLengthCheck()
    {
        assertTrue(NonRefContigFilter.checkContigLengths(mInputHeader, mRefDictionary));

        SAMFileHeader mismatchHeader = new SAMFileHeader();
        mismatchHeader.setSequenceDictionary(new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", 9999))));
        assertFalse(NonRefContigFilter.checkContigLengths(mismatchHeader, mRefDictionary));
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
