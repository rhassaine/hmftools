package com.hartwig.hmftools.bamtools.checker;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MISMATCHES_AND_DELETIONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import java.util.List;

import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.SequenceUtil;

public class AlignmentTagCalculatorTest
{
    private static final String CHR = "chr1";

    // positions 1-40 are real sequence, 41-60 are masked with N, 61-80 are real sequence
    private static final String REF_BASES = "ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA" + "NNNNNNNNNNNNNNNNNNNN" + "GGCCTTAAGGCCTTAAGGCC";

    private final MockRefGenome mRefGenome;
    private final SAMFileHeader mHeader;

    public AlignmentTagCalculatorTest()
    {
        mRefGenome = new MockRefGenome(true);
        mRefGenome.RefGenomeMap.put(CHR, REF_BASES);
        mRefGenome.ChromosomeLengths.put(CHR, REF_BASES.length());

        mHeader = new SAMFileHeader();
        mHeader.setSequenceDictionary(new SAMSequenceDictionary(List.of(new SAMSequenceRecord(CHR, REF_BASES.length()))));
    }

    private SAMRecord createRead(int position, final String readBases, final String cigar)
    {
        SAMRecord read = new SAMRecord(mHeader);
        read.setReadName("READ_01");
        read.setReferenceName(CHR);
        read.setAlignmentStart(position);
        read.setCigarString(cigar);
        read.setReadBases(readBases.getBytes());
        read.setBaseQualities(new byte[readBases.length()]);
        read.setMappingQuality(60);
        read.setAttribute(NUM_MUTATONS_ATTRIBUTE, 99);
        read.setAttribute(MISMATCHES_AND_DELETIONS_ATTRIBUTE, "stale");
        return read;
    }

    private static String refBases(int posStart, int posEnd) { return REF_BASES.substring(posStart - 1, posEnd); }

    @Test
    public void testMatchAndMismatches()
    {
        AlignmentTagCalculator calculator = new AlignmentTagCalculator(mRefGenome);

        SAMRecord read = createRead(1, refBases(1, 20), "20M");
        calculator.recalculate(read);
        assertEquals(0, read.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE).intValue());
        assertEquals("20", read.getStringAttribute(MISMATCHES_AND_DELETIONS_ATTRIBUTE));

        // mismatches at read index 5 (ref G) and the last base (ref T)
        String bases = refBases(1, 20);
        bases = bases.substring(0, 5) + "A" + bases.substring(6, 19) + "A";
        read = createRead(1, bases, "20M");
        calculator.recalculate(read);
        assertEquals(2, read.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE).intValue());
        assertEquals("5G13T0", read.getStringAttribute(MISMATCHES_AND_DELETIONS_ATTRIBUTE));
    }

    @Test
    public void testIndelsAndClipping()
    {
        AlignmentTagCalculator calculator = new AlignmentTagCalculator(mRefGenome);

        // 3S soft-clip, 5M at 11-15, 2I, 4M at 16-19, 3D at 20-22, 6M at 23-28
        String bases = "TTT" + refBases(11, 15) + "GG" + refBases(16, 19) + refBases(23, 28);
        SAMRecord read = createRead(11, bases, "3S5M2I4M3D6M");
        calculator.recalculate(read);

        assertEquals(5, read.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE).intValue());
        assertEquals("9^" + refBases(20, 22) + "6", read.getStringAttribute(MISMATCHES_AND_DELETIONS_ATTRIBUTE));
    }

    @Test
    public void testMaskedReferenceBases()
    {
        AlignmentTagCalculator calculator = new AlignmentTagCalculator(mRefGenome);

        // a read aligned across the start of the N-masked region: 5 real bases then 15 against N
        String bases = refBases(36, 40) + "ACGTACGTACGTACG";
        SAMRecord read = createRead(36, bases, "20M");
        calculator.recalculate(read);

        assertEquals(15, read.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE).intValue());
        assertEquals("5N0N0N0N0N0N0N0N0N0N0N0N0N0N0N0", read.getStringAttribute(MISMATCHES_AND_DELETIONS_ATTRIBUTE));
        assertFalse(SequenceUtil.basesEqual((byte)'A', (byte)'N'));
    }

    @Test
    public void testMatchesHtsjdkCalculation()
    {
        // htsjdk's version needs the whole contig's bases, this one works from a window, so they should agree
        AlignmentTagCalculator calculator = new AlignmentTagCalculator(mRefGenome, 12);

        String[][] reads = {
                { "3", "TTGAAACGTTGCAACG", "16M" },
                { "11", "TTT" + refBases(11, 15) + "GG" + refBases(16, 19) + refBases(23, 28), "3S5M2I4M3D6M" },
                { "30", refBases(30, 40) + "ACGTA", "16M" },
                { "58", "AAC" + refBases(61, 70), "13M" },
                { "62", "AA" + refBases(62, 63) + refBases(66, 75), "2S2M2D10M" } };

        for(String[] readData : reads)
        {
            SAMRecord read = createRead(Integer.parseInt(readData[0]), readData[1], readData[2]);
            SAMRecord htsjdkRead = createRead(Integer.parseInt(readData[0]), readData[1], readData[2]);

            calculator.recalculate(read);
            SequenceUtil.calculateMdAndNmTags(htsjdkRead, REF_BASES.getBytes(), true, true);

            assertEquals(htsjdkRead.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE), read.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE));
            assertEquals(htsjdkRead.getStringAttribute(MISMATCHES_AND_DELETIONS_ATTRIBUTE), read.getStringAttribute(MISMATCHES_AND_DELETIONS_ATTRIBUTE));
        }
    }

    @Test
    public void testUnrecalculableReadsUnchanged()
    {
        AlignmentTagCalculator calculator = new AlignmentTagCalculator(mRefGenome);

        // alignment beyond the contig end
        SAMRecord read = createRead(75, "ACGTACGTAC", "10M");
        calculator.recalculate(read);
        assertEquals(99, read.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE).intValue());

        // fewer read bases than the CIGAR
        read = createRead(1, refBases(1, 15), "20M");
        calculator.recalculate(read);
        assertEquals(99, read.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE).intValue());
    }
}
