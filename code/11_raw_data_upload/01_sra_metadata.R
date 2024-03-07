library(tidyverse)
library(here)
library(sessioninfo)

sample_info_path = here('processed-data', '01_SPEAQeasy', 'data_info.csv')
design_description_text = paste(
    "Assays were completed using 22 individual blocks of postmortem human DLPFC tissue collected across anterior (Ant), middle (Mid), and posterior (Post) positions. These were a subset of the same tissue blocks used for Visium (10x Genomics) and 3’ gene expression snRNA-seq assays described in Huuki-Myers et al. (doi: 10.1101/2023.02.15.528722). Tissue sections for the majority of assays were cryosectioned on the same day for each run (i.e. a round of 3-4 tissue blocks balanced across anterior-posterior DLPFC axis) to minimize tissue loss that occurs when obtaining a flat cutting face on the tissue block. In a Leica 3050 cryostat, blocks were allowed to equilibrate for 30 minutes prior to mounting with optimal cutting temperature (OCT) medium. Excess OCT was scored from each side of the block using a chilled razor blade to minimize interference of OCT in RNA extraction. Each block was trimmed to achieve a flat cutting face and then several 10um serial sections were collected for RNAScope/immunofluorescence (IF) assays across the 3-4 tissue blocks in that round. In particular, eight slides containing 4 tissue sections each were collected per round, with one tissue section from each block on a given slide. Following collection of serial sections for RNAScope/IF, the cutting thickness was adjusted to 100um, and ten serial sections (~1mm of tissue) were collected for single nucleus RNA sequencing (snRNA-seq). Processing and analysis of snRNA-seq data was reported in Huuki-Myers et al. (doi: 10.1101/2023.02.15.528722). Immediately following tissue collection for snRNA-seq, six 100um serial sections (~600um tissue) were collected for bulk RNA extraction. Six additional 100um serial sections were collected for fractionated (nuclear/cytosolic) RNA extraction. Collected tissue was stored in Eppendorf tubes at -80°C until use.",
    "Total RNA was extracted from tissue aliquots (2 extractions per block for 19 tissue blocks, n=38) using the Qiagen RNeasy mini kit (RNeasy Mini Kit, Cat No. 74104, Qiagen, Hilden, Germany) . A modified version of Qiagen’s “Purification of Total RNA from Animal Tissues” protocol from the 2020 version of the RNeasy Mini Handbook was used. Briefly, tissue cryosections were homogenized via wide bore pipette in 0.7 mL of Trizol. Next, 0.14 mL of chloroform was added, and the aqueous phase of the gradient was removed and transferred to a new tube. An equal volume of 70% ethanol was added, and then the mixture was put onto an RNeasy mini column. At this point, RNA was extracted according to the manufacturer's instructions with DNAse digestion treatment. RNA quantity was measured using a Qubit 4 fluorometer (Qubit 4 fluorometer; Qubit dnDNA HS Assay Kit, Cat No. Q32854 Invitrogen, Eugene, OR, United States). RNA quality was assessed using an Agilent RNA Nano kit on a BioAnalyzer instrument (RNA 6000 Nano Kit, Agilent, Santa Clara, CA, United States). Libraries were subsequently prepared and sequenced at Psomagen. For each sample, 100-500ng of RNA from the same tube was used to prepare a “RiboZeroGold” library with the TruSeq Stranded Total RNA with Ribo-Zero Gold Library Prep kit (Illumina) and “PolyA” library with TruSeq Stranded mRNA Library Prep Kit (Illumina) according to manufacturer’s instructions. Libraries were sequenced on an Illumina Novaseq 6000 targeting 80 million reads per sample. ERCC spike in sequences were included in all samples except the initial pilot round (n = 24).",
    "Fractionated RNA extraction was performed on tissue aliquots using the Cytoplasmic and Nuclear RNA Purification kit (Norgen Biotek, Cat. No. 21000, ON, Canada) (cyto n=38, nuc n=37) according to the “Animal Tissues” protocol in the manufacturer’s manual (PI21000-19, section B). Briefly, reagent J was added to the tissue, which was homogenized via a wide bore pipette. Lysate was spun resulting in a supernatant and a pellet. The supernatant was removed and used for the cytoplasmic fraction, and the pellet was retained for the nuclear fraction. For cytoplasmic RNA purification, buffer SK and 100% ethanol were added to the supernatant, which was then transferred to a spin column, and centrifuged. Flow through was discarded and on column DNA removal was completed (RNase-Free DNase I Kit, Cat No. 25710, Norgen Biotek, ON, Canada) For nuclear RNA purification, the pellet was resuspended in buffer SK and 100% ethanol was added. Lysate was then passed through a 25 gauge needle 5 times and added to a different spin column. Following centrifugation,the flow through was discarded. For both fractions, the spin columns were washed twice with Wash Solution A before drying the membranes. Cytoplasmic and nuclear RNA were eluted from each column in Elution Buffer E. Both fractions of RNA were stored at -80 degrees Celsius. RNA quality and quantity was measured as described above for bulk RNA extraction. RiboZeroGold and PolyA libraries were subsequently prepared and sequenced at Psomagen as described above. In total n=113 bulk RNA-seq samples were generated (19 * 2 * 3 = 114 minus one) as a sample failed during library preparation due to an insufficient amount of starting material)."
)
required_cols = c(
    'sample_name', 'library_ID', 'title', 'library_strategy', 'library_source',
    'library_selection', 'library_layout', 'platform', 'instrument_model',
    'design_description', 'filetype', 'filename', 'filename2'
)
out_path = here('processed-data', '11_raw_data_upload', 'metadata.tsv')

meta_df = sample_info_path |>
    read.csv() |>
    as_tibble() |>
    rename(
        sample_name = SAMPLE_ID,
        filename = fastq1,
        filename2 = fastq2
    ) |>
    mutate(
        library_ID = sample_name,
        title = paste(
            library_prep, 'RNA-seq from the', location, 'human DLPFC:', Dataset,
            'dataset'
        ),
        library_strategy = 'RNA-Seq',
        library_source = 'TRANSCRIPTOMIC',
        library_selection = case_when(
            library_type == 'polyA' ~ 'PolyA',
            library_type == 'RiboZeroGold' ~ 'Inverse rRNA',
            TRUE ~ NA
        ),
        library_layout = 'paired',
        platform = 'ILLUMINA',
        instrument_model = 'Illumina NovaSeq 6000',
        design_description = factor(design_description_text),
        filetype = 'fastq'
    )

#   Re-order columns to match SRA's expectations and write to TSV
stopifnot(all(required_cols %in% colnames(meta_df)))
meta_df |>
    select(all_of(required_cols)) |>
    write_tsv(out_path)


#   Also check disk usage (to make sure we comply with SRA requirements)
disk_usage = sapply(
    c(meta_df$filename, meta_df$filename2),
    function(x) {
        as.integer(system(sprintf('du -k %s | cut -f 1', x), intern = TRUE))
    }
)
message(
    sprintf(
        "FASTQ files occupy a total of %sGB.", round(sum(disk_usage) / 1e6, 1)
    )
)

#   Individual files must be less than 100GB
stopifnot(all(disk_usage < 100e6))

session_info()
