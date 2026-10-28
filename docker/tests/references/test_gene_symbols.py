"""buildDB/make_gene_symbols.py (release 0.3.0, 3f) on small made-up inputs in the formats of the real files
(HGNC, MGI and MirGeneDB checked with --inspect on 2026-10-26; miRBase still assumed). Run from the repository root:
python -m pytest docker/tests/references"""
import gzip
import os
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / 'buildDB'))
import make_gene_symbols as mgs  # noqa: E402

BED = ('ENST01|BTG3\t10\t500\tBTG3\tprotein_coding\t+\n'
       'ENST02.3|ENSG00000252404\t10\t150\tENSG00000252404\tscaRNA\t+\n'
       '21\t100\t900\tBTG3_intron\tpre-mRNA\t-\n'
       '1\t100\t900\tENSG00000252404_intron\tpre-mRNA\t-\n'
       'Hsa-Mir-21_5p\t10\t32\tHsa-Mir-21_5p\tmiRNA\t+\n'
       'Hsa-Mir-21_3p*\t10\t32\tHsa-Mir-21_3p*\tmiRNA-star\t+\n'
       'Hsa-Mir-21_pre\t10\t70\tHsa-Mir-21_pre\tmiRNA-pre\t+\n'
       'Hsa-Let-7-P1d_5p|Hsa-Let-7-P2a1_5p\t10\t32\tHsa-Let-7-P1d_5p|Hsa-Let-7-P2a1_5p\tmiRNA\t+\n'
       'Homo_sapiens_tRNA-Ala-AGC-1\t10\t85\tHomo_sapiens_tRNA-Ala-AGC-1\ttRNA\t+\n'
       '1\t5\t80\ttRNA-Ala-AGC-1-2\tpre-tRNA\t-\n'
       'piR-hsa-9\t10\t40\tpiR-hsa-9\tpiRNA\t+\n')
GTF = ('#!genome-build GRCh38\n'
       '21\tensembl\tgene\t1\t2\t.\t+\t.\tgene_id "ENSG01"; gene_name "BTG3"; gene_biotype "protein_coding";\n'
       '21\tensembl\ttranscript\t1\t2\t.\t+\t.\tgene_id "ENSG01"; transcript_id "ENST01"; gene_name "BTG3";\n'
       '1\tensembl\tgene\t1\t2\t.\t+\t.\tgene_id "ENSG00000252404"; gene_biotype "scaRNA";\n'
       '1\tensembl\ttranscript\t1\t2\t.\t+\t.\tgene_id "ENSG00000252404"; transcript_id "ENST02";\n')
HGNC = ('hgnc_id\tsymbol\tname\tensembl_gene_id\tmirbase\tgtrnadb\n'
        'HGNC:1132\tBTG3\tBTG anti-proliferation factor 3\tENSG01\t\t\n'
        'HGNC:9\tSCARNA99\tsmall Cajal body-specific RNA\tENSG00000252404\t\t\n'
        'HGNC:31586\tMIR21\tmicroRNA 21\t\tMI0000077\t\n'
        'HGNC:31479\tMIRLET7D\tmicroRNA let-7d\t\tMI0000065\t\n'
        'HGNC:31477\tMIRLET7A1\tmicroRNA let-7a-1\t\tMI0000060\t\n'
        'HGNC:34500\tTRA-AGC1-1\ttRNA-Ala (anticodon AGC) 1-1\t\t\ttRNA-Ala-AGC-1-1\n'
        'HGNC:34501\tTRA-AGC1-2\ttRNA-Ala (anticodon AGC) 1-2\t\t\ttRNA-Ala-AGC-1-2\n')
GFF = ('##gff-version 3\n'
       'chr17\t.\tpre_miRNA\t1\t70\t.\t+\t.\tID=Hsa-Mir-21_pre;Alias=MI0000077\n'
       'chr9\t.\tpre_miRNA\t1\t70\t.\t+\t.\tID=Hsa-Let-7-P1d_pre;Alias=MI0000065\n'
       'chr9\t.\tmiRNA\t1\t22\t.\t+\t.\tID=Hsa-Let-7-P1d_5p;Alias=MIMAT0000065\n'
       'chr9\t.\tpre_miRNA\t1\t70\t.\t+\t.\tID=Hsa-Let-7-P2a1_pre;Alias=MI0000060\n'
       'chr9\t.\tmiRNA\t1\t22\t.\t+\t.\tID=Hsa-Let-7-P2a1_5p;Alias=MIMAT0000062\n'
       'chr17\t.\tmiRNA\t1\t22\t.\t+\t.\tID=Hsa-Mir-21_5p;Alias=MIMAT0000076\n'
       'chr17\t.\tmiRNA\t40\t62\t.\t+\t.\tID=Hsa-Mir-21_3p*;Alias=MIMAT0004494\n')
MATURE = ('>hsa-miR-21-5p MIMAT0000076 Homo sapiens miR-21-5p\nUAGC\n>hsa-miR-21-3p MIMAT0004494 Homo sapiens miR-21-3p\nCAAC\n'
          '>hsa-let-7d-5p MIMAT0000065 Homo sapiens let-7d-5p\nAGAG\n>hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p\nUGAG\n')
HAIRPIN = '>hsa-mir-21 MI0000077 Homo sapiens miR-21 stem-loop\nUGUC\n'


def write(tmp_path, name, text, gz=False):
    path = tmp_path / name
    if gz:
        with gzip.open(path, 'wt') as f:
            f.write(text)
    else:
        path.write_text(text)
    return str(path)


def test_every_join_route(tmp_path):
    transcript_gene, by_name = mgs.read_gtf(write(tmp_path, 'g.gtf.gz', GTF, gz=True))
    rows = mgs.symbol_rows('human', mgs.read_bed(write(tmp_path, 'human.bed', BED)), transcript_gene, by_name,
                           mgs.read_hgnc(write(tmp_path, 'hgnc.txt', HGNC)),
                           mgs.read_mirgenedb_gff(write(tmp_path, 'm.gff', GFF)),
                           mgs.read_mirbase_fasta(write(tmp_path, 'mature.fa', MATURE)),
                           mgs.read_mirbase_fasta(write(tmp_path, 'hairpin.fa', HAIRPIN)))
    assert rows == [
        ('BTG3', 'HGNC:1132', 'BTG3', ''),                                  # Ensembl transcript -> gene -> HGNC
        ('ENSG00000252404', 'HGNC:9', 'SCARNA99', ''),                      # unnamed gene; versioned transcript id
        ('BTG3_intron', 'HGNC:1132', 'BTG3', ''),                           # pre-mRNA by name, protein_coding gene
        ('ENSG00000252404_intron', 'HGNC:9', 'SCARNA99', ''),               # pre-mRNA of an unnamed gene: by its id
        ('Hsa-Mir-21_5p', 'HGNC:31586', 'MIR21', 'hsa-miR-21-5p'),          # arm: precursor MI -> HGNC; MIMAT name
        ('Hsa-Mir-21_3p*', 'HGNC:31586', 'MIR21', 'hsa-miR-21-3p'),         # star arm
        ('Hsa-Mir-21_pre', 'HGNC:31586', 'MIR21', 'hsa-mir-21'),            # precursor: hairpin name
        ('Hsa-Let-7-P1d_5p|Hsa-Let-7-P2a1_5p', 'HGNC:31479;HGNC:31477', 'MIRLET7D;MIRLET7A1',
         'hsa-let-7d-5p;hsa-let-7a-5p'),                                    # merged arms: each part
        ('Homo_sapiens_tRNA-Ala-AGC-1', 'HGNC:34500;HGNC:34501', 'TRA-AGC1-1;TRA-AGC1-2', ''),   # family: all loci
        ('tRNA-Ala-AGC-1-2', 'HGNC:34501', 'TRA-AGC1-2', ''),               # one locus
    ]                                                                       # piRNA: no row


def test_mgi_report_columns(tmp_path):
    report = write(tmp_path, 'MRK_ENSEMBL.rpt',
                   'MGI:88180\tBtg3\tB cell translocation gene 3\t40.9\t16\tENSMUSG00000022863\tENSMUST01\n'
                   'not a row\n')
    index = mgs.read_mgi(report, mgs.SOURCES['mouse']['mgi_columns'])
    assert index['ensembl'] == {'ENSMUSG00000022863': [('MGI:88180', 'Btg3')]}


def test_wrong_hgnc_columns_stop_the_build(tmp_path):
    with pytest.raises(SystemExit, match='columns not found'):
        mgs.read_hgnc(write(tmp_path, 'hgnc.txt', 'hgnc_id\tsymbol\n'))


def test_placeholder_sources_are_refused():
    placeholder = dict(mgs.SOURCES['human'], status='PLACEHOLDER: not checked')
    with pytest.raises(SystemExit, match='placeholders'):
        mgs.build('human', placeholder, {}, '.')
    assert all(spec['status'].startswith('verified') for spec in mgs.SOURCES.values())   # v1 sources, 2026-10-27


def test_downloads_go_to_a_folder_per_species(tmp_path, monkeypatch):
    """The downloads are named by source (gtf.gz, mirgenedb_gff, ...): one subfolder per species, so that a shared
    --work folder cannot give mouse the human files."""
    seen = []
    monkeypatch.setattr(mgs, 'fetch', lambda spec, work, skip_failed=False: seen.append(work) or {})
    monkeypatch.setattr(mgs, 'inspect', lambda files: None)
    for species in ('human', 'mouse'):
        mgs.main(['--species', species, '--work', str(tmp_path), '--inspect'])
    assert seen == [str(tmp_path / 'human'), str(tmp_path / 'mouse')]


def test_provenance_paths_hold_no_user_folders(tmp_path):
    """The references BED is recorded relative to the repository, or by its file name only (the release leak scan
    refuses user and computer names, 2026-10-27)."""
    assert mgs.shown_path(os.path.join(mgs.REPO, 'docker', 'assets', 'human.bed')) == 'docker/assets/human.bed'
    assert mgs.shown_path(str(tmp_path / 'human.bed')) == 'human.bed'
