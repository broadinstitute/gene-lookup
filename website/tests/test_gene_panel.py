"""Tests for gene-panel / gene-list search on the index page.

Covers the two ways to search many genes at once:
  - typing a comma-separated list into the search box
  - uploading a file via the (+) button: a TXT/TSV gene list, a BED region file,
    or any of these gzip-compressed

and the "found vs. missing" summary banner (including the "more..." modal).

These tests exercise the live BigQuery backend, mirroring test_basic_search.py.
"""

import gzip

import pytest


# Environmental console errors to ignore: Google Analytics is often unreachable in CI/sandbox
# environments, which surfaces as "Failed to load resource" for google-analytics.com.
def real_errors(console_errors):
    """Return console errors excluding blocked-analytics network failures."""
    return [
        e for e in console_errors
        if "failed to load resource" not in e.lower() and "google-analytics" not in e.lower()
    ]


def summary_text(page):
    """Return the visible gene-panel summary banner text, waiting for it to appear."""
    page.wait_for_selector("#gene-panel-summary", state="visible", timeout=30000)
    return page.inner_text("#gene-panel-summary").replace("\n", " ")


def wait_summary_contains(page, needle):
    """Wait until the summary banner is visible and contains `needle`."""
    page.wait_for_function(
        """needle => {
            const el = document.querySelector('#gene-panel-summary');
            return el && getComputedStyle(el).display !== 'none' && el.innerText.includes(needle);
        }""",
        arg=needle,
        timeout=30000,
    )


def wait_results_contain(page, needle):
    """Wait until the results table body contains `needle`."""
    page.wait_for_function(
        """needle => {
            const el = document.querySelector('#results-table tbody');
            return el && el.innerText.includes(needle);
        }""",
        arg=needle,
        timeout=30000,
    )


def upload(page, path):
    """Set the hidden gene-file input and wait for the chip to render."""
    page.set_input_files("#gene-file-input", str(path))
    page.wait_for_selector("#gene-file-chip-container", state="visible", timeout=10000)


def test_typed_gene_list_all_found(index_page, console_errors):
    """A comma-separated list of real genes shows an all-found summary and results."""
    index_page.fill("#search-query", "SMN1, CFTR, DMD")
    index_page.click("#search-button")
    wait_summary_contains(index_page, "All 3 genes were found")
    assert "All 3 genes were found in the database" in summary_text(index_page)
    wait_results_contain(index_page, "SMN1")
    assert real_errors(console_errors) == [], f"JS errors: {real_errors(console_errors)}"


def test_typed_gene_list_reports_missing(index_page, console_errors):
    """A gene ending in a digit still splits correctly; a fake gene is reported missing."""
    index_page.fill("#search-query", "SMN1, FAKEGENE123, CFTR")
    index_page.click("#search-button")
    wait_summary_contains(index_page, "FAKEGENE123")
    text = summary_text(index_page)
    assert "1 out of 3 genes were not found in the database" in text
    assert "FAKEGENE123" in text
    assert real_errors(console_errors) == [], f"JS errors: {real_errors(console_errors)}"


def test_gene_file_upload(index_page, console_errors, tmp_path):
    """Uploading a TSV gene list restricts results and reports found vs. missing."""
    f = tmp_path / "panel.tsv"
    f.write_text("gene_symbol\tnote\nBRCA1\tx\nCFTR\ty\nNOTAREALGENE\tz\n")
    upload(index_page, f)

    chip = index_page.inner_text("#gene-file-chip-container")
    assert "panel.tsv" in chip and "3 genes" in chip

    index_page.click("#search-button")
    wait_summary_contains(index_page, "NOTAREALGENE")
    assert "1 out of 3 genes were not found in the database" in summary_text(index_page)
    wait_results_contain(index_page, "BRCA1")
    assert real_errors(console_errors) == [], f"JS errors: {real_errors(console_errors)}"


def test_txt_regions_upload(index_page, console_errors, tmp_path):
    """A TXT whose first column has hg38 regions returns overlapping genes and reports coverage."""
    f = tmp_path / "regions.txt"
    # chr17:43044295-43125483 (1-based) overlaps BRCA1; chr1:100-200 overlaps nothing.
    f.write_text("chr17:43044295-43125483\nchr1:100-200\n")
    upload(index_page, f)
    assert "2 regions" in index_page.inner_text("#gene-file-chip-container")

    index_page.click("#search-button")
    wait_summary_contains(index_page, "chr1:100-200")
    assert "1 out of 2 regions had no overlapping gene" in summary_text(index_page)
    wait_results_contain(index_page, "BRCA1")
    assert real_errors(console_errors) == [], f"JS errors: {real_errors(console_errors)}"


def test_drop_file_on_main_search_box(index_page, console_errors, tmp_path):
    """Dropping a file directly on the main search box (not just in the upload modal) behaves the
    same as dropping it into the modal's dropzone."""
    page = index_page
    data_transfer = page.evaluate_handle(
        """() => {
            const dt = new DataTransfer()
            const file = new File(['BRCA1\\nCFTR\\n'], 'dropped_on_search_box.tsv', { type: 'text/tab-separated-values' })
            dt.items.add(file)
            return dt
        }"""
    )
    page.dispatch_event("#gene-main-search-dropzone", "drop", {"dataTransfer": data_transfer})
    page.wait_for_selector("#gene-file-chip-container", state="visible", timeout=10000)
    chip = page.inner_text("#gene-file-chip-container")
    assert "dropped_on_search_box.tsv" in chip and "2 genes" in chip

    page.click("#search-button")
    wait_summary_contains(page, "All 2 genes were found")
    wait_results_contain(page, "BRCA1")
    assert real_errors(console_errors) == [], f"JS errors: {real_errors(console_errors)}"


def test_txt_mixed_genes_and_regions(index_page, console_errors, tmp_path):
    """A TXT mixing a gene and a region is treated as 'items' (gene found, region not)."""
    f = tmp_path / "mixed.txt"
    f.write_text("BRCA1\nchr1:100-200\n")
    upload(index_page, f)
    assert "2 items" in index_page.inner_text("#gene-file-chip-container")

    index_page.click("#search-button")
    wait_summary_contains(index_page, "chr1:100-200")
    text = summary_text(index_page)
    assert "1 out of 2 items were not found" in text
    wait_results_contain(index_page, "BRCA1")
    assert real_errors(console_errors) == [], f"JS errors: {real_errors(console_errors)}"


def test_gene_file_chip_remove(index_page, tmp_path):
    """Removing the chip (x) clears the uploaded gene list."""
    f = tmp_path / "panel.txt"
    f.write_text("BRCA1\nCFTR\n")
    upload(index_page, f)
    assert index_page.is_visible("#gene-file-chip-container")

    index_page.click("#gene-file-chip-remove")
    assert not index_page.is_visible("#gene-file-chip-container")


def test_bed_file_upload(index_page, console_errors, tmp_path):
    """Uploading a BED file returns overlapping genes and reports region coverage."""
    f = tmp_path / "regions.bed"
    # chr17:43,044,295-43,125,483 (1-based) overlaps BRCA1; chr1:100-200 overlaps nothing.
    f.write_text("track name=test\nchr17\t43044294\t43125483\tBRCA1\nchr1\t100\t200\tnothing\n")
    upload(index_page, f)

    assert "2 regions" in index_page.inner_text("#gene-file-chip-container")

    index_page.click("#search-button")
    wait_summary_contains(index_page, "chr1:100-200")
    text = summary_text(index_page)
    assert "1 out of 2 regions had no overlapping gene" in text
    assert "chr1:100-200" in text
    wait_results_contain(index_page, "BRCA1")
    assert real_errors(console_errors) == [], f"JS errors: {real_errors(console_errors)}"


def test_bed_empty_interval_skipped(index_page, tmp_path):
    """A zero-length BED interval (start == end) is skipped, not counted as a region."""
    f = tmp_path / "empty.bed"
    # chr17 interval overlaps BRCA1 (valid); the chr1 line is an empty half-open interval.
    f.write_text("chr17\t43044294\t43125483\tBRCA1\nchr1\t100\t100\tempty\n")
    upload(index_page, f)
    # Only the valid interval is kept; the empty one is dropped.
    assert "1 region" in index_page.inner_text("#gene-file-chip-container")


def test_gzip_gene_file_upload(index_page, tmp_path):
    """A gzip-compressed gene list is decompressed and parsed."""
    f = tmp_path / "panel.txt.gz"
    with gzip.open(f, "wt") as fh:
        fh.write("SMN1\nTTN\nZZZFAKE\n")
    upload(index_page, f)
    assert "3 genes" in index_page.inner_text("#gene-file-chip-container")


def test_missing_genes_more_modal(index_page, tmp_path):
    """More than 20 missing genes truncate inline and list one-per-line in the modal."""
    f = tmp_path / "big.txt"
    f.write_text("\n".join(f"FAKE{i}" for i in range(25)) + "\nBRCA1\n")
    upload(index_page, f)

    index_page.click("#search-button")
    wait_summary_contains(index_page, "more (")

    index_page.click("#gene-panel-more-missing")
    index_page.wait_for_selector("#gene-panel-missing-modal", state="visible", timeout=5000)
    modal_text = index_page.input_value("#gene-panel-missing-modal-text")
    lines = [line for line in modal_text.split("\n") if line.strip()]
    assert len(lines) == 25, f"Expected 25 missing genes, got {len(lines)}"
    assert "FAKE24" in lines


def test_panel_plus_keyword_wording(index_page, console_errors, tmp_path):
    """Panel + keyword uses 'did not match your search' wording, not 'not found in the database'."""
    f = tmp_path / "panel.tsv"
    f.write_text("BRCA1\nCFTR\n")
    upload(index_page, f)
    index_page.fill("#search-query", "cystic fibrosis")
    index_page.click("#search-button")
    # CFTR matches "cystic fibrosis"; BRCA1 does not -> reported as not matching.
    wait_summary_contains(index_page, "did not match your search")
    text = summary_text(index_page)
    assert "did not match your search" in text
    assert "not found in the database" not in text
    assert "BRCA1" in text
    assert real_errors(console_errors) == [], f"JS errors: {real_errors(console_errors)}"


def test_nongene_file_rejected_no_widening(index_page, tmp_path):
    """A file whose first column is free-text phrases is rejected, not accepted as an empty panel."""
    f = tmp_path / "phrases.txt"
    f.write_text("spinal muscular atrophy\nlimb girdle muscular dystrophy\n")
    index_page.click("#gene-file-upload-button")
    index_page.wait_for_selector("#gene-file-upload-modal", state="visible", timeout=5000)
    index_page.set_input_files("#gene-file-input", str(f))
    index_page.wait_for_selector("#gene-file-upload-error", state="visible", timeout=5000)
    # No panel is attached, so a subsequent search would NOT be silently widened to the whole table.
    assert not index_page.is_visible("#gene-file-chip-container")


def test_txt_region_thousands_commas(index_page, console_errors, tmp_path):
    """A region using comma thousands-separators in a TXT is parsed, not truncated/dropped."""
    f = tmp_path / "region.txt"
    f.write_text("chr17:43,044,295-43,125,483\n")
    upload(index_page, f)
    assert "1 region" in index_page.inner_text("#gene-file-chip-container")
    index_page.click("#search-button")
    wait_results_contain(index_page, "BRCA1")
    assert real_errors(console_errors) == [], f"JS errors: {real_errors(console_errors)}"


def test_typed_mixed_gene_and_disease_no_banner(index_page):
    """A typed list mixing a gene and a lowercase disease word is not treated as a gene panel."""
    index_page.fill("#search-query", "SMN1, ataxia")
    index_page.click("#search-button")
    index_page.wait_for_selector("#results-table tbody tr", timeout=30000)
    # 'ataxia' is not gene-shaped, so no found/missing banner should appear (results still load).
    assert not index_page.is_visible("#gene-panel-summary")


def test_typed_all_missing_shows_banner(index_page):
    """A typed panel where every gene is absent still shows the found/missing banner."""
    index_page.fill("#search-query", "FAKEGENE123, FAKEGENE456")
    index_page.click("#search-button")
    wait_summary_contains(index_page, "were not found in the database")
    text = summary_text(index_page)
    assert "2 out of 2 genes were not found in the database" in text
    assert "FAKEGENE123" in text and "FAKEGENE456" in text


def test_typed_versioned_ensembl_matches_results(index_page):
    """A versioned ENSG id in a typed panel matches the results table (not just the banner)."""
    # ENSG00000012048 = BRCA1; the ".15" version suffix must be stripped for the results query too.
    index_page.fill("#search-query", "ENSG00000012048.15, CFTR")
    index_page.click("#search-button")
    wait_results_contain(index_page, "BRCA1")
    wait_summary_contains(index_page, "All 2 genes were found in the database")


def test_chip_remove_clears_banner(index_page, tmp_path):
    """Removing the uploaded-file chip clears the (now-orphaned) found/missing banner."""
    f = tmp_path / "panel.tsv"
    f.write_text("BRCA1\nCFTR\nNOTAREALGENE\n")
    upload(index_page, f)
    index_page.click("#search-button")
    wait_summary_contains(index_page, "NOTAREALGENE")
    index_page.click("#gene-file-chip-remove")
    assert not index_page.is_visible("#gene-panel-summary")


def test_multi_file_upload_unions_genes(index_page, console_errors, tmp_path):
    """Uploading several files searches the UNION of their genes, deduped across files."""
    f1 = tmp_path / "panel1.txt"
    f1.write_text("SMN1\nCFTR\n")
    f2 = tmp_path / "panel2.txt"
    f2.write_text("CFTR\nDMD\n")  # CFTR overlaps f1 -> must be counted once in the union
    index_page.set_input_files("#gene-file-input", [str(f1), str(f2)])
    index_page.wait_for_selector("#gene-file-chip-container", state="visible", timeout=10000)

    chip = index_page.inner_text("#gene-file-chip-container")
    assert "2 files" in chip, f"chip should name 2 files, got: {chip}"
    assert "3 genes" in chip, f"union of SMN1/CFTR/DMD should be 3 genes, got: {chip}"

    index_page.click("#search-button")
    wait_summary_contains(index_page, "All 3 genes were found")
    assert "All 3 genes were found in the database" in summary_text(index_page)
    wait_results_contain(index_page, "DMD")
    assert real_errors(console_errors) == [], f"JS errors: {real_errors(console_errors)}"


def test_multi_file_upload_mixes_genes_and_regions(index_page, console_errors, tmp_path):
    """A gene-list file and a BED file uploaded together union into genes + regions."""
    genes = tmp_path / "genes.txt"
    genes.write_text("BRCA1\n")
    regions = tmp_path / "regions.bed"
    # chr7:117480025-117668665 (1-based) overlaps CFTR.
    regions.write_text("chr7\t117480024\t117668665\tCFTR_region\n")
    index_page.set_input_files("#gene-file-input", [str(genes), str(regions)])
    index_page.wait_for_selector("#gene-file-chip-container", state="visible", timeout=10000)

    chip = index_page.inner_text("#gene-file-chip-container")
    assert "2 files" in chip, f"chip should name 2 files, got: {chip}"
    assert "2 items" in chip, f"1 gene + 1 region should be 2 items, got: {chip}"

    index_page.click("#search-button")
    wait_results_contain(index_page, "BRCA1")
    assert real_errors(console_errors) == [], f"JS errors: {real_errors(console_errors)}"


def test_multi_file_upload_skips_unreadable_file(index_page, tmp_path):
    """An unreadable file among several is skipped and reported; the readable ones still apply."""
    good = tmp_path / "good.txt"
    good.write_text("BRCA1\nCFTR\n")
    bad = tmp_path / "bad.gz"
    bad.write_bytes(b"this is not valid gzip content\n")  # .gz name but not gzip -> read error
    # Open the upload modal first: the skipped-file notice lives inside it (as when a user browses/drops).
    index_page.click("#gene-file-upload-button")
    index_page.wait_for_selector("#gene-file-upload-modal", state="visible", timeout=5000)
    index_page.set_input_files("#gene-file-input", [str(good), str(bad)])

    # The readable file's genes are still loaded (chip names only the file that was read).
    index_page.wait_for_selector("#gene-file-chip-container", state="visible", timeout=10000)
    chip = index_page.inner_text("#gene-file-chip-container")
    assert "good.txt" in chip and "2 genes" in chip, f"chip: {chip}"

    # The skipped file is reported (modal stays open with the notice).
    index_page.wait_for_selector("#gene-file-upload-error", state="visible", timeout=5000)
    err = index_page.inner_text("#gene-file-upload-error")
    assert "Skipped 1 file" in err and "bad.gz" in err, f"error notice: {err}"
