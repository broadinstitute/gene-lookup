"""Tests for the per-gene detail page (gene.html).

The gene page reads the gene identifier from the URL hash
(#hgnc_gene_id=HGNC:1100 or #ensembl_gene_id=ENSG...), queries the live
BigQuery backend, and renders the gene's details. These tests cover the basic
load paths that previously had no coverage:
  - loading by HGNC id
  - loading by Ensembl gene id
  - the error shown when no gene identifier is supplied

These tests exercise the live BigQuery backend, mirroring test_basic_search.py.
"""


def real_errors(console_errors):
    """Return console errors excluding blocked-CDN / analytics network failures.

    Google Analytics and third-party CDNs are often unreachable in CI/sandbox
    environments, which surfaces as "Failed to load resource" console errors that
    are unrelated to the page's own logic.
    """
    return [
        e for e in console_errors
        if "failed to load resource" not in e.lower()
        and "google-analytics" not in e.lower()
        and "net::err" not in e.lower()
    ]


def wait_for_title_contains(page, needle):
    """Wait until #locus-page-title (set after the gene loads) contains `needle`."""
    page.wait_for_function(
        """needle => {
            const el = document.querySelector('#locus-page-title');
            return el && el.innerText.includes(needle);
        }""",
        arg=needle,
        timeout=30000,
    )


def test_gene_page_loads_by_hgnc_id(page, console_errors, base_url):
    """gene.html#hgnc_gene_id=HGNC:1100 loads BRCA1 with no error and no JS errors."""
    page.goto(f"{base_url}/gene.html#hgnc_gene_id=HGNC:1100")
    wait_for_title_contains(page, "BRCA1")

    assert not page.is_visible("#error-message"), (
        f"Unexpected error on gene page: {page.inner_text('#error-text')}"
    )
    assert real_errors(console_errors) == [], f"JS errors: {real_errors(console_errors)}"


def test_gene_page_loads_by_ensembl_id(page, console_errors, base_url):
    """gene.html#ensembl_gene_id=ENSG00000139618 loads BRCA2 with no error and no JS errors."""
    page.goto(f"{base_url}/gene.html#ensembl_gene_id=ENSG00000139618")
    wait_for_title_contains(page, "BRCA2")

    assert not page.is_visible("#error-message"), (
        f"Unexpected error on gene page: {page.inner_text('#error-text')}"
    )
    assert real_errors(console_errors) == [], f"JS errors: {real_errors(console_errors)}"


def test_gene_page_missing_identifier_shows_error(page, console_errors, base_url):
    """gene.html with no gene identifier in the URL shows a helpful error message."""
    page.goto(f"{base_url}/gene.html")
    page.wait_for_selector("#error-message", state="visible", timeout=30000)

    assert "Missing gene identifier" in page.inner_text("#error-text")
    assert real_errors(console_errors) == [], f"JS errors: {real_errors(console_errors)}"


def test_gene_page_not_found_shows_error(page, console_errors, base_url):
    """A well-formed but non-existent HGNC id reports 'Gene not found', not a crash."""
    page.goto(f"{base_url}/gene.html#hgnc_gene_id=HGNC:99999999")
    page.wait_for_selector("#error-message", state="visible", timeout=30000)

    assert "not found" in page.inner_text("#error-text").lower()
    # The page intentionally console.error()s the not-found condition; allow that expected log
    # while still catching any OTHER unexpected JS errors.
    unexpected = [e for e in real_errors(console_errors) if "gene not found" not in e.lower()]
    assert unexpected == [], f"Unexpected JS errors: {unexpected}"
