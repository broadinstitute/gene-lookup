"""Tests for export/download functionality on the index page.

Searches for a query, then exports to each available format (TSV, BED, JSON) and verifies
the generated file downloads successfully, is non-empty, and (for the unambiguous formats)
that its row count matches the total number of search results.
"""

import gzip
import json
import re
import urllib.request

import pytest

from test_basic_search import real_errors, do_search


# Export format values from the dropdown menu.
EXPORT_FORMATS = ["tsv", "bed", "json"]

# Formats that export directly without first showing the column-options dialog.
FORMATS_WITHOUT_DIALOG = {"bed"}

# A query with a single result keeps the export fast and the row count predictable.
EXPORT_QUERY = "BRCA1"


def get_total_results(page):
    """Parse the total result count from 'Showing X-Y out of Z'."""
    text = page.locator(".total-results").first.inner_text()
    match = re.search(r"out of\s+([\d,]+)", text)
    return int(match.group(1).replace(",", "")) if match else 0


def do_export(page, fmt):
    """Trigger an export and return the list of captured download URLs.

    Monkey-patches anchor click to capture the GCS URL instead of navigating away, so the
    exported file can be fetched directly rather than relying on a real browser download.
    """
    page.evaluate("""() => {
        window.__capturedDownloadUrls = [];
        const origClick = HTMLAnchorElement.prototype.click;
        HTMLAnchorElement.prototype.click = function() {
            if (this.href && this.href.includes('storage')) {
                window.__capturedDownloadUrls.push(this.href);
            } else {
                origClick.call(this);
            }
        };
    }""")

    # Accept the "large export split into N files" confirm() if it ever appears.
    page.on("dialog", lambda dialog: dialog.accept())

    # Selecting a value fires the dropdown's onChange, which starts the export.
    page.evaluate(
        "fmt => $('.export-dropdown').first().dropdown('set selected', fmt)", fmt
    )

    if fmt not in FORMATS_WITHOUT_DIALOG:
        page.wait_for_selector("#export-column-dialog", state="visible", timeout=10000)
        page.click("#export-dialog-download-button")

    page.wait_for_function(
        "() => window.__capturedDownloadUrls && window.__capturedDownloadUrls.length > 0",
        timeout=60000,
    )
    return page.evaluate("() => window.__capturedDownloadUrls")


def fetch_export_content(url):
    """Fetch exported file content from a GCS URL, handling gzip decompression.

    Converts storage.cloud.google.com URLs to storage.googleapis.com for direct public access.
    """
    url = url.replace("storage.cloud.google.com/", "storage.googleapis.com/")
    with urllib.request.urlopen(url) as resp:
        data = resp.read()
    # Detect gzip by its magic bytes (0x1f 0x8b) rather than the URL suffix: signed GCS URLs
    # append query params (e.g. ...bed.gz?X-Goog-Algorithm=...), so a `.gz` suffix check misses.
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    return data.decode("utf-8")


@pytest.mark.parametrize("fmt", EXPORT_FORMATS)
def test_export_format(index_page, console_errors, fmt):
    """Export to a given format and verify the download is non-empty and error-free."""
    page = index_page
    do_search(page, EXPORT_QUERY)
    total_results = get_total_results(page)
    assert total_results > 0, "Expected at least one result to export"

    download_urls = do_export(page, fmt)
    assert len(download_urls) > 0, f"No download URL captured for format '{fmt}'"

    for url in download_urls:
        content = fetch_export_content(url)
        assert content.strip(), f"Downloaded file is empty for format '{fmt}'"
        # The searched gene must appear in every format's export.
        assert "BRCA1" in content, f"Exported '{fmt}' file does not contain the searched gene"

        if fmt == "json":
            assert len(json.loads(content)) == total_results
        elif fmt == "bed":
            # BED has no header row, so every non-empty line is a data row.
            lines = [ln for ln in content.strip().split("\n") if ln.strip()]
            assert len(lines) == total_results
        else:  # tsv has a header row followed by one line per result
            lines = [ln for ln in content.strip().split("\n") if ln.strip()]
            assert len(lines) == total_results + 1

    assert real_errors(console_errors) == [], (
        f"JavaScript errors during export '{fmt}': {real_errors(console_errors)}"
    )
