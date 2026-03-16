"""Tests for basic search functionality on the index page.

Each example query type is tested: gene name, HGNC ID, Ensembl gene ID,
genomic region, and disease/phenotype text. Tests confirm that results load
without JS errors. After basic search validation, all columns are enabled and
each is sorted ascending and descending.
"""

import re

import pytest


# Queries with expected minimum result counts
SEARCH_QUERIES = [
    ("BRCA1", 1),
    ("HGNC:1100", 1),
    ("ENSG00000012048", 1),
    ("chr17:43044295-43125483", 1),
    ("cardiomyopathy", 1),
]

# Queries expected to have more than one page of results
MULTI_PAGE_QUERIES = ["cardiomyopathy"]


def do_search(page, query):
    """Fill the search box, click Search, and wait for results."""
    page.fill("#search-query", query)
    page.click("#search-button")
    page.wait_for_selector("#results-table tbody tr", timeout=30000)


def wait_for_results_refresh(page):
    """Wait for a search/sort/page-change to complete."""
    page.wait_for_selector("#loading-container", state="visible", timeout=5000)
    page.wait_for_selector("#results-table tbody tr", timeout=30000)


def assert_no_errors(page, console_errors, context):
    """Assert no error message is shown and no JS errors occurred."""
    assert not page.is_visible("#error-message"), (
        f"Error message displayed ({context}): {page.inner_text('#error-text')}"
    )
    assert console_errors == [], (
        f"JavaScript errors ({context}): {console_errors}"
    )


def enable_all_columns(page):
    """Open column settings dialog, check all boxes, and click Apply."""
    page.click("#results-table-settings-button")
    page.wait_for_selector("#results-table-column-settings-dialog", state="visible", timeout=5000)

    for cb in page.query_selector_all(
        ".results-table-column-settings-checkboxes input[type='checkbox']"
    ):
        if not cb.is_checked():
            cb.check(force=True)

    page.click("#results-table-column-settings-apply-button")
    page.wait_for_selector("#results-table tbody tr", timeout=30000)


def get_sortable_column_names(page):
    """Return list of column data-column attribute values for all sortable headers."""
    return [
        h.get_attribute("data-column")
        for h in page.query_selector_all("#results-table th.results-table-sortable-header")
    ]


def get_total_pages(page):
    """Parse total page count from the pagination 'of N' text."""
    text = page.locator(".bottom-pagination .total-pages").inner_text()
    match = re.search(r"of\s+([\d,]+)", text)
    return int(match.group(1).replace(",", "")) if match else 1


def click_next_page(page):
    """Click the next-page button using bottom pagination and wait for refresh."""
    page.locator(".bottom-pagination .next-page").click(force=True)
    wait_for_results_refresh(page)


def click_last_page(page):
    """Click the last-page button using bottom pagination and wait for refresh."""
    page.locator(".bottom-pagination .last-page").click(force=True)
    wait_for_results_refresh(page)


@pytest.mark.parametrize("query,min_results", SEARCH_QUERIES)
def test_search_returns_results(index_page, console_errors, query, min_results):
    """Search for a query and verify results appear with no JS errors."""
    do_search(index_page, query)

    rows = index_page.query_selector_all("#results-table tbody tr")
    assert len(rows) >= min_results, (
        f"Expected at least {min_results} result(s) for '{query}', got {len(rows)}"
    )
    assert_no_errors(index_page, console_errors, f"search '{query}'")


def test_search_brca1_has_expected_gene(index_page, console_errors):
    """Search BRCA1 and verify the result contains BRCA1 gene info."""
    do_search(index_page, "BRCA1")

    row_text = index_page.inner_text("#results-table tbody")
    assert "BRCA1" in row_text, f"Expected 'BRCA1' in results, got: {row_text[:200]}"
    assert_no_errors(index_page, console_errors, "search 'BRCA1'")


@pytest.mark.parametrize("query,min_results", SEARCH_QUERIES)
def test_all_columns_and_sorting(index_page, console_errors, query, min_results):
    """Enable all columns, verify they load, then sort each column both directions."""
    page = index_page
    do_search(page, query)
    enable_all_columns(page)
    assert_no_errors(page, console_errors, f"all columns for '{query}'")

    columns = get_sortable_column_names(page)
    assert len(columns) > 0, "No sortable columns found"

    for column_name in columns:
        selector = f"th.results-table-sortable-header[data-column='{column_name}']"

        # Click to sort (first click sets the column's default sort direction)
        page.click(selector)
        wait_for_results_refresh(page)
        assert_no_errors(page, console_errors, f"sort '{column_name}' click 1 for '{query}'")

        # Click again to toggle sort direction
        page.click(selector)
        wait_for_results_refresh(page)
        assert_no_errors(page, console_errors, f"sort '{column_name}' click 2 for '{query}'")


@pytest.mark.parametrize("query", MULTI_PAGE_QUERIES)
def test_pagination_all_columns(index_page, console_errors, query):
    """Enable all columns, navigate to pages 2, 3, and the last page without errors."""
    page = index_page
    do_search(page, query)
    enable_all_columns(page)
    assert_no_errors(page, console_errors, f"all columns page 1 for '{query}'")

    # Verify there are multiple pages
    total_pages = get_total_pages(page)
    assert total_pages > 1, (
        f"Expected more than 1 page of results for '{query}', got {total_pages}"
    )

    # Navigate to page 2
    click_next_page(page)
    assert_no_errors(page, console_errors, f"all columns page 2 for '{query}'")

    # Navigate to page 3 (if it exists)
    if total_pages >= 3:
        click_next_page(page)
        assert_no_errors(page, console_errors, f"all columns page 3 for '{query}'")

    # Jump to last page (if not already there)
    if total_pages > 3:
        click_last_page(page)
        assert_no_errors(page, console_errors, f"all columns last page ({total_pages}) for '{query}'")
