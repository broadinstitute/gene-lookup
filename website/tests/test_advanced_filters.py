"""Tests for the advanced-filter panel on the index page.

Covers applying a standard advanced filter (source/inheritance) and confirms the
filter is actually reflected in the SQL sent to the backend. Includes a regression
test that a set filter keeps applying after the advanced-filter panel is collapsed
(the panel is a visibility toggle; collapsing it must not silently drop filters on a
subsequent sort/pagination search).

These tests inspect the POST body sent to the query proxy rather than asserting on
result counts, so they are deterministic and do not depend on the live data volume.
They still require the live backend to be reachable for the initial search.
"""

import json


def is_proxy_post(request):
    """True for a POST to the BigQuery query proxy."""
    return "query_gene_lookup_db" in request.url and request.method == "POST"


def request_sql(request):
    """Extract the `sql` field from a proxy POST request body."""
    return json.loads(request.post_data)["sql"]


def real_errors(console_errors):
    """Return console errors excluding blocked-CDN / analytics network failures."""
    return [
        e for e in console_errors
        if "failed to load resource" not in e.lower()
        and "google-analytics" not in e.lower()
        and "net::err" not in e.lower()
    ]


def do_search(page, query):
    """Fill the search box, click Search, and wait for results."""
    page.fill("#search-query", query)
    page.click("#search-button")
    page.wait_for_selector("#results-table tbody tr", timeout=30000)


def test_source_filter_adds_clause(index_page, console_errors):
    """Checking the OMIM source filter adds an OMIM clause to the query SQL."""
    page = index_page
    do_search(page, "cardiomyopathy")

    # Open the panel if it's not already visible
    if not page.is_visible("#advanced-search-form"):
        page.click("#filter-button")
    page.locator('input[name="keepOMIM"]').check(force=True)

    with page.expect_request(is_proxy_post) as req_info:
        page.click("#search-button")
    sql = request_sql(req_info.value)

    assert "LIKE '%OMIM%'" in sql, f"OMIM source filter not reflected in query SQL: {sql}"
    page.wait_for_selector("#results-table tbody tr", timeout=30000)
    assert real_errors(console_errors) == [], f"JS errors: {real_errors(console_errors)}"


def test_source_filter_persists_after_panel_collapse(index_page, console_errors):
    """A set source filter still applies to a later (sort-triggered) search after the panel is collapsed.

    Regression test: collapsing the advanced-filter panel is only a visibility toggle and must not
    silently drop still-checked standard filters on the next search (previously it did, while custom
    filters kept applying).
    """
    page = index_page
    do_search(page, "cardiomyopathy")

    # Open the panel if it's not already visible
    if not page.is_visible("#advanced-search-form"):
        page.click("#filter-button")
    page.locator('input[name="keepOMIM"]').check(force=True)
    with page.expect_request(is_proxy_post):
        page.click("#search-button")
    page.wait_for_selector("#results-table tbody tr", timeout=30000)

    # Collapse the panel WITHOUT unchecking the filter, then trigger a re-search by sorting.
    if page.is_visible("#advanced-search-form"):
        page.click("#filter-button")
    with page.expect_request(is_proxy_post) as req_info:
        page.locator("th.results-table-sortable-header").first.click()
    sql = request_sql(req_info.value)

    assert "LIKE '%OMIM%'" in sql, (
        f"source filter was dropped after collapsing the advanced-filter panel; SQL: {sql}"
    )
    assert real_errors(console_errors) == [], f"JS errors: {real_errors(console_errors)}"


def test_search_with_no_results_shows_message(index_page, console_errors):
    """A query that matches nothing shows the 'No matching results found' message, not a crash."""
    page = index_page
    page.fill("#search-query", "xyzzynotarealdiseaseterm")
    page.click("#search-button")
    page.wait_for_selector(".no-results", state="visible", timeout=30000)

    assert "No matching results found" in page.inner_text(".no-results")
    assert real_errors(console_errors) == [], f"JS errors: {real_errors(console_errors)}"
