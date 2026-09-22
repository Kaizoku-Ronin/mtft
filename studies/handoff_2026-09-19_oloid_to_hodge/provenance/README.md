# Provenance and preservation

The three earlier delivered study ZIPs were retrieved in their current saved versions and retained unchanged in `original_archives/`. Their expanded files are under `studies/`. Original report prose, numeric outputs, fixtures, licenses, and manifests have not been edited to make this consolidation. In particular, the scope of references to “this bundle” inside a study is that original individual study, not the new outer handoff.

The reviewed distribution `mtft-0.32.0.tar.gz` was copied from the source-audit workspace and checked against the SHA-256 recorded in the studies. The archive, rather than an edited or installed checkout, is included. The source and copied MTFT fixtures retain their original licensing; the graded study includes the MIT license. Supplied PDFs and images retain their own source attribution and rights. The handoff does not impose a new blanket license on third-party materials.

The new navigation documents, session summary, oloid note, integrity checker, and combined manifest were prepared during this packaging request. The intermediate derivation is preserved as a working note, subordinate to the completed reports.

The packager checked file identity and integrity. It did not rerun the three research workloads, the full MTFT suite, or an independent review of every theorem. The completed studies already distinguish exact statements from numerical diagnostics and assumptions. Those classifications are preserved in the claim ledger.

`PACKAGING_REPORT.json` records the upstream identity checks and archive membership comparisons. `INVENTORY.csv` at the outer root lists the delivered content and byte counts. `MANIFEST.json` and `SHA256SUMS.txt` provide the final cryptographic checksums. A manifest cannot contain its own hash; the text checksum list also covers the manifest, excluding only itself.
