Repository: BlazeQC — Mojo-based FASTQ QC tool

Short goal for the assistant
- Help implement, refactor, and test features in a small Mojo codebase that performs FASTQ parsing and QC reporting (HTML + plots via Python matplotlib).

Key architecture and components (what to read first)
- Entry points: `QC.mojo` and `blazeseq/record.mojo` — `QC.mojo` contains a small `main` used for quick runs; `record.mojo` and `blazeqc/stats.mojo` contain the core data flow for parsing and aggregating stats.
- Parser & data model: `blazeseq/kseq.mojo` and `blazeseq/record.mojo` define `FastqRecord`, low-level buffer types and record iterators.
- Analysis & output: `blazeqc/stats.mojo` implements analyzers (BasepairDistribution, QualityDistribution, AdapterContent, etc.) and `blazeqc/html_maker.mojo` builds the HTML report and panels.
- Helpers & constants: `blazeqc/helpers.mojo` has utility functions that mix Mojo and Python interop (image encoding, numpy conversions), and `blazeqc/CONSTS.mojo` stores schemas and config-like constants.

Development workflows and important commands
- This project uses Mojo source files and relies on a local Python environment for plotting. Typical quick run (from repo root): run the compiled Mojo binary or execute the `QC.mojo`/`record.mojo` main when available in your Mojo environment. (No build script in repo — prefer your local Mojo toolchain to compile/run.)
- Python dependency path: many modules call Python with an explicit path: `.pixi/envs/default/lib/python3.12/site-packages/`. If running outside the original environment, update `py_lib` aliases in `blazeqc/*` modules to point to your Python site-packages (or add that path to PYTHONPATH at runtime).
- To produce HTML output locally: call the `main` in `record.mojo` with a FASTQ filename. The program tallies records and writes `<filename>_blazeseq.html`.

Project-specific conventions & patterns
- Mixed Mojo↔Python interop: plotting and some heavy numeric ops are done by importing Python modules at runtime (matplotlib, numpy, seaborn). Mojo code converts tensors to numpy arrays with helpers in `blazeqc/helpers.mojo` — update these helpers first if you change plotting logic.
- Memory & performance idioms: the code uses manual UnsafePointer and Tensor grow/copy helpers (`grow_tensor`, `grow_matrix`, `cpy_tensor`). When changing data structures prefer reusing the same patterns to avoid introducing hidden allocation bugs.
- Hand-rolled containers: Mojo Dict internals are used (e.g., `_find_index`, direct `_entries` access in `PerTileQuality`) — be careful when refactoring; these low-level accesses are deliberate performance optimizations.
- Encoding/hash conventions: adapters and kmers use compact integer hashes (`_seq_to_hash`, bit-packed kmer hashing). Keep bit-width assumptions in mind (see `bits` generic param and limits like 64//bits).

Things to watch / gotchas (examples from the code)
- Hardcoded Python site-packages path: many files set `py_lib` to `.pixi/envs/default/...`. Make that configurable before running outside the dev environment.
- Assumed quality schema offsets: `QualitySchema` offsets (33 vs 64) are used widely (see `blazeseq/quality_schama.mojo` and `blazeqc/stats.mojo` for `_guess_schema`). Changing schema handling affects many plots.
- Low-level string/buffer handling: `blazeseq/kseq.mojo` uses UnsafePointer and manual find/resize logic. Introduce higher-level helpers only after ensuring tests cover edge cases (newlines, CRLF, empty records).
- Large buffers to Python boundary: converting large tensors to numpy with elementwise loops is currently used and may be slow for huge datasets — prefer batching or memory views if performance becomes an issue.

How to change code safely
- When adding features that change record shape or statistics, update the corresponding `Analyser` implementations in `blazeqc/stats.mojo` and the HTML insertion points in `blazeqc/html_maker.mojo` to produce `result_panel` entries.
- For plotting changes, edit `blazeqc/helpers.mojo::encode_img_b64` and the `py_lib` path first to ensure Python modules are importable.
- Preserve existing public APIs: `FastqRecord`, `FullStats`, and `result_panel` are used across modules — keep their method signatures stable or add backwards-compatible overloads.

Files to reference when implementing tasks (examples)
- FASTQ parsing and record model: `blazeseq/kseq.mojo`, `blazeseq/record.mojo`
- Aggregation and report generation: `blazeqc/stats.mojo`, `blazeqc/html_maker.mojo`
- Helpers and constants: `blazeqc/helpers.mojo`, `blazeqc/CONSTS.mojo`, `blazeqc/config.mojo`

If something is missing
- Ask for the preferred Mojo toolchain and Python environment paths to produce runnable examples. I can also update `py_lib` usage to read from an env var or repository config if you want.

Please review and tell me which sections you want expanded (examples, run commands, or more file-level notes) and I will iterate.