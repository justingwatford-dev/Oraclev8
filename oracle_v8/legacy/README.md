# Legacy per-storm drivers (provenance only)

Superseded by `oracle_v8/run_storm.py` + `oracle_v8/production_config.py` in V8.7.
Retained for historical reproducibility of pre-unification runs.  **Do not use for
new results** — they carry storm-specific configs that drift from the shared
production config.  Each requires `ORACLE_ALLOW_LEGACY=1` to run.
